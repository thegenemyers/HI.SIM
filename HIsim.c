/*******************************************************************************************
 *
 *  Synthetic DNA shotgun dataset simulator
 *     From a supplied reference genome in the form of a Dazzler .dam, sample reads of
 *     mean length -m from a log-normal length distribution with standard deviation -s,
 *     but ignore reads of length less than -x.  Collect enough reads to cover the genome
 *     -c times.   Introduce -e fraction errors into each read where the ratio of insertions,
 *     deletions, and substitutions are set by defined constants INS_RATE and DEL_RATE
 *     within generate.c.  The fraction -f controls the rate at which reads are picked from
 *     the forward and reverse strands which defaults to 50%.  If -C is set then assume the
 *     scaffolds are circular.
 *
 *     The -r parameter seeds the random number generator for the generation of the genome
 *     so that one can reproducbile produce the same underlying genome to sample from.  If
 *     missing, then the job id of the invocation seeds the generator.  The output is sent
 *     to the standard output (i.e. it is a pipe).  The output is in fasta format (i.e. it is
 *     a UNIX pipe).  The output is in Pacbio .fasta format suitable as input to fasta2DB.
 *
 *     The genome is considered a sequence of *scaffolds* (these are reconstituted from the
 *     Dazzler's internal encoding of a .dam), where the gaps are filled with a random
 *     sequence that follows the base distribution of the contigs of the genome.  The program
 *     then samples these filled in scaffolds for reads.  If the -C optioin is set then the
 *     program assumes each scaffold is a circular sequence.
 *
 *     The -M option requests that the scaffold and coordinates from which each read has
 *     been sampled are written to the indicated file, one line per read, ASCII encoded.
 *     This "map" file essentially tells one where every read belongs in an assembly and
 *     is very useful for debugging and testing purposes.  If a read pair is say b,e then
 *     if b < e the read was sampled from [b,e] in the forward direction, and from [e,b]
 *     in the reverse direction otherwise.
 *
 *  Author:  Gene Myers
 *  Date  :  June 2021
 *
 ********************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <ctype.h>

#undef   DEBUG_HAPLO
#undef   DEBUG_HAPLO_SCAN
#undef   DEBUG_HAPLO_GEN
#undef   DEBUG_HAPLO_SNPS
#undef   DEBUG_HAPLO_MAKE
#undef   SHOW_HAPLO
#undef   DEBUG_SHOTGUN
#undef   DEBUG_TABLE
#undef   DEBUG_MUTATE
#undef   DEBUG_OPS

#include "gene_core.h"

static char *Usage[] =
  { " <genome>[.fast[aq]] <model>[.model] [-p<tree>] [-o<name>] [-ehf]",
    " [-c<double(50.)>] [-m<int>] [-s<int>] [-x<int(0)>]",
    " [-vCU] [-w<int(100)>] [-r<int>]"
  };

static int    VERBOSE;    // -v option
static int    CIRCULAR;   // -C option
static int    UPPER;      // -U option
static int    RMEAN;      // -m option
static int    RSDEV;      // -s option
static double COVERAGE;   // -c option
static int    RSHORT;     // -x option
static int    WIDTH;      // -w option
static uint32 SEED;       // -r option
static char  *OUT;        // -o option
static int    ERRINFO;    // -e option
static int    HAPINFO;    // -h option
static int    FASTOUT;    // -f option
static FILE  *READ_OUT;   // OUT.fasta if -f, stdout otherwise
static FILE  *ERR_OUT;    // OUT.err   if -e

static uint64 GENER;      // Generator state

#define PARETO_MIN     1.0     //  Pareto distribution parameters
#define PARETO_SHAPE   1.025

#define SV_DELETION  .50   //  Proportion of each type of SV
#define SV_RANDOM    .75
#define SV_TANDEM   1.00

#define FLIP_RATE    .5  //  Forward/reverse sampling ratio

static char dna[4] = { 'a', 'c', 'g', 't' };


/*******************************************************************************************
 *
 *  Distribution Samplers
 *
 ********************************************************************************************/


//  Psuedo-randome number generator

static uint64 myrand48_a = 0x5deece66dull;
static uint64 myrand48_c = 0xbull;

static inline double erand()
{ uint64 temp;

  GENER = (GENER * myrand48_a + myrand48_c) & 0xffffffffffffull;
  temp = 0x3ff0000000000000ull | (GENER << 4);
  return (*((double *) &temp) - 1.0);
}

static inline void eseed(uint32 seedval)
{ GENER = ((((uint64) seedval) << 16) | 0x330e); }


//  Bin search distribution table

static int bin_search(int len, double *tab, double y)
{ int l, m, r;

  // Searches tab[0..len] for min { r : y < tab[r] }.
  //   Assumes y < 1 and tab[len] = 1.
  //   So returned index is in [0,len].

  l = 0;
  r = len;
  while (l < r)
    { m = (l+r) >> 1;
      if (y < tab[m])
        r = m;
      else
        l = m+1;
    }
  return (r);
}

static inline int sample_exponential(double mean)
{ return ((int) (-log(1.-erand())*mean)); }

static inline double sample_pareto(double min, double shape)
{ return ((int) (min / powl(1-erand(),1./shape))); }



/*******************************************************************************************
 *
 *  Read fasta or fastq sequences
 *
 ********************************************************************************************/

typedef struct
  { uint64  hlen;    // length of header line
    char   *header;  // header line
    uint64  slen;    // length of DNA sequence
    char   *seq;     // DNA sequence
  } Entry;

static Entry *(*Fetch)(FILE *, Entry *);  //  Either  points at Get_Fasta_Entry or Get_Fastq_Entry
                                          //    depending on assembly file type.

//  If here == NULL scan next line using fixed-size buffers, return ptr to 1st 1000 bps of
//    line and set *plen to the length of the line scanned.
//  If here != NULL load next line into here, return a ptr to it and set *plen to the length
//    of the line loaded.

static char *Read_Line(int64 *plen, char *here, FILE *in)
{ static char buffer1[1000];
  static char buffer2[1000];
  char *where;
  int   len, rln;
 
  if (here == NULL)
    where = buffer1;
  else
    where = here;

  if (fgets(where,1000,in) == NULL)
    return (NULL);
 
  len = rln = strlen(where);
  while (where[rln-1] != '\n')
    { if (here == NULL)
        where = buffer2;
      else
        where += rln;
      if (fgets(where,1000,in) == NULL)
        { fprintf(stderr,"%s: Last line of file does not end with new-line\n",Prog_Name);
          exit (1);
        }
      rln = strlen(where);
      len += rln;
    }
  where[--rln] = '\0';
  len -= 1;

  *plen = len;
  if (here == NULL)
    return (buffer1);
  else
    return (here);
}

//  Read next fastq entry returning it in supplied Entry record.   If on input either
//    ->seq or ->header are NULL, then do not load the data, just set the lengths,
//    otherwise data is loaded starting at the given pointers.

static Entry *Get_Fastq_Entry(FILE *in, Entry *entry)
{ int64   hlen, slen, qlen;
  char   *line;

  line = Read_Line(&hlen,entry->header,in);
  if (line == NULL)
     return (NULL);

  if (line[0] != '@')
    { fprintf(stderr,"%s: Entry header does not start with an @-sign\n",Prog_Name);
      exit (1);
    }

  line = Read_Line(&slen,entry->seq,in);
  if (line == NULL)
    { fprintf(stderr,"%s: Incomplete fastq entry\n",Prog_Name);
      exit (1);
    }

  line = Read_Line(&qlen,NULL,in);
  if (line == NULL)
    { fprintf(stderr,"%s: Incomplete fastq entry\n",Prog_Name);
      exit (1);
    }
  if (line[0] != '+')
    { fprintf(stderr,"%s: Divider line does not start with a +-sign\n",Prog_Name);
      exit (1);
    }

  line = Read_Line(&qlen,NULL,in);
  if (line == NULL)
    { fprintf(stderr,"%s: Incomplete fastq entry\n",Prog_Name);
      exit (1);
    }
  if (slen != qlen)
    { fprintf(stderr,"%s: QV line does not have the same length as sequence line\n",Prog_Name);
      exit (1);
    }

  entry->slen = slen;
  entry->hlen = hlen;
  return (entry);
}

//  Read next fasta entry returning it in supplied Entry record.   If on input either
//    ->seq or ->header are NULL, then do not load the data, just set the lengths,
//    otherwise data is loaded starting at the given pointers.

static Entry *Get_Fasta_Entry(FILE *in, Entry *entry)
{ int64 hlen, slen, m;
  char *line;
  int   c;

  line = Read_Line(&hlen,entry->header,in);
  if (line == NULL)
    return (NULL);

  if (line[0] != '>')
    { fprintf(stderr,"%s: First line does not start with a >-sign\n",Prog_Name);
      exit (1);
    }

  line = Read_Line(&slen,entry->seq,in);
  if (line == NULL)
    { fprintf(stderr,"%s: Sequence missing after last header\n",Prog_Name);
      exit (1);
    }
  if (line[0] == '>')
    { fprintf(stderr,"%s: Sequence missing between two headers\n",Prog_Name);
      exit (1);
    }

  c = fgetc(in);
  while (c != EOF)
    { ungetc(c,in);
      if (c == '>')
        break;

      if (entry->seq == NULL)
        line = Read_Line(&m,NULL,in);
      else
        line = Read_Line(&m,entry->seq+slen,in);
      slen += m;
      c = fgetc(in);
    }

  entry->hlen = hlen;
  entry->slen = slen;
  return (entry);
}


/*******************************************************************************************
 *
 *  Read in error model produced by FKmodel
 *
 ********************************************************************************************/

typedef struct
  { float  all;
    float  ins;
    float  op[9];
  } E_Rates;

typedef struct
  { float all;
    float op[6];
  } M_Rates;

static int KMAX[4];

static E_Rates  HepTab[0x4000];
static M_Rates *Mic3Tab[0x40];
static M_Rates *Mic2Tab[0x10];
static M_Rates *Mic1Tab[0x04];

static char  Bell_Adapter[100] = "ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT";
static int   Bell_Len;
static char *Bell_End;

static int Hex[0x1000];
static int HxB[0x1000];

#if defined(DEBUG_OPS) || defined(DEBUG_MUTATE)

static int MICINS[4] = { '.', 'H', 'Z', 'T' };
static int MICDEL[4] = { '.', 'h', 'z', 't' };

#endif

static int Number[256] =
    { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1,  0, -1,  1, -1, -1, -1,  2, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1,  3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1,  0, -1,  1, -1, -1, -1,  2, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1,  3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
    };

static int64  RMean;
static int64  RsDev;
static double Rscale;
static double Roffset;

static int Rlimit;
static int Elimit;
static int Rbuck;
static int Ebuck;
static int AvErr;

static double  *Rlength;
static double **Rerror;

static double sample_read_length(double *erate)
{ double x, *e;
  int    f, rlen;

  x = erand();
  f = bin_search(Rlimit,Rlength,x);    // Bin. search

#ifdef DEBUG
  printf("Search %g -> %d",x,f);
#endif

  // Linear interpolate between table points

  if (f == 0)
    rlen = ( x / Rlength[f] ) * Rbuck;
  else
    rlen = ( f + (x-Rlength[f-1]) / (Rlength[f] - Rlength[f-1]) ) * Rbuck;

  x = erand();
  e = Rerror[rlen/Ebuck];
  f = bin_search(Elimit,e,x);

#ifdef DEBUG
  printf("Search %g -> %d",x,f);
#endif

  if (f == 0)
    *erate = x / e[f];
  else
    *erate = f + (x-e[f-1]) / (e[f] - e[f-1]);

  return (Roffset + Rscale*rlen);
}

static void Load_Error_Model(char *name)
{ FILE *efile;
  int   krange;

  if (VERBOSE)
    fprintf(stderr,"\n  Loading error model %s\n",name);
  fflush(stderr);

  { char *pwd, *root;

    pwd   = PathTo(name);
    root  = Root(name,".model");
    efile = fopen(Catenate(pwd,"/",root,".model"),"r");
    if (efile == NULL)
      { fprintf(stderr,"%s: Cannot open FK error model %s with extension .model\n",
                       Prog_Name,name);
        exit (1);
      }
    free(root);
    free(pwd);
  }

  { M_Rates *mic1, *mic2, *mic3;
    int      i, kmer;

    if (fread(&kmer,sizeof(int),1,efile) != 1)
      goto error_exit;

    krange = kmer/2 - 6;

    mic1 = ((M_Rates *) Malloc(sizeof(M_Rates)* 4*krange,"Homo-Table"));
    mic2 = ((M_Rates *) Malloc(sizeof(M_Rates)*16*krange,"Di-Table"));
    mic3 = ((M_Rates *) Malloc(sizeof(M_Rates)*64*krange,"Tri-Table"));

    if ((int) fread(HepTab,sizeof(E_Rates),0x4000,efile) != 0x4000)
      goto error_exit;
    if ((int) fread(mic1,sizeof(M_Rates),krange*4,efile) != krange*4)
      goto error_exit;
    if ((int) fread(mic2,sizeof(M_Rates),krange*16,efile) != krange*16)
      goto error_exit;
    if ((int) fread(mic3,sizeof(M_Rates),krange*64,efile) != krange*64)
      goto error_exit;

    mic1 -= 2;
    mic2 -= 4;
    mic3 -= 6;

    KMAX[1] = 2+krange;
    KMAX[2] = 4+krange;
    KMAX[3] = 6+krange;

    for (i = 0; i < 0x04; i++)
      Mic1Tab[i] = mic1 + krange*i;

    for (i = 0; i < 0x10; i++)
      Mic2Tab[i] = mic2 + krange*i;

    for (i = 0; i < 0x40; i++)
      Mic3Tab[i] = mic3 + krange*i;
  }

  { int i, x;

    bzero(Hex,sizeof(int)*0x1000);

    for (x = 0; x < 4096; x += 16)
      { Hex[x|0x0] = Hex[x|0x5] = Hex[x|0xa] = Hex[x|0xf] = 1;

        HxB[x|0x0] = 0;
        HxB[x|0x5] = 1;
        HxB[x|0xa] = 2;
        HxB[x|0xf] = 3;
      }

    for (x = 0; x < 4096; x += 256)
      { Hex[x|0x11] = Hex[x|0x22] = Hex[x|0x33] = Hex[x|0x44] = 2;
        Hex[x|0x66] = Hex[x|0x77] = Hex[x|0x88] = Hex[x|0x99] = 2;
        Hex[x|0xbb] = Hex[x|0xcc] = Hex[x|0xdd] = Hex[x|0xee] = 2;

        HxB[x|0x11] = HxB[x|0x22] = HxB[x|0x33] = 0;
        HxB[x|0x44] = HxB[x|0x66] = HxB[x|0x77] = 1;
        HxB[x|0x88] = HxB[x|0x99] = HxB[x|0xbb] = 2;
        HxB[x|0xcc] = HxB[x|0xdd] = HxB[x|0xee] = 3;
      }

    for (x = 1; x < 63; x++)
      if (x != 0x15 && x != 0x2a)
        { Hex[x << 6 | x] = 3;
          HxB[x << 6 | x] = (x >> 4);
        }

    for (i = 0; Bell_Adapter[i] != 0; i++)
      Bell_Adapter[i] = Number[(int) Bell_Adapter[i]];
    Bell_Len = i;
    Bell_End = Bell_Adapter+i;
  }

#ifdef DEBUG_TABLE
  { int      i, j, k, n;
    int      mic_size[4] = { 1, 4, 16, 64 };
    M_Rates *mtab;

    for (i = 0; i < 0x4000; i++)
      { printf("\n");
        for (j = 12; j >= 0; j -= 2)
          printf("%c",dna[(i>>j)&0x3]);
        printf("\n");
        for (j = 0; j < 4; j++)
          printf("  I %c = %.3f\n",dna[j],100.*HepTab[i].op[j]);
        for (j = 0; j < 4; j++)
          printf("  S %c = %.3f\n",dna[j],100.*HepTab[i].op[j+4]);
        printf("  D   = %.3f\n",100.*HepTab[i].op[8]);
        printf("  +   = %.3f\n",100.*HepTab[i].all);
        printf("  +I  = %.3f\n",100.*HepTab[i].ins);
      }

    for (n = 1; n <= 3; n++)
      for (i = 0; i < mic_size[n]; i++)
        { if (n == 1)
            mtab = Mic1Tab[i];
          else if (n == 2)
            mtab = Mic2Tab[i];
          else
            mtab = Mic3Tab[i];
          printf("\n%d ",n);
          for (k = 2*(n-1); k >= 0; k -= 2)
            printf("%c",dna[(i>>k)&0x3]);
          printf("\n");
          for (k = 2*n; k < krange+2*n; k++)
            { printf("  D %2d:",k);
              for (j = 0; j < n; j++)
                printf(" %.3f",100.*mtab[k].op[j]);
              printf("\n");
            }
          for (k = 2*n; k < krange+2*n; k++)
            { printf("  I %2d:",k);
              for (j = 0; j < n; j++)
                printf(" %.3f",100.*mtab[k].op[j+n]);
              printf("\n");
            }
        }

    fflush(stdout);
  }
#endif

  { int i, emax;

    fread(&RMean,sizeof(int64),1,efile);
    fread(&RsDev,sizeof(int64),1,efile);

    fread(&Rlimit,sizeof(int),1,efile);
    fread(&Rbuck,sizeof(int),1,efile);
    fread(&Ebuck,sizeof(int),1,efile);
    fread(&Elimit,sizeof(int),1,efile);
    fread(&AvErr,sizeof(int),1,efile);

    emax = ((Rlimit+1)*Rbuck-1)/Ebuck;

    Rlength = Malloc(sizeof(double)*(Rlimit+1),"Read length distribution");
    Rerror  = Malloc(sizeof(double *)*(emax+1),"Read error distribution");
    if (Rlength == NULL || Rerror == NULL)
      exit (1);

    Rerror[0] = Malloc(sizeof(double)*(emax+1)*(Elimit+1),"Read error distribution");
    for (i = 1; i <= emax; i++)
      Rerror[i] = Rerror[i-1] + (Elimit+1);

    fread(Rlength,sizeof(double),Rlimit+1,efile);
    fread(Rerror[0],sizeof(double),(emax+1)*(Elimit+1),efile);
  }

#ifdef DEBUG_TABLE
  { int    i, e, emax;
    double last;

    printf("\nHistogram of Real Lengths: Mean = %lld SDev = %lld Averr = %.2f%%\n",
           RMean,RsDev,AvErr/100.);

    last = 0.;
    for (i = 0; i <= Rlimit; i++)
      { if (Rlength[i] != last)
          printf(" %5d: %.5f\n",i*Rbuck,Rlength[i]);
        last = Rlength[i];
      }
    printf("\n");

    emax = ((Rlimit+1)*Rbuck-1)/Ebuck;

    printf("\nHistograms of Error Rates\n");
    for (i = 0; i <= emax; i++)
      { printf("%5d - %5d:\n",i*Ebuck,(i+1)*Ebuck);
        last = 0.;
        for (e = 0; e <= Elimit; e++)
          { if (Rerror[i][e] != last)
              printf("  %3d: %.5f\n",e,Rerror[i][e]);
            last = Rerror[i][e];
         }
      }

    fflush(stdout);
  }
#endif

  fclose(efile);

  return;

error_exit:
  fprintf(stderr,"%s: Read of .model file %s failed\n",Prog_Name,name);
  exit (1);
}


/*******************************************************************************************
 *
 *  Mutate read
 *     Produce a sequence of mutations in ops and output a CIGAR string reflecting it.
 *     Return the difference in the lengths of the mutated vs. source sequence.
 *
 ********************************************************************************************/

#if defined(DEBUG_OPS) || defined(DEBUG_MUTATE)

static char *op_text[9] = { "I a", "I c", "I g", "I t", "S a", "S c", "S g", "S t", "D" };

#endif

typedef struct
  { int    pos;
    uint8  kind;
    uint8  data;
    int16  del;
  } Edit;

#define NORMAL 0
#define MICRO  1

static int mutate_read(uint8 *seq, int len, Edit *ops, int emax, double erate)
{ static int delta[] = { 1, 1, 1, 1, 0, 0, 0, 0, -1 };
  static int dbase[] = { 0, 1, 2, 3, 0, 1, 2, 3,  0 }; 

  int      e1, e2, e3;
  E_Rates *h;
  M_Rates *g;
  double   u, v, w, t;
  uint32   x, y, z;
  int      i, m, k;
  int      j, q, n;
  int      ma, loc, etop;
  double   ex;

  e1 = seq[len];
  e2 = seq[len+1];
  e3 = seq[len+2];
  seq[len] = Bell_Adapter[0];
  seq[len+1] = Bell_Adapter[1];
  seq[len+2] = Bell_Adapter[2];

  x = 0;
  for (i = 0; i < 5; i++)
    x = (x << 2) | Bell_End[i-5];
  x = (x << 2) | seq[0];
  y = ((((x & 0x3ff) << 2) | seq[1]) << 2) | seq[2];

  u = erand();
#ifdef DEBUG_MUTATE
  printf("\nran %g\n",u);
#endif

#ifdef DEBUG_MUTATE
  printf("%5d: %c\n",i,dna[seq[0]]);
  printf("%5d: %c\n",i,dna[seq[1]]);
  printf("%5d: %c\n",i,dna[seq[2]]);
#endif

  etop = 0;
  loc = 0;
  m   = Hex[x];
  k   = -2*m;
  for (i = 0; i < len; i++)
    { z = ((x << 2) | seq[i+1]) & 0xfff;
      y = ((y << 2) | seq[i+3]) & 0x3fff;

      h = HepTab+y;
      v = h->all;
      t = u/erate;
#ifdef DEBUG_MUTATE
      printf("%5d: %c %04x %03x %03x %1d :: %.7f %.7f %.7f",i,dna[seq[i+3]],y,x,z,m,u,t,v);
#endif
      if (m > 0 && Hex[z] != m)
        { k = i-k;
          if (k > KMAX[m])
            k = KMAX[m];
          if (m == 1)
            g = Mic1Tab[x&0x3] + k;
          else if (m == 2)
            g = Mic2Tab[x&0xf] + k;
          else
            g = Mic3Tab[x&0x3f] + k;
          w  = g->all;
          ma = HxB[x];
          v -= h->op[8] + h->op[ma];
#ifdef DEBUG_MUTATE
          printf(" %.7f :: %d.%d",w,m,k);
#endif

          if (erate*(v+w) > 1)
            t = u*(v+w);
          if (t < w)
            { j = 0;
              q = 2*m;
              n = 0;
              while (1)
                { for (; j < q; j++)
                    { v = g->op[j];
                      if (t < v)
                        { if (j < m)
                            n += j+1;
                          else
                            n -= (j-m)+1;
                          break;
                        }
                      else
                        t -= v;
                    }
                  u = erand();
#ifdef DEBUG_MUTATE
                  printf(" [%d/%d] <%g,%g>",j,n,u,g->op[0]);
#endif
                  if (j == m-1)
                    { k += m;
                      if (k <= KMAX[m])
                        g += m;
                      if (u >= g->all)
                        { j = 0;
                          q = m;
                          t = u;
                          continue;
                        }
                    }
                  else if (j == 2*m-1 && j < q)
                    { k -= m;
                      if (k >= 2*m)
                        g -= m;
                      if (u >= g->all)
                        { j = m;
                          q = 2*m;
                          t = u;
                          continue;
                        }
                    }
                  if (n > 0)
                    { if (etop >= emax)
                        return (-1);
                      ops[etop].pos  = i;
                      ops[etop].kind = MICRO;
                      ops[etop].del  = n;
                      ops[etop].data = m;
                      etop += 1;
#ifdef DEBUG_OPS
                      printf(" %5d %c %d\n",i,MICINS[m],n);
#endif
#ifdef DEBUG_MUTATE
                      printf(" [%5d %c %d]",i,MICINS[m],n);
#endif
                    }
                  else
                    { if (etop >= emax)
                        return (-1);
                      ops[etop].pos  = i;
                      ops[etop].kind = MICRO;
                      ops[etop].del  = n;
                      ops[etop].data = m;
                      etop += 1;
#ifdef DEBUG_OPS
                      printf(" %5d %c %d\n",i,MICDEL[m],n);
#endif
#ifdef DEBUG_MUTATE
                      printf(" [%5d %c %d]",i,MICDEL[m],n);
#endif
                    }
                  loc = i+1;
                  break;
                }
            }
          else
            { t -= w;
              if (t < v)
                { ex = h->op[ma];
                  h->op[ma] = 0.;
                repeat1:
                  for (j = 0; j < 8; j++)
                    { v = h->op[j];
                      if (t < v)
                        { if (etop >= emax)
                            return (-1);
                          ops[etop].pos  = i;
                          ops[etop].kind = NORMAL;
                          ops[etop].del  = delta[j];
                          ops[etop].data = dbase[j];
                          etop += 1;
                          loc = i+1;
#ifdef DEBUG_OPS
                          printf(" %5d %s *\n",i,op_text[j]);
#endif
#ifdef DEBUG_MUTATE
                          printf(" <%5d %s>",i,op_text[j]);
#endif
                          break;
                        }
                      else
                        t -= v;
                    }
                  h->op[ma] = ex;
                  u = erand();
                  if (j < 4 && u < h->ins*erate)
                    { t = u;
                      goto repeat1;
                    }
                }
              else
                { v += w;
                  v *= erate;
                  u = (u-v)/(1.-v);
                }
            }
        }
      else if (t < v)
        { if (erate*v > 1)
            t = u*v;
        repeat2:
          for (j = 0; j < 9; j++)
            { v = h->op[j];
              if (t < v)
                { if (etop >= emax)
                    return (-1);
                  ops[etop].pos  = i;
                  ops[etop].kind = NORMAL;
                  ops[etop].del  = delta[j];
                  ops[etop].data = dbase[j];
                  etop += 1;
                  loc = i+1;
#ifdef DEBUG_OPS
                  printf(" %5d %s **\n",i,op_text[j]);
#endif
#ifdef DEBUG_MUTATE
                  printf(" [%5d %s]",i,op_text[j]);
#endif
                  break;
                }
              else
                t -= v;
            }
          u = erand();
          if (j < 4 && u < h->ins*erate)
            { t = u; 
              goto repeat2;
            }
        } 
      else
        { v *= erate;
          u = (u-v)/(1.-v);
        }

#ifdef DEBUG_MUTATE
      printf("\n");
#endif

      if (Hex[z] != m)
        { m = Hex[z];
          k = (i+1)-2*m;
        }
      x = z;
    }

  seq[len]   = e1;
  seq[len+1] = e2;
  seq[len+2] = e3;

  return (etop);
}


/*******************************************************************************************
 *
 *  Read source genome, return in Genome record with all contigs and scaffolds cat'd together
 *
 ********************************************************************************************/

typedef struct
  { int64   sfnum;   //  # of scaffold
    int64   nbase;   //  # of bases
    int64  *sflen;   //  sflen[i] = length of scaffold i in [0,sfnum)
    uint8 **scafs;   //  scafs[i] = ptr to 2-bit compressed sequence of scaffold i
  } Genome;

static char *fmer[256], _fmer[1280];

static void setup_fmer_table()
{ char *t;
  int   i, l3, l2, l1, l0;

  i = 0;
  t = _fmer;
  for (l3 = 0; l3 < 4; l3++)
   for (l2 = 0; l2 < 4; l2++)
    for (l1 = 0; l1 < 4; l1++)
     for (l0 = 0; l0 < 4; l0++)
       { fmer[i] = t;
         *t++ = dna[l3];
         *t++ = dna[l2];
         *t++ = dna[l1];
         *t++ = dna[l0];
         *t++ = 0;
         i += 1;
       }
}

static void Free_Genome(Genome *gene)
{ free(gene->sflen);
  free(gene->scafs[0]);
  free(gene->scafs);
  free(gene);
}

static int64 Size_Of_Genome(Genome *gene)
{ int64 nbps, size;
  int   i;

  nbps = 0;
  for (i = 0; i < gene->sfnum; i++)
    nbps  += (gene->sflen[i]+3)/4;
  size = sizeof(Genome) + gene->sfnum*(sizeof(int64)+sizeof(uint8 *)) + nbps;
  return ((size-1)/0x100000+1);
}

static void Print_Genome(Genome *gene, FILE *file)
{ int64 i, j, k, b;

  setup_fmer_table();

  // WIDTH must be >= 4

  for (i = 0; i < gene->sfnum; i++)
    { int64  len  = gene->sflen[i];
      uint8 *base = gene->scafs[i];
      int64  lenw = len-WIDTH;

      fprintf(file,"> Scaffold %lld of %lld (len = %lld)\n",i+1,gene->sfnum,len);

      k = 0;
      while (k < len)
        { b = (k>>2);
          if (k < lenw)
            { fprintf(file,"%s",fmer[base[b]] + (k&0x3));
              k += WIDTH;
            }
          else
            { fprintf(file,"%.*s",(int) (len-k),fmer[base[b]] + (k&0x3));
              k = len;
            }
          j = b+1;
          for (b = (k>>2); j < b; j++)
            fprintf(file,"%s",fmer[base[j]]);
          if (j <= b)
            { b = (k&0x3);
              if (b > 0)
                fprintf(file,"%.*s",(int) b,fmer[base[j]]);
            }
          fprintf(file,"\n");
        }
    }
  fflush(file);
}

Genome *Load_Genome(char *name, char **core)
{ Genome    *db;
  Entry      entry;
  FILE      *gene;

  if (VERBOSE)
    fprintf(stderr,"\n  Loading genome %s\n",name);
  fflush(stderr);

  //  Auto-complete extensions on name, set Fetch acording to fasta or fastq

  { char *pwd, *root;
    int   i;
    char *suffix[] = { ".fasta", ".fa", ".fastq", ".fq" };

    pwd = PathTo(name);
    for (i = 0; i < 4; i++)
      { root = Root(name,suffix[i]);
        gene = fopen(Catenate(pwd,"/",root,suffix[i]),"r");
        if (gene != NULL)
          break;
        free(root);
      }
    if (gene == NULL)
      { fprintf(stderr,"%s: Cannot open genome %s with extension .f[ast][aq]\n",
                       Prog_Name,name);
        exit (1);
      }
    *core = root;
    if (i < 2)
      Fetch = Get_Fasta_Entry;
    else
      Fetch = Get_Fastq_Entry;
  }

  //  Scan genome wihtout loading to get # scaffolds & # bps

  entry.seq    = NULL;
  entry.header = NULL;

  { int64   sfnum, *sflen;
    uint8 **scafs, *bases;
    int64   slen, smax, nbps;
 
    smax  = 0;
    nbps  = 0;
    sfnum = 0;
    while (Fetch(gene,&entry) != NULL)
      { slen = entry.slen;
        if (slen > smax)
          smax = slen;
        sfnum += 1;
        nbps  += (slen+3)/4;
      }

    //  Allocate now that you know sizes

    db    = Malloc(sizeof(Genome),"Allocating genome");
    sflen = Malloc(sizeof(int64)*sfnum,"Allocating genome");
    scafs = Malloc(sizeof(uint8 *)*sfnum,"Allocating genome");
    bases = Malloc(sizeof(uint8)*nbps,"Allocating genome");
    if (db == NULL || sflen == NULL || scafs == NULL || bases == NULL)
      exit (1);

    db->sfnum    = sfnum;
    db->sflen    = sflen;
    db->scafs    = scafs;
    db->scafs[0] = bases;

    entry.seq = Malloc(smax+4,"Allocating genome");
    if (entry.seq == NULL)
      exit (1);
  }

  //  In a second scan load and compress the scaffolds

  rewind(gene);

  { int64   sfnum, *sflen, nbase;
    uint8 **scafs, *bases;
    char   *seq;
    int64   slen, i, j;
    int     u;

    seq   = entry.seq;
    bases = db->scafs[0];
    sflen = db->sflen;
    scafs = db->scafs;

    nbase = 0;
    sfnum = 0;
    while (Fetch(gene,&entry) != NULL)
      { scafs[sfnum] = bases;
        slen = entry.slen;

        //  Cull N's and convert bases to numbers

        j = 0;
        for (i = 0; i < slen; i++)
          { u = Number[(int) seq[i]];
            if (u >= 0)
              seq[j++] = u;
           }
        sflen[sfnum++] = slen = j;

        nbase += slen;

        //  Pack scaffold into 2-bits per base

        seq[slen] = seq[slen+1] = seq[slen+2] = 0;
        for (i = 0; i < slen; i += 4)
          *bases++ = (seq[i] << 6) | (seq[i+1] << 4) | (seq[i+2] << 2) | seq[i+3];
      }

    db->nbase = nbase;

    free(entry.seq);

    //  Pack compressed array of bases, if it moves adjust all scaf ptrs

    bases = Realloc(scafs[0],bases-scafs[0],"Refitting packed genome");
    if (bases != scafs[0])
      { int64 move = bases-scafs[0];
        for (i = 0; i < sfnum; i++)
          scafs[i] += move;
      } 
  }

  if (VERBOSE)
    fprintf(stderr,"    %lld scaffolds, %lldMbp (memory = %lldMB)\n",
                   db->sfnum,(db->nbase-1)/1000000+1,Size_Of_Genome(db));

  return (db);
}


/*******************************************************************************************
 *
 *  Routines to create, mutate, and display haplotypes
 *
 ********************************************************************************************/

typedef struct _block
  { int64   beg;    //  start position of block in ancestral genome
    int64   cum;    //  # of bases before start of this block in the haplotype
    uint32 *snp;    //  ptr to 1st SNP in this block
  } Block;

typedef struct
  { Block  **blocks;
    uint32  *snps;
    Genome  *gene;
    int      max_blk;
    double   rate;
  } Haplotype;

static void Free_Haplotype(Haplotype *hap)
{ free(hap->blocks[0]);
  free(hap->blocks);
  free(hap->snps);
  free(hap);
}

static int64 Size_Of_Haplotype(Haplotype *hap)
{ int64 size;

  size = sizeof(Haplotype) + (hap->gene->sfnum+1)*sizeof(Block *)
       + ((hap->blocks[hap->gene->sfnum] - hap->blocks[0])+1) * sizeof(Block)
       + (hap->blocks[hap->gene->sfnum]->snp - hap->blocks[0]->snp) * sizeof(uint32);
  return ((size-1)/0x100000+1);
}

static void Print_Haplotype(Haplotype *hap, FILE *file)
{ int     i;
  int64   off, beg, end;
  uint32 *s;
  Block  *b, *e;

  off = 0;
  fprintf(file,"\nHaplotype %g:\n",hap->rate);
  for (i = 0; i < hap->gene->sfnum; i++)
    { b = hap->blocks[i];
      e = hap->blocks[i+1];
      fprintf(file,"  Scaffold %d, len = %lld, blocks = %ld, snps = %ld\n",
                   i+1,e->cum-b->cum,e-b,e->snp-b->snp);
      for ( ; b < e; b++)
        { beg = b->beg-off;
          end = beg + (b[1].cum - b->cum);
          fprintf(file,"    %10lld - %10lld  [%10lld]\n",beg,end,b[1].cum);
          for (s = b->snp; s < b[1].snp; s++)
            fprintf(file,"      %9u +%1d\n",((*s)>>2)+1,(*s)&0x3);
        }
      off += (4 - (hap->gene->sflen[i] & 0x3)) & 0x3;
    }
  fflush(file);
}

static void Output_Haplotype(Haplotype *hap, FILE *file)
{ int     i;
  int64   off, beg, end;
  uint32 *s;
  Block  *b, *e;

  off = 0;
  fprintf(file,"1 3 hap\n");
  for (i = 0; i < hap->gene->sfnum; i++)
    { b = hap->blocks[i];
      e = hap->blocks[i+1];
      fprintf(file,"c %ld %lld\n", e-b, e->cum - b->cum);
      for ( ; b < e; b++)
        { beg = b->beg-off;
          end = beg + (b[1].cum - b->cum);
          fprintf(file,"B %10lld %10lld\n",beg,end);
          if (b->snp < b[1].snp)
            { fprintf(file,"P %ld",b[1].snp-b->snp);
              for (s = b->snp; s < b[1].snp; s++)
                fprintf(file," %u",((*s)>>2)+1);
              fprintf(file,"\n");
              fprintf(file,"S %ld",b[1].snp-b->snp);
              for (s = b->snp; s < b[1].snp; s++)
                fprintf(file," %1d",(*s)&0x3);
              fprintf(file,"\n");
            }
        }
      off += (4 - (hap->gene->sflen[i] & 0x3)) & 0x3;
    }
  fflush(file);
}

static int get_sequence(uint8 *pack, int64 beg, int64 len, uint8 *seq)
{ int64 p, s;
  int   k;

  p = (beg >> 2);
  k = 6-2*(beg&0x3);
  for (s = 0; s < len; s++)
    { seq[s] = (pack[p]>>k) & 0x3; 
      if (k == 0)
        { p += 1; k = 6; }
      else
        k -= 2;
    }
  return ((int) len);
}

static void mutate_block(uint8 *seq, uint32 *snps, int len)
{ uint32 w, p, o;
  int    i;

  for (i = 0; i < len; i++)
    { w = snps[i];
      p = (w >> 2);
      o = (w & 0x3);
      seq[p] = (seq[p]+o) & 0x3;
    }
}

static void show_block(uint8 *seq, Block *blk)
{ int64   b, e;
  uint32 *c, *d;
  int     s, nu;

  b = blk->beg;
  e = b + (blk[1].cum-blk->cum);
  c = blk->snp;
  d = blk[1].snp;
  if (c < d)
    nu = ((*c++)>>2);
  else
    nu = -1;
  printf("%*s",(int) (b%WIDTH),"");
  for (s = 0; b < e; s++)
    { if (s == nu)
        { printf("%c",dna[seq[s]]-32);
          if (c < d)
            nu = ((*c++)>>2)-1;
          else
            nu = -1;
        }
      else
        printf("%c",dna[seq[s]]);
      if (((++b) % WIDTH) == 0)
        printf("\n");
    }
  if ((b % WIDTH) != 0)
    printf("\n");
}

static Genome *Haplotype_Sequence(Haplotype *hap)
{ int64   sfnum  = hap->gene->sfnum; 
  Block **blocks = hap->blocks;
  uint8  *sbase  = hap->gene->scafs[0];

  Genome *db;
  int64  *sflen;
  uint8 **scafs, *bases, *seq;
  int64   p, nbps;
  Block  *b;
  int     i, r, s, len;

#ifndef DEBUG_HAPLO_MAKE
  (void) show_block;
#endif

  nbps = 0;
  for (i = 0; i < sfnum; i++)
    nbps += ((blocks[i+1]->cum - blocks[i]->cum)+3)/4;

  db    = Malloc(sizeof(Genome),"Allocating haplotype sequence");
  sflen = Malloc(sizeof(int64)*sfnum,"Allocating haplotype sequence");
  scafs = Malloc(sizeof(uint8 *)*sfnum,"Allocating haplotype sequence");
  bases = Malloc(sizeof(uint8)*nbps,"Allocating haplotype sequence");
  seq   = Malloc(hap->max_blk,"Allocating haplotype sequencea");

  if (db == NULL || sflen == NULL || scafs == NULL || bases == NULL || seq == NULL)
      exit (1);

  for (i = 0; i < sfnum; i++)
    sflen[i] = blocks[i+1]->cum - blocks[i]->cum;

  p = 0;
  r = 6;
  for (i = 0; i < sfnum; i++)
    { scafs[i] = bases+p;

      for (b = blocks[i]; b < blocks[i+1]; b++)
        { len = get_sequence(sbase,b->beg,b[1].cum-b->cum,seq);
          mutate_block(seq,b->snp,b[1].snp-b->snp);
#ifdef DEBUG_HAPLO_MAKE
          show_block(seq,b);
#endif

          s = 0;
          while (r != 6 && s < len)
            { bases[p] |= (seq[s++] << r);
              if (r > 0)
                r -= 2;
              else
                { r  = 6;
                  p += 1;
                }
            }
          while (s+4 <= len)
            { bases[p++] = (seq[s] << 6) | (seq[s+1] << 4) | (seq[s+2] << 2) | seq[s+3];
              s += 4;
            }
          bases[p] = 0;
          while (s < len)
            { bases[p] |= (seq[s++] << r);
              if (r > 0)
                r -= 2;
              else
                { r  = 6;
                  p += 1;
                }
            }
        }
      if (r != 6)
        { p += 1; r = 6; }
    }

  free(seq);

  db->sfnum = sfnum;
  db->sflen = sflen;
  db->scafs = scafs;
  db->nbase = blocks[sfnum]->cum - blocks[0]->cum;

  return (db);
}

static Haplotype *init_haplotype(Genome *gene)
{ Haplotype *hap;
  Block     *blocks, **bptrs;
  uint32    *snps;
  int        i, numb, mblk;
  int64      j, acum, len;

  //  Must cut each scaffold into blocks of size no longer than 1Gbp
  //    so that snps can be recorded in a uint32

  numb = 0;
  for (i = 0; i < gene->sfnum; i++)
    numb += (gene->sflen[i]-1)/0x40000000+1;

  blocks = Malloc(sizeof(Block)*(numb+1),"Allocating Haplotype");
  bptrs  = Malloc(sizeof(Block *)*(gene->sfnum+1),"Allocating Haplotype");
  snps   = Malloc(sizeof(uint32),"Allocating Haplotype");
  hap    = Malloc(sizeof(Haplotype),"Allocating Haplotype");
  if (blocks == NULL || bptrs == NULL || snps == NULL || hap == NULL)
    exit (1);

  acum = 0;
  numb = 0;
  mblk = 0;
  for (i = 0; i < gene->sfnum; i++)
    { bptrs[i] = blocks+numb;
      for (j = 0; j < gene->sflen[i]; j += 0x40000000)
        { blocks[numb].beg = 4*(gene->scafs[i]-gene->scafs[0])+j;
          if (j+0x40000000 > gene->sflen[i])
            len = gene->sflen[i] - j;
          else
            len = 0x40000000;
          blocks[numb].cum = acum;
          blocks[numb].snp = snps;
          acum += len;
          numb += 1;
          if (len > mblk)
            mblk = len;
        }
    }
  blocks[numb].cum   = acum;
  blocks[numb].snp   = snps;
  bptrs[gene->sfnum] = blocks+numb;

  hap->gene    = gene;
  hap->blocks  = bptrs;
  hap->snps    = snps;
  hap->max_blk = mblk;
  hap->rate    = 0.;
  
  return (hap);
}

static Haplotype *mutate_haplotype(double rate, Haplotype *hap)
{ Block **blocks = hap->blocks;
  int64   sfnum  = hap->gene->sfnum;

  Block  *C, *D, *I;
  int64   cpt, cbe, ipt;
  int64   rlen, glen;
  int64   dist, len, mlen;
  int     m, isdel, nmut;
  double  x;
  uint32 *csp, *isp, a, rel;

  int64      nint, nlen, isnp, nsnp, nxtcum;
  int        nscaf, mblk;
  Haplotype *nap;
  Block     *nlocks, **nbptrs;
  uint32    *nsnpts;

  uint64 state = GENER;

  glen = hap->blocks[sfnum]->cum;

#ifdef DEBUG_HAPLO
  printf("\nMutating at %g%%\n",100.*rate);
  printf("    %lld scaffolds totaling = %lld bps\n",sfnum,glen);
#endif

  //  In a first pass, generate mutation lengths until have rate*glen total
  //    recording how many mutations are needed.

  rlen = rate *glen;
  mlen = 0;
  nmut = 0;
  while (mlen < rlen)
    { erand();
      mlen += sample_pareto(PARETO_MIN,PARETO_SHAPE);
      if (erand() > SV_RANDOM)
        erand();
      nmut += 1;
    }

#ifdef DEBUG_HAPLO
  printf("    Introducing %d SVs totalling %lld bases\n",nmut,mlen);
#endif

  //  In the next scan, determine number of intervals in result (nint),
  //    total bps in result (ncum), and number of inherited SNPs (isnp);

  GENER = state;

  nlen = 0;
  nint = 0;
  isnp = 0;

  rlen = glen;
  mlen = 0;

  C   = blocks[0];
  D   = C+1;
  cpt = C->cum;
  csp = C->snp;
  cbe = D->cum;
  for (m = nmut; m > 0; m--)

    { dist = sample_exponential((1.*rlen)/m);
      len  = sample_pareto(PARETO_MIN,PARETO_SHAPE);
      if (m == 1)
        len = (rate*glen)-mlen;
      if (dist > rlen)
        dist = rlen;
      x = erand();
      if (x < SV_DELETION)
        { isdel = 1;
          if (dist+len > rlen)
            len = rlen-dist;
          if (len <= 0)
            continue;
#ifdef DEBUG_HAPLO_SCAN
          printf("  Del %lld (%lld)\n",cpt+dist,len); fflush(stdout);
#endif
        }
      else
        { isdel = 0;
          if (x > SV_RANDOM)
            { if (len > glen)
                len = glen;
              ipt = erand()*(glen-len);
#ifdef DEBUG_HAPLO_SCAN
              printf("  Ins %lld %lld (%lld)",cpt+dist,ipt,len); fflush(stdout);
#endif
            }
          else
            { ipt = cpt + dist; 
              if (ipt <= 0)
                continue;
              if (len > ipt)
                len = ipt;
              ipt -= len;
#ifdef DEBUG_HAPLO_SCAN
              printf("  Dup %lld %lld (%lld)",cpt+dist,ipt,len); fflush(stdout);
#endif
            }
          I = blocks[0] + ((blocks[sfnum]-blocks[0])*ipt)/glen;
          while (I->cum <= ipt)
            I += 1;
          while (I->cum > ipt) 
            I -= 1;
          isp = I->snp;
          while (isp < I[1].snp && (*isp)>>2 < ipt-I->cum)
            isp += 1;
#ifdef DEBUG_HAPLO_SCAN
          printf("  [%lld,%lld]\n",I->cum,I[1].cum); fflush(stdout);
#endif
        }
    
      mlen += len;
      rlen -= dist;

      while (1)
        { if (cpt + dist > cbe)
            { nlen += cbe-cpt;
              nint += 1;
              isnp += D->snp-csp;
              csp   = D->snp;

              dist -= cbe-cpt; 
              C   = D++;
              cpt = cbe;
              cbe = D->cum;
            }
          else
            break;
        }

      if (dist > 0)
        { nlen += dist;
          nint += 1;
          cpt  += dist;
          while (csp < C[1].snp && (*csp>>2) < cpt-C->cum)
            { csp += 1;
              isnp += 1;
            }
        }

      if (isdel)
        { rlen -= len;
          while (cpt + len >= cbe && cpt + len < glen)
            { len -= cbe-cpt;
              C = D++;
              cpt = cbe;
              cbe = D->cum;
              csp = C->snp;
            }
          cpt += len;
          while (csp < D->snp && (*csp>>2) < cpt-C->cum)
            csp += 1;
        }
      else
        { while (ipt+len > I[1].cum)
            { int64 n = I[1].cum - ipt;
              nint += 1;
              nlen += n;
              I += 1;
              isnp += I->snp - isp;
              isp   = I->snp;

              ipt += n;
              len -= n;
            }
          nlen += len;
          nint += 1;
          while (isp < I[1].snp && (*isp)>>2 < len)
            { isp  += 1;
              isnp += 1;
            }
          if (cpt >= cbe && cpt < glen)
            { C = D++;
              cpt = cbe;
              cbe = D->cum;
              csp = C->snp;
            }
        }
    }
  while (rlen > 0)
    { nint += 1;
      nlen += cbe-cpt;
      isnp += D->snp-csp;
      csp   = D->snp;

      rlen -= cbe-cpt;

      C = D++;
      cpt = cbe;
      cbe = D->cum;
    }

  //  Add in room for new SNPs and generate data structure for haplotype now
  //    that all sizes are known.

  nsnp = rate*nlen; 

#ifdef DEBUG_HAPLO
  printf("    Resulting in %lld blocks totalling %lld bases with %lld inherited snps\n",nint,nlen,isnp);
#endif

  nlocks = Malloc(sizeof(Block)*(nint+1),"Allocating Haplotype");
  nbptrs = Malloc(sizeof(Block *)*(sfnum+1),"Allocating Haplotype");
  nsnpts = Malloc(sizeof(uint32)*(isnp+nsnp),"Allocating Haplotype");
  nap    = Malloc(sizeof(Haplotype),"Allocating Haplotype");
  if (nlocks == NULL || nbptrs == NULL || nsnpts == NULL || nap == NULL)
    exit (1);

  //  In the next pass, mutate the input haplotype using the same pseudo-random sequence
  //    as in previous passes so everything fits perfectly into the new allocated structure.

  GENER = state;

  nlen = 0;
  nint = 0;
  isnp = nsnp;   //  leave room for the new SNPs added in a final pass over the new haplotype

  rlen = glen;
  mlen = 0;

  nbptrs[0] = nlocks;
  nscaf     = 1;
  nxtcum    = blocks[nscaf]->cum;

  C   = blocks[0];
  D   = C+1;
  cpt = C->cum;
  csp = C->snp;
  cbe = D->cum;
  for (m = nmut; m > 0; m--)

    { dist = sample_exponential((1.*rlen)/m);
      len  = sample_pareto(PARETO_MIN,PARETO_SHAPE);
      if (m == 1)
        len = (rate*glen)-mlen;
      if (dist > rlen)
        dist = rlen;
      x = erand();
      if (x < SV_DELETION)
        { isdel = 1;
          if (dist+len > rlen)
            len = rlen-dist;
          if (len <= 0)
            continue;
        }
      else
        { isdel = 0;
          if (x > SV_RANDOM)
            { if (len > glen)
                len = glen;
              ipt = erand()*(glen-len);
            }
          else
            { ipt = cpt + dist; 
              if (ipt <= 0)
                continue;
              if (len > ipt)
                len = ipt;
              ipt -= len;
            }
          I = blocks[0] + ((blocks[sfnum]-blocks[0])*ipt)/glen;
          while (I->cum <= ipt)
            I += 1;
          while (I->cum > ipt) 
            I -= 1;
          isp = I->snp;
          while (isp < I[1].snp && (*isp)>>2 < ipt-I->cum)
            isp += 1;
        }
    
      mlen += len;
      rlen -= dist;

      while (1)
        { if (cpt + dist > cbe)
            { if (cpt >= nxtcum)
                { nbptrs[nscaf++] = nlocks+nint;
                  nxtcum = blocks[nscaf]->cum;
                }
              rel = cpt-C->cum;
              nlocks[nint].beg = C->beg + rel;
              nlocks[nint].cum = nlen;
              nlocks[nint].snp = nsnpts+isnp;
              while (csp < D->snp)
                { uint32 w = *csp++;
                  nsnpts[isnp++] = (((w >> 2) - rel) << 2) | (w & 0x3);
                }
              nlen += cbe-cpt;
              nint += 1;
#ifdef DEBUG_HAPLO_GEN
              printf("  %5lld P: %lld (%lld) %lld {%lld}\n",
                     nint,cpt,cbe-cpt,nlocks[nint-1].beg,isnp);
#endif

              dist -= cbe-cpt; 
              C = D++;
              cpt = cbe;
              cbe = D->cum;
            }
          else
            break;
        }

      if (dist > 0)
        { if (cpt >= nxtcum)
            { nbptrs[nscaf++] = nlocks+nint;
              nxtcum = blocks[nscaf]->cum;
            }
          rel = cpt-C->cum;
          nlocks[nint].beg = C->beg + rel;
          nlocks[nint].cum = nlen;
          nlocks[nint].snp = nsnpts+isnp;
          cpt += dist;
          while (csp < C[1].snp && (*csp>>2) < cpt-C->cum)
            { uint32 w = *csp++;
              nsnpts[isnp++] = (((w >> 2) - rel) << 2) | (w & 0x3);
            }
          nlen += dist;
          nint += 1;
#ifdef DEBUG_HAPLO_GEN
          printf("  %5lld M: %lld (%lld) %lld {%lld}\n",nint,cpt-dist,dist,nlocks[nint-1].beg,isnp);
#endif
        }

      if (isdel)
        { rlen -= len;
          while (cpt + len >= cbe && cpt + len < glen)
            { len -= cbe-cpt;
              C = D++;
              cpt = cbe;
              cbe = D->cum;
              csp = C->snp;
            }
          cpt += len;
          while (csp < D->snp && (*csp>>2) < cpt-C->cum)
            csp += 1;
        }
      else
        { while (ipt+len > I[1].cum)
            { int64 n = I[1].cum - ipt;
              rel = ipt-I->cum;
              nlocks[nint].beg = I->beg + rel;
              nlocks[nint].cum = nlen;
              nlocks[nint].snp = nsnpts + isnp;
              nint += 1;
              nlen += n;
              I += 1;
              while (isp < I->snp)
                { uint32 w = *isp++;
                  nsnpts[isnp++] = (((w >> 2) - rel) << 2) | (w & 0x3);
                }
#ifdef DEBUG_HAPLO_GEN
              printf("  %5lld L: %lld (%lld) %lld {%lld}\n",nint,ipt,n,nlocks[nint-1].beg,isnp);
#endif

	      ipt += n;
              len -= n;
            }
          rel = ipt-I->cum;
          nlocks[nint].beg = I->beg + rel;
          nlocks[nint].cum = nlen;
          nlocks[nint].snp = nsnpts + isnp;
          nlen += len;
          nint += 1;
          while (isp < I[1].snp && (*isp)>>2 < len)
            { uint32 w = *isp++;
              nsnpts[isnp++] = (((w >> 2) - rel) << 2) | (w & 0x3);
            }
#ifdef DEBUG_HAPLO_GEN
          printf("  %5lld I: %lld (%lld) %lld {%lld}\n",nint,ipt,len,nlocks[nint-1].beg,isnp);
#endif

          if (cpt >= cbe && cpt < glen)
            { C = D++;
              cpt = cbe;
              cbe = D->cum;
              csp = C->snp;
            }
        }
    }
  while (rlen > 0)
    { if (cpt >= nxtcum)
        { nbptrs[nscaf++] = nlocks+nint;
          nxtcum = blocks[nscaf]->cum;
        }
      rel = cpt-C->cum;
      nlocks[nint].beg = C->beg + rel;
      nlocks[nint].cum = nlen;
      nlocks[nint].snp = nsnpts + isnp;
      nint += 1;
      nlen += cbe-cpt;
      while (csp < D->snp)
        { uint32 w = *csp++;
          nsnpts[isnp++] = (((w >> 2) - rel) << 2) | (w & 0x3);
        }
#ifdef DEBUG_HAPLO_GEN
      printf("  %5lld E: %lld (%lld) %lld {%lld}\n",nint,cpt,cbe-cpt,nlocks[nint-1].beg,isnp);
#endif

      rlen -= cbe-cpt;
      C = D++;
      cpt = cbe;
      cbe = D->cum;
      csp = C->snp;
    }
  nbptrs[sfnum] = nlocks+nint;
  nlocks[nint].cum = nlen;
  nlocks[nint].snp = nsnpts+isnp;

  //  In a pass over the new haplotype, add nsnps new SNPs to its blocks

#ifdef DEBUG_HAPLO
  printf("    Adding %lld new nsnps for a total of %lld snps\n",nsnp,isnp);
#endif

  C = nbptrs[0];
  D = C+1;
  csp  = C->snp;
  isnp = 0;
  C->snp = nsnpts;
  cpt  = 0;
  mblk = 0;
  rlen = nlen-nsnp;
  for (m = nsnp; m > 0; m--)
    { dist = sample_exponential((1.*rlen)/m);
      if (dist > rlen)
        dist = rlen;
      a = ((uint32) (erand()*3.)) + 1;
      if (a >= 4)
        a = 3;
      cpt += dist;
#ifdef DEBUG_HAPLO_SNPS
      printf("Next location (%lld) +%lld -> %lld an %d",rlen/m,dist,cpt,a);
#endif
      while (cpt >= D->cum)
        { while (csp < D->snp)
            nsnpts[isnp++] = *csp++;  
          if (D->cum - C->cum > mblk)
            mblk = D->cum - C->cum;
          C = D++;
          C->snp = nsnpts + isnp;
        }
      while (csp < D->snp && (*csp>>2) < (cpt-C->cum))
        nsnpts[isnp++] = *csp++;  
      if (csp < D->snp && (*csp>>2) == (cpt-C->cum))
        csp += 1;
#ifdef DEBUG_HAPLO_SNPS
      printf(" in [%lld,%lld] %lld\n",C->cum,D->cum,cpt-D->cum);
#endif
      nsnpts[isnp++] = (((uint32) (cpt-C->cum)) << 2) | a;
      cpt  += 1;
      rlen -= dist; 
#ifdef DEBUG_HAPLO_SNPS
      printf("  rlen = %lld, isnp = %lld, psnp = %ld\n",rlen,isnp,csp-nsnpts);
#endif
    }
  while (C < nbptrs[sfnum])
    { while (csp < D->snp)
        nsnpts[isnp++] = *csp++;  
      if (D->cum - C->cum > mblk)
        mblk = D->cum - C->cum;
      C = D++;
      C->snp = nsnpts + isnp;
    }
#ifdef DEBUG_HAPLO_SNPS
  printf("  rlen = %lld, isnp = %lld, psnp = %ld\n",rlen,isnp,csp-nsnpts);
#endif

  nap->gene    = hap->gene;
  nap->blocks  = nbptrs;
  nap->snps    = nsnpts;
  nap->max_blk = mblk;
  nap->rate    = hap->rate + rate;

#ifdef SHOW_HAPLO
  Print_Haplotype(nap);
#else
  (void) Print_Haplotype;
#endif

  return (nap);
}


/*******************************************************************************************
 *
 *  Parse the ploidy tree
 *
 ********************************************************************************************/

typedef struct _node
  { double  rate;
    struct _node *sib;
    struct _node *sub;
  } Node;

static int         Nhaps;
static Haplotype **Haps;

static char *Scan;
static int   Error;

static char *Error_Messages[] =
  { "Out of memory",            // 0
    "Expecting real number",    // 1
    "Expecting , or )"          // 2
  };

static void free_tree(Node *v)
{ if (v->sib != NULL)
    free_tree(v->sib);
  if (v->sub != NULL)
    free_tree(v->sub);
  free(v);
}

static void skip_white()
{ while (isspace(*Scan))
    Scan += 1;
}

static Node *scan_rate()
{ Node  *n;
  double v;
  char  *eptr;

  v = strtod(Scan,&eptr);
  if (eptr <= Scan)
    { Error = 1;
      return (NULL);
    }
  Scan = eptr;
  n = Malloc(sizeof(Node),"Allocating ploidy tree");
  if (n == NULL)
    { Error = 0;
      return (NULL);
    }
  n->sub = n->sib = NULL;
  n->rate = v/100.;
  return (n);
}

static Node *scan_tree(int close)
{ Node *val, *nxt;

  val = scan_rate();
  if (val == NULL)
    return (NULL);
  skip_white();
  while (*Scan != close)
    { if (*Scan == '(')
        { Scan += 1;
          val->sub = scan_tree(')');
          if (val->sub == NULL)
            { free_tree(val);
              return (NULL);
            }
          Scan += 1;
          skip_white();
        }
      else if (*Scan == ',')
        { Scan += 1;
          Nhaps += 1;
          nxt = scan_rate();
          if (nxt == NULL)
            { free_tree(val);
              return (NULL);
            }
          nxt->sib = val;
          val = nxt;
          skip_white();
        }
      else
        { free_tree(val);
          Error = 2;
          return (NULL);
        }
    }
  nxt = NULL;          //  Reverse sibling list
  while (val != NULL)
    { Node *prv;
      prv = val->sib;
      val->sib = nxt;
      nxt = val;
      val = prv;
    } 
  return (nxt);
}

static void show_tree(Node *tree, int level)
{ Node *v;
 
  for (v = tree; v != NULL; v = v->sib)
    { printf("%*s%g\n",level,"",v->rate);
      if (v->sub != NULL)
        show_tree(v->sub,level+2);
    }
}


/*******************************************************************************************
 *
 *  Interpret the ploidy tree and build the given number of haplotypes, returning
 *    an array of pointers to them along with their number.
 *
 ********************************************************************************************/

static void gen_haps(Node *v, Haplotype *hap)
{ Haplotype *rez;

  for ( ; v != NULL; v = v->sib)
    { rez = mutate_haplotype(v->rate,hap);
      if (v->sub != NULL)
        gen_haps(v->sub,rez);
      else
        { if (VERBOSE)
            fprintf(stderr,"*");
          fflush(stderr);
          Haps[Nhaps++] = rez;
        }
    }
  Free_Haplotype(hap);
}

static Haplotype **Gen_Haplotypes(Genome *root, char *ploidy, int *nhaps)
{ Node *tree;
  Haplotype *hap;

  Nhaps = 1;
  Scan = ploidy;
  skip_white();
  if (*Scan == 0)
    tree = NULL;
  else
    { tree = scan_tree(0);
      if (tree == NULL)
        { if (Error == 0)
            fprintf(stderr,"%s: Out of memory parsing expression\n",Prog_Name);
          else
            { fprintf(stderr,"%s: Expression syntax error (%d):\n\n",Prog_Name,Error);
              fprintf(stderr,"    %s\n",ploidy);
              fprintf(stderr,"%*s^ %s\n",(int) ((Scan-ploidy)+4),"",Error_Messages[Error]);
            }
          exit (1);
        }
#ifdef DEBUG_HAPLO
      printf("\nPloidy Tree:\n");
      show_tree(tree,4);
#else
      (void) show_tree;
#endif
    }

  if (VERBOSE)
    fprintf(stderr,"\n  Generating %d haplotypes: ",Nhaps);
  fflush(stderr);

  Haps = (Haplotype **) Malloc(sizeof(Haplotype *)*Nhaps,"Allocating haplotypes");
  if (Haps == NULL)
    exit (1);

  Nhaps = 0;
  hap = init_haplotype(root);
  if (tree != NULL)
    gen_haps(tree,hap);
  else
    { if (VERBOSE)
        fprintf(stderr,"*");
      fflush(stderr);
      Haps[Nhaps++] = hap;
    }
  if (VERBOSE)
    fprintf(stderr,"\n");
  fflush(stderr);

  if (VERBOSE)
    { int p;
      int bmax, smax, lmax, mmax;

      bmax = smax = lmax = mmax = 0;
      for (p = 0; p < Nhaps; p++)
        { Block *b, *e;
         int     x;

          b = Haps[p]->blocks[0];
          e = Haps[p]->blocks[root->sfnum];
          x = Number_Digits(e-b);
          if (x > bmax)
            bmax = x;
          x = Number_Digits(e->snp-b->snp);
          if (x > smax)
            smax = x;
          x = Number_Digits(((e->cum-b->cum)-1)/1000000+1);
          if (x > lmax)
            lmax = x;
          x = Number_Digits(Size_Of_Haplotype(Haps[p]));
          if (x > mmax)
            mmax = x;
        }

      for (p = 0; p < Nhaps; p++)
        { Block *b, *e;

          b = Haps[p]->blocks[0];
          e = Haps[p]->blocks[root->sfnum];
          fprintf(stderr,"    %2d: %*ld blocks, %*ld snps, %*lldMbp (memory = %*lldMB)\n",
                         p+1,bmax,e-b,smax,e->snp-b->snp,lmax,((e->cum-b->cum)-1)/1000000+1,
                         mmax,Size_Of_Haplotype(Haps[p]));
        }
    }

  *nhaps = Nhaps;
  return (Haps);
}


/*******************************************************************************************
 *
 *  Interpret the ploidy tree and build the given number of haplotypes, returning
 *  Generate reads (a) whose lengths are log-normal distributed with mean *RMEAN* and
 *    standard deviation *RSDEV*, and (b) that are never shorter than *RSHORT*.  Each
 *    read is a randomly sampled interval of the haplotype *gene* that has errors
 *    introduced into it according to the user-supplied error model and in either strand
 *    direction with equal probability.  If the -C option is set then each scaffold is
 *    assumed to be circular and reads can be sampled that span the origin.   Reads are
 *    generated until the sum of the lengths of the reads is greater thant coverage times
 *    the sum of the lengths of the scaffolds in the scaffold.  The reads are output as
 *    fasta entries with the trace information in the header as per the documentation.
 *
 ********************************************************************************************/

//  Complement (in the DNA sense) string *s*.

static void complement(int64 elen, uint8 *s)
{ uint8 *t;
  int    c;

  t = s + (elen-1);
  while (s <= t)
    { c = *s;
      *s = 3-*t;
      *t = 3-c;
      s += 1;
      t -= 1;
    }
}

static int64 Shotgun(Genome *gene, int ploid, double prate)
{ static char normal_op[] = { 'D', 'S', 'I' };
  static char micro_del[] = { 0, 'h', 'z', 't' }; 
  static char micro_ins[] = { 0, 'H', 'Z', 'T' }; 

  int64      glen;
  uint8     *scaf;
  int64      nbeg, rtag;
  int64      nreads;
  int64      totbp, genbp;
  int64      tooshort;
  int        omax, smax, emax, elen;
  uint8     *oseq, *sseq;
  Edit      *ops;
  double     erate;
  int64      emark;
  int        i;

  omax = RMEAN + 5*RSDEV;
  oseq = Malloc(omax+3,"Allocating read buffer");

  emax = (RMEAN + 5*RSDEV)/sizeof(Edit);
  ops  = Malloc(emax*sizeof(Edit),"Allocating edit script");

  smax = RMEAN + 5*RSDEV;
  sseq = Malloc(smax+3,"Allocating mutated read buffer");

  tooshort = RMEAN / COVERAGE;
  if (tooshort < RSHORT)
    tooshort = RSHORT; 

  if (ERRINFO)
    { emark = ftello(ERR_OUT);
      fprintf(ERR_OUT,"h 1000000000 %g\n",prate);
    }

  genbp  = 0;
  nreads = 0;
  for (i = 0; i < gene->sfnum; i++)
    { scaf = gene->scafs[i];
      glen = gene->sflen[i];

      if (glen <= RMEAN+RSDEV)
        { fprintf(stderr,"%s: Scaffold length is less than mean read length + 1 std. deviation!\n",
                         Prog_Name);
          exit (1);
        }

      if (CIRCULAR)
        nbeg  = 0;
      else
        nbeg  = -RMEAN;
      totbp = COVERAGE*(glen-nbeg);
      rtag = 0;
      while (totbp > 0)
        { int64 len, rbeg, rend, del;

          nbeg += sample_exponential((RMEAN*(glen-nbeg))/totbp);

          do
            len = sample_read_length(&erate);
          while (len < RSHORT);

          rbeg = nbeg;
          rend = nbeg + len;
          if (CIRCULAR)
            { if (rend > glen)
                { rend = rend % glen;
                  rtag += 1;
                  if (rtag >= COVERAGE)
                    break;
                }
              if (glen-rbeg < tooshort)
                break;
            }
          else
            { if (rend > glen)
                { rend = glen;
                  rtag += 1;
                  if (rtag >= COVERAGE)
                    break;
                }
              if (rbeg < 0)
                rbeg = 0;
              len = rend-rbeg;
              if (len < tooshort)
                { if (nbeg < 0)
                    continue;
                  else
                    break;
                }
            }

          totbp -= len;

          if (len > omax)
            { omax = 1.2*len + 1000;
              oseq = Realloc(oseq,omax+3,"Allocating read buffer");
              if (oseq == NULL)
                exit (1);
            }

          get_sequence(scaf,rbeg,len,oseq);

          if (erand() >= FLIP_RATE)    //  Complement the string with probability FLIP_RATE.
            { complement(len,oseq);
              fprintf(READ_OUT,">Sim %d %d - %lld %lld\n",ploid,i+1,rbeg,rend);
              if (ERRINFO)
                fprintf(ERR_OUT,"S %d %d 1 %lld %lld\n",ploid,i+1,rbeg,rend);
#ifdef DEBUG_OPS
              if (READ_OUT != stdout)
                printf(">Sim %d %d - %lld %lld\n",ploid,i+1,rbeg,rend);
#endif
            }
          else
            { fprintf(READ_OUT,">Sim %d %d + %lld %lld\n",ploid,i+1,rbeg,rend);
              if (ERRINFO)
                fprintf(ERR_OUT,"S %d %d 0 %lld %lld\n",ploid,i+1,rbeg,rend);
#ifdef DEBUG_OPS
              if (READ_OUT != stdout)
                printf(">Sim %d %d - %lld %lld\n",ploid,i+1,rbeg,rend);
#endif
            }

          //  Generate errors and output CIGAR string

          while ((elen = mutate_read(oseq,len,ops,emax,erate/AvErr)) < 0)
            { emax = 1.2*emax+100;
              ops  = Realloc(ops,sizeof(Edit)*emax,"Allocating edit script");
              if (ops == NULL)
                exit (1);
            }

          { int i;

            del = 0;
            for (i = 0; i < elen; i++)
              del += ops[i].del;
            if (len+del > smax)
              { smax = 1.2*(len+del) + 1000;
                sseq = Realloc(sseq,omax+3,"Allocating read buffer");
              }
          }

          //  Create the erroneous read according to edits in ops

          { int   i, j, l,  m, n, u;
            Edit *co;

            l = 0;
            j = 0;
            for (i = 0; i < elen; i++)
              { co = ops+i;
                while (l < co->pos)
                  sseq[j++] = oseq[l++];
                if (co->kind == NORMAL)
                  { if (co->del > 0)
                      { if (l <= co->pos)
                          sseq[j++] = oseq[l++];
                        sseq[j++] = co->data;
                      }
                    else if (co->del == 0)
                      { sseq[j++] = co->data;
                        l += 1;
                      }
                  }
                else
                  { if (co->del < 0)
                      j += co->del+1;
                    else
                      { n = co->del;
                        m = co->data;
                        sseq[j++] = oseq[l++];
                        while (n >= m)
                          { for (u = l-m; u < l; u++)
                              sseq[j++] = oseq[u]; 
                            n -= m;
                          }
                        for (u = l-m; n > 0; u++)
                          { sseq[j++] = oseq[u]; 
                            n -= 1;
                          }
                      }
                  }
              }
            while (l < len)
              sseq[j++] = oseq[l++];
   
            sseq[j] = 4;
          }

          //  If -e output edit script

          if (ERR_OUT)
            { int   i, j;
              Edit *co;

              fprintf(ERR_OUT,"O %d ",elen);
              for (i = 0; i < elen; i++)
                { co = ops+i;
                  if (co->kind == NORMAL)
                    fprintf(ERR_OUT,"%c",normal_op[co->del+1]);
                  else if (co->del < 0)
                    fprintf(ERR_OUT,"%c",micro_del[co->data]);
                  else
                    fprintf(ERR_OUT,"%c",micro_ins[co->data]);
                }
              fprintf(ERR_OUT,"\n");

              fprintf(ERR_OUT,"L %d",2*elen+1);
              j = -1;
              for (i = 0; i < elen; i++)
                { co = ops+i;
                  if (co->del > 0)
                    fprintf(ERR_OUT," %d",co->pos-j);
                  else if (co->del < 0)
                    fprintf(ERR_OUT," %d",(co->pos-j)+co->del);
                  else
                    fprintf(ERR_OUT," %d",(co->pos-j)-1);
                  if (co->kind == NORMAL)
                    fprintf(ERR_OUT," %d",co->data);
                  else
                    fprintf(ERR_OUT," %d",abs(co->del));
                  j = co->pos;
                }
              fprintf(ERR_OUT," %lld\n",len-j);
            }

          //  Output the simulated read

          if (UPPER)
            Upper_Read((char *) sseq);
          else
            Lower_Read((char *) sseq);

          { int j;

            len += del;
            for (j = 0; j+WIDTH < len; j += WIDTH)
              fprintf(READ_OUT,"%.*s\n",WIDTH,sseq+j);
            if (j < len)
              fprintf(READ_OUT,"%s\n",sseq+j);
          }

          nreads += 1;
        }

      if (CIRCULAR)
        genbp += COVERAGE*glen - totbp;
      else
        genbp += COVERAGE*(glen+RMEAN) - totbp;
    }

  fseek(ERR_OUT,emark,SEEK_SET);
  fprintf(ERR_OUT,"0%10lld",nreads);
  fseek(ERR_OUT,0,SEEK_END);

  free(oseq);

  return (genbp);
}


/*******************************************************************************************
 *
 *  Main
 *
 ********************************************************************************************/

int main(int argc, char *argv[])
{ Genome     *source;
  char       *PLOIDY;
  int         nhaps;
  Haplotype **haps;

  //  Process command line

  { int    i, j, k;
    int    flags[128];
    char  *eptr;

    ARG_INIT("jass");

    COVERAGE  = 50.;
    RMEAN     = -1;
    RSDEV     = -1;
    RSHORT    = 0;
    SEED      = getpid();
    READ_OUT  = stdout;
    WIDTH     = 100;
    PLOIDY    = "";
    OUT       = NULL;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("ehfvCU");
            break;
          case 'c':
            ARG_REAL(COVERAGE)
            if (COVERAGE < 0.)
              { fprintf(stderr,"%s: Coverage must be non-negative (%g)\n",Prog_Name,COVERAGE);
                exit (1);
              }
            break;
          case 'm':
            ARG_POSITIVE(RMEAN,"Mean read length")
            break;
          case 'o':
            OUT = argv[i]+2;
            break;
          case 'p':
            PLOIDY = argv[i]+2;
            break;
          case 'r':
            ARG_NON_NEGATIVE(SEED,"Random number seed")
            break;
          case 's':
            ARG_NON_NEGATIVE(RSDEV,"Read length standard deviation")
            break;
          case 'w':
            ARG_NON_NEGATIVE(WIDTH,"Line width")
            break;
          case 'x':
            ARG_NON_NEGATIVE(RSHORT,"Read length minimum")
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    CIRCULAR = flags['C'];
    UPPER    = flags['U'];
    VERBOSE  = flags['v'];

    ERRINFO  = flags['e'];
    HAPINFO  = flags['h'];
    FASTOUT  = 1 - flags['f'];

    if (argc != 3)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage[0]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[1]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[2]);
        fprintf(stderr,"\n");
        fprintf(stderr,"      -p: haplotype evolutionary tree (in percent)\n");
        fprintf(stderr,"      -o: root for output files (default genome file root)\n");
        fprintf(stderr,"      -h: generate haplotype ground truth in <o>.hap#\n");
        fprintf(stderr,"      -e: generate read ground truth in <o>.error\n");
        fprintf(stderr,"      -f: output reads to <o>.fasta (default stdout)\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -c: coverage of genome\n");
        fprintf(stderr,"      -m: average read length (log normal distribution).\n");
        fprintf(stderr,"      -s: standard deviation of read lengths (log normal)\n");
        fprintf(stderr,"      -x: ignore reads below this length\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -v: verbose output\n");
        fprintf(stderr,"      -C: assume genome is circular (default is linear)\n");
        fprintf(stderr,"      -U: Use upper case for DNA (default is lower case).\n");
        fprintf(stderr,"      -r: Random number generator seed (default is process id).\n");
        fprintf(stderr,"      -w: Print -w bp per line (default is 80).\n");
        exit (1);
      }
  }

  (void) Print_Genome;

  //  Seed the random # generator, load the error model and genome, and
  //    then generate haplotypes

  eseed(SEED);

  Load_Error_Model(argv[2]);

  if (RMEAN < 0)
    RMEAN = RMean;
  if (RSDEV < 0)
    RSDEV = RsDev;
  Rscale  = RSDEV / RsDev;
  Roffset = RMEAN - Rscale * RMean;

  { char *root;

    source = Load_Genome(argv[1],&root);
    if (OUT == NULL)
      OUT = root;
    else
      free(root);
    if (FASTOUT)
      READ_OUT = stdout;
    else
      READ_OUT = fopen(Catenate(OUT,".fasta","",""),"w");
    if (ERRINFO)
      ERR_OUT = fopen(Catenate(OUT,".err","",""),"w");
  }

  haps = Gen_Haplotypes(source,PLOIDY,&nhaps);

  //  If the -t option set then output haplotype trace to appropriately named files

  if (HAPINFO)
    { int   p;
      FILE *f;

      for (p = 0; p < nhaps; p++)
        { f = fopen(Catenate(OUT,Numbered_Suffix(".hap",p+1,""),"",""),"w");
          if (f == NULL)
            { fprintf(stderr,"%s: Cannot create %s.%d.hap for writing\n",Prog_Name,OUT,p+1);
              exit (1);
            }
          Output_Haplotype(haps[p],f);
          fclose(f);
        }
    }

  //  Generate each haplotype sequence in turn and shotgun reaads from it with the supplied
  //    error model to the requested coverage

  { int     p;
    Genome *ahap;
    int64   totbp;

    if (VERBOSE)
      { fprintf(stderr,"\n  Begin Shotgun at coverage %.1fX\n",COVERAGE);
        if (CIRCULAR)
          fprintf(stderr,"      Assume molecules are circular\n");
        else
          fprintf(stderr,"      Assume molecules are linear\n");
        fprintf(stderr,"      Model-based read lengths with mean %d, sigma %d, and min cutoff %d\n",
                       RMEAN,RSDEV,RSHORT);
      }

    COVERAGE /= nhaps;
    for (p = 0; p < nhaps; p++)
      { if (VERBOSE)
          fprintf(stderr,"    Gen'ing sequence of haplotype %d\n",p+1);
        ahap = Haplotype_Sequence(haps[p]);
        if (VERBOSE)
          { fprintf(stderr,"    Sampling hap %d to depth %.1fX\n",p+1,COVERAGE);
            fflush(stderr);
          }
        totbp = Shotgun(ahap,p+1,haps[p]->rate);
        if (VERBOSE)
          { fprintf(stderr,"      Actually sampled %.1fX\n",(1.*totbp)/ahap->nbase);
            fflush(stderr);
          }
        Free_Genome(ahap);
      }

    if (VERBOSE)
      fprintf(stderr,"  End Shotgun\n\n");
  }

  if (ERRINFO)
    fclose(ERR_OUT);
  if (!FASTOUT)
    fclose(READ_OUT);

  exit (0);
}
