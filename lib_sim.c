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
#define  DEBUG_OPS

#include "gene_core.h"
#include "lib_sim.h"

#define WIDTH 100


/*******************************************************************************************
 *
 *  Read in error model produced by FKmodel
 *
 ********************************************************************************************/

typedef struct
  { float  all;
    float  op[9];
  } E_Rates;

typedef struct
  { float all;
    float op[6];
  } M_Rates;

typedef struct
  { float  all;
    float  mic[4];
  } Seven;

typedef struct
  { int    kmax[4];
    Seven  heptab[0x4000];
    float *mic3tab[0x40];
    float *mic2tab[0x10];
    float *mic1tab[0x04];
  } _Error_Model;

static char  Bell_Adapter[100] = "ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT";
static int   Bell_Len;
static char *Bell_End;

static int Hex[0x1000];
static int HxB[0x1000];

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

Error_Model *Load_Error_Model(char *name)
{ _Error_Model *epro;
  FILE            *efile;
  int              krange;

  { char *pwd, *root;

    pwd   = PathTo(name);
    root  = Root(name,".emod");
    efile = fopen(Catenate(pwd,"/",root,".emod"),"r");
    if (efile == NULL)
      { fprintf(stderr,"%s: Cannot open FK error model %s with extension .emod\n",
                       Prog_Name,name);
        return (NULL);
      }
    free(root);
    free(pwd);
  }

  { M_Rates  mrate;
    E_Rates  erate;
    float   *mic1, *mic2, *mic3;
    int      i, j, kmer;

    if (fread(&kmer,sizeof(int),1,efile) != 1)
      goto error_exit;

    krange = kmer/2 - 6;

    epro = (_Error_Model *) Malloc(sizeof(_Error_Model),"Error Profiler");
    mic1 = (float *) Malloc(sizeof(float)*4*krange,"Homo-Table");
    mic2 = (float *) Malloc(sizeof(float)*16*krange,"Di-Table");
    mic3 = (float *) Malloc(sizeof(float)*64*krange,"Tri-Table");

    for (i = 0; i < 0x4000; i++)
      { if (fread(&erate,sizeof(E_Rates),1,efile) != 1)
          goto error_exit;
        epro->heptab[i].all = erate.all;
        for (j = 0; j < 4; j++)
          epro->heptab[i].mic[j] = erate.all - (erate.op[8] + erate.op[j]);
      }

    mic1 -= 2;
    mic2 -= 4;
    mic3 -= 6;

    epro->kmax[1] = 2+krange;
    epro->kmax[2] = 4+krange;
    epro->kmax[3] = 6+krange;

    for (i = 0; i < 0x04; i++)
      epro->mic1tab[i] = mic1 + krange*i;

    for (i = 0; i < 0x10; i++)
      epro->mic2tab[i] = mic2 + krange*i;

    for (i = 0; i < 0x40; i++)
      epro->mic3tab[i] = mic3 + krange*i;

    for (i = 0; i < 4; i++)
      for (j = 2; j < epro->kmax[1]; j++)
        { if (fread(&mrate,sizeof(M_Rates),1,efile) != 1)
            goto error_exit;
          epro->mic1tab[i][j] = mrate.op[1]/(1.-mrate.op[1])
                              + mrate.op[0]/(1.-mrate.op[0]);
        }

    for (i = 0; i < 16; i++)
      for (j = 4; j < epro->kmax[2]; j++)
        { if (fread(&mrate,sizeof(M_Rates),1,efile) != 1)
            goto error_exit;
          epro->mic2tab[i][j] = (mrate.op[2] + mrate.op[3])/(1.-mrate.op[2])
                              + (mrate.op[0] + mrate.op[1])/(1.-mrate.op[1]);
        }

    for (i = 0; i < 64; i++)
      for (j = 6; j < epro->kmax[3]; j++)
        { if (fread(&mrate,sizeof(M_Rates),1,efile) != 1)
            goto error_exit;
          epro->mic3tab[i][j] = (mrate.op[3] + mrate.op[4] + mrate.op[5])/(1.-mrate.op[5])
                              + (mrate.op[0] + mrate.op[1] + mrate.op[2])/(1.-mrate.op[2]);
        }

    fclose(efile);
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

    printf("\n");
    for (i = 0; i < 0x4000; i++)
      { for (j = 12; j >= 0; j -= 2)
          printf("%c",dna[(i>>j)&0x3]);
        printf(" %5.3f :: ",100.*epro->heptab[i].all);
        for (j = 0; j < 4; j++)
          if (epro->heptab[i].all == epro->heptab[i].mic[j])
            printf("  --- ");
          else
            printf(" %5.3f",100.*epro->heptab[i].mic[j]);
        printf("\n");
      }

    for (n = 1; n <= 3; n++)
      { printf("\n");
        for (i = 0; i < mic_size[n]; i++)
          { if (n == 1)
              mtab = Mic1Tab[i];
            else if (n == 2)
              mtab = Mic2Tab[i];
            else
              mtab = epro->mic3tab[i];
            for (k = 2*(n-1); k >= 0; k -= 2)
              printf("%c",dna[(i>>k)&0x3]);
            printf(" %5.3f\n",mtab[k]);
          }
      }
  }
#endif

  return ((Error_Model *) epro);

error_exit:
  fprintf(stderr,"%s: Read of .emod file %s failed\n",Prog_Name,name);
  exit (1);
}

float *Read_Error_Profile(Error_Model *_epro, uint8 *seq, int len, float *profile)
{ _Error_Model *epro = (_Error_Model *) _epro;

  int      e1, e2, e3;
  Seven   *h;
  float   *g;
  uint32   x, y, z;
  int      i, m, k;
  int      ma;

  if (profile == NULL)
    { profile = (float *) Malloc(sizeof(float)*len,"Read Profile");
      if (profile == NULL)
        return (NULL);
    }

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

  m   = Hex[x];
  k   = -2*m;
  for (i = 0; i < len; i++)
    { z = ((x << 2) | seq[i+1]) & 0xfff;
      y = ((y << 2) | seq[i+3]) & 0x3fff;

      h = epro->heptab+y;
      if (m > 0 && Hex[z] != m)
        { k = i-k;
          if (k > epro->kmax[m])
            k = epro->kmax[m];
          if (m == 1)
            g = epro->mic1tab[x&0x3];
          else if (m == 2)
            g = epro->mic2tab[x&0xf];
          else
            g = epro->mic3tab[x&0x3f];
          ma = HxB[x];
          profile[i] = g[k] + h->mic[ma];
        }
      else
        profile[i] = h->all;

      if (Hex[z] != m)
        { m = Hex[z];
          k = (i+1)-2*m;
        }
      x = z;
    }

  seq[len]   = e1;
  seq[len+1] = e2;
  seq[len+2] = e3;

  return (profile);
}


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
 *  Read source genome, return in Genome record with all contigs and scaffolds cat'd together
 *
 ********************************************************************************************/

typedef struct
  { int64   sfnum;   //  # of scaffold
    int64   nbase;   //  # of bases
    int64  *sflen;   //  sflen[i] = length of scaffold i in [0,sfnum)
    uint8 **scafs;   //  scafs[i] = ptr to 2-bit compressed sequence of scaffold i
  } _Genome;

static char *fmer[256], _fmer[1280];

static char dna[4] = { 'a', 'c', 'g', 't' };

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

void Free_Genome(Genome *_gene)
{ _Genome *gene = (_Genome *) _gene;
  free(gene->sflen);
  free(gene->scafs[0]);
  free(gene->scafs);
  free(gene);
}

int64 Size_Of_Genome(Genome *_gene)
{ _Genome *gene = (_Genome *) _gene;
  int64 nbps, size;
  int   i;

  nbps = 0;
  for (i = 0; i < gene->sfnum; i++)
    nbps  += (gene->sflen[i]+3)/4;
  size = sizeof(Genome) + gene->sfnum*(sizeof(int64)+sizeof(uint8 *)) + nbps;
  return ((size-1)/0x100000+1);
}

void Print_Genome(Genome *_gene, FILE *file)
{ _Genome *gene = (_Genome *) _gene;
  int64 i, j, k, b;

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

Genome *Load_Genome(char *name)
{ _Genome   *db;
  Entry      entry;
  FILE      *gene;

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

    db    = Malloc(sizeof(_Genome),"Allocating genome");
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
        slen = j;
        sflen[sfnum++] = slen;

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

  return ((Genome *) db);
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
    _Genome *gene;
    int      max_blk;
  } _Haplotype;

void Free_Haplotype(Haplotype *_hap)
{ _Haplotype *hap = (_Haplotype *) _hap;
  free(hap->blocks[0]);
  free(hap->blocks);
  free(hap->snps);
  free(hap);
}

int64 Size_Of_Haplotype(Haplotype *_hap)
{ _Haplotype *hap = (_Haplotype *) _hap;
  int64 size;

  size = sizeof(Haplotype) + (hap->gene->sfnum+1)*sizeof(Block *)
       + ((hap->blocks[hap->gene->sfnum] - hap->blocks[0])+1) * sizeof(Block)
       + (hap->blocks[hap->gene->sfnum]->snp - hap->blocks[0]->snp) * sizeof(uint32);
  return ((size-1)/0x100000+1);
}

void Print_Haplotype(Haplotype *_hap, FILE *file)
{ _Haplotype *hap = (_Haplotype *) _hap;
  int     i;
  int64   off, beg, end;
  uint32 *s;
  Block  *b, *e;

  off = 0;
  fprintf(file,"\nHaplotype:\n");
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

Haplotype *Load_Haplotype(FILE *file, Genome *_gene)
{ _Genome    *gene = (_Genome *) _gene;
  Block      *blocks;
  Block     **bptrs;
  uint32     *snps;
  _Haplotype *hap;
  int64       numb, nums;

  Block  *cblk;
  uint32 *csnp;
  int64   ntot, mblk, off;

  int         line;
  char        type, name[3];
  int64       nbases, beg, end;
  int         nlen, nblock, nsnp, nsnp2, loc, snp, last;
  int         i, b, s;

  line = 1;
  if (fscanf(file,"1 %d ",&nlen) != 1)
    { fprintf(stderr,"%s: ONE-Code header missing?\n",Prog_Name);
      exit (1);
    }
  if (nlen != 3)
    { fprintf(stderr,"%s: ONE-Code type name length is not 3\n",Prog_Name);
      exit (1);
    }
  if (fscanf(file," %s\n",name) != 1)
    { fprintf(stderr,"%s: ONE-Code type name missing?\n",Prog_Name);
      exit (1);
    }
  if (strcmp(name,"hap") != 0)
    { fprintf(stderr,"%s: ONE-Code type name is not 'hap'\n",Prog_Name);
      exit (1);
    }

  nums = 0;
  numb = 1;
  for (i = 0; i < gene->sfnum; i++)
    { line += 1;
      if (fscanf(file,"c %d %lld\n",&nblock,&nbases) != 2)
        { fprintf(stderr,"%s: Expecting proper c-line\n",Prog_Name);
          exit (1);
        }
      numb += nblock;
      for (b = 0; b < nblock; b++)
        { line += 1;
          if (fscanf(file,"B %lld %lld\n",&beg,&end) != 2)
            { fprintf(stderr,"%s: Expecting proper B-line\n",Prog_Name);
              exit (1);
            }
          if (fscanf(file," %c",&type) != 1)
            break;
          if (type == 'P')
            { line += 1;
              if (fscanf(file," %d",&nsnp) != 1)
                { fprintf(stderr,"%s: Expecting length of P-line\n",Prog_Name);
                  exit (1);
                }
              nums += nsnp;
              last = -1;
              for (s = 0; s < nsnp; s++)
                { if (fscanf(file," %d",&loc) != 1)
                    { fprintf(stderr,"%s: Expecting a snp location\n",Prog_Name);
                      exit (1);
                    }
                  if (loc < 1 || end-beg < loc)
                    { fprintf(stderr,"%s: Snp location is not withing its block\n",Prog_Name);
                      exit (1);
                    }
                  if (loc <= last)
                    { fprintf(stderr,"%s: Snp location is not in order\n",Prog_Name);
                      exit (1);
                    }
                  last = loc;
                }
              line += 1;
              if (fscanf(file,"\nS %d",&nsnp2) != 1)
                { fprintf(stderr,"%s: Expecting start of S-line\n",Prog_Name);
                  exit (1);
                }
              if (nsnp2 != nsnp)
                { fprintf(stderr,"%s: P- and S-line list don't have the same length\n",Prog_Name);
                  exit (1);
                }
              for (s = 0; s < nsnp; s++)
                { if (fscanf(file," %d",&snp) != 1)
                    { fprintf(stderr,"%s: Expecting a snp offset\n",Prog_Name);
                      exit (1);
                    }
                  if (snp < 1 || snp > 3)
                    { fprintf(stderr,"%s: Snp offest is not in [1,3]\n",Prog_Name);
                      exit (1);
                    }
                }
              fscanf(file,"\n");
            }
          else
            ungetc(type,file);
        }
    }
 
  blocks = Malloc(sizeof(Block)*(numb+1),"Allocating Haplotype");
  bptrs  = Malloc(sizeof(Block *)*(gene->sfnum+1),"Allocating Haplotype");
  snps   = Malloc(sizeof(uint32)*nums,"Allocating Haplotype");
  hap    = Malloc(sizeof(_Haplotype),"Allocating Haplotype");
  if (hap == NULL || snps == NULL || bptrs == NULL || blocks == NULL)
    exit (1);

  rewind(file);

  cblk = blocks;
  csnp = snps;
  ntot = 0;
  mblk = 0;
  off  = 0;

  fscanf(file,"1 %d ",&nlen);
  fscanf(file," %s\n",name);
  for (i = 0; i < gene->sfnum; i++)
    { fscanf(file,"c %d %lld\n",&nblock,&nbases);
      bptrs[i] = cblk;
      for (b = 0; b < nblock; b++)
        { fscanf(file,"B %lld %lld\n",&beg,&end);
          cblk->beg = beg+off;
          cblk->snp = csnp;
          cblk->cum = ntot;
          ntot += end-beg;
          if (end-beg > mblk)
            mblk = end-beg;
          fscanf(file," %c",&type);
          if (type == 'P')
            { fscanf(file," %d",&nsnp);
              for (s = 0; s < nsnp; s++)
                { fscanf(file," %d",&loc);
                  csnp[s] = (loc-1) << 2;
                }
              fscanf(file,"\nS %d",&nsnp2);
              for (s = 0; s < nsnp; s++)
                { fscanf(file," %d",&snp);
                  csnp[s] |= snp;
                }
              fscanf(file,"\n");
              csnp += nsnp;
            }
          else
            ungetc(type,file);
          cblk += 1;
        }
      off += (4 - (gene->sflen[i] & 0x3)) & 0x3;
    }
  cblk->cum = ntot;
  cblk->snp = csnp;
  bptrs[gene->sfnum] = cblk;

  hap->snps    = snps;
  hap->blocks  = bptrs;
  hap->gene    = gene;
  hap->max_blk = mblk;

  return ((Haplotype *) hap);
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

Genome *Haplotype_Sequence(Haplotype *_hap)
{ _Haplotype *hap = (_Haplotype *) _hap;
  int64   sfnum  = hap->gene->sfnum; 
  Block **blocks = hap->blocks;
  uint8  *sbase  = hap->gene->scafs[0];

  _Genome *db;
  int64   *sflen;
  uint8  **scafs, *bases, *seq;
  int64    p, nbps;
  Block   *b;
  int      i, r, s, len;

#ifndef DEBUG_HAPLO_MAKE
  (void) show_block;
#endif

  nbps = 0;
  for (i = 0; i < sfnum; i++)
    nbps += ((blocks[i+1]->cum - blocks[i]->cum)+3)/4;

  db    = Malloc(sizeof(_Genome),"Allocating haplotype sequence");
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

  return ((Genome *) db);
}


/*******************************************************************************************
 *
 *  Routines to read and manipulate read ground truth
 *
 ********************************************************************************************/

typedef struct
  { int      hap;
    int      scaf;
    int      orient;
    int64    rbeg;
    int64    rend;
    char    *ops;
    int     *length;
  } Script;

typedef struct
  { int64    nreads;
    float   *rate;    // [0..nhaps)
    Script  *edit;    // [0..nreads)
  } _Read_Truth;

void Free_Read_Truth(Read_Truth *_truth)
{ _Read_Truth *truth = (_Read_Truth *) _truth;
  free(truth->edit[0].ops);
  free(truth->edit[0].length);
  free(truth->edit);
  free(truth->rate);
  free(truth);
}

Read_Truth *Load_Read_Truth(FILE *file)
{ _Read_Truth *truth;
  int          line, nlen, nhap, nops, nlens, hap;
  int         *lens;
  char         name[3], *ops;
  int64        nreads, totrds, totops;
  float        rate, *rates;
  Script      *reads;
  int          i, j, h;

  line = 1;
  if (fscanf(file,"1 %d ",&nlen) != 1)
    { fprintf(stderr,"%s: ONE-Code header missing?\n",Prog_Name);
      exit (1);
    }
  if (nlen != 3)
    { fprintf(stderr,"%s: ONE-Code type name length is not 3\n",Prog_Name);
      exit (1);
    }
  if (fscanf(file," %s\n",name) != 1)
    { fprintf(stderr,"%s: ONE-Code type name missing?\n",Prog_Name);
      exit (1);
    }
  if (strcmp(name,"err") != 0)
    { fprintf(stderr,"%s: ONE-Code type name is not 'err'\n",Prog_Name);
      exit (1);
    }

  totrds = 0;
  totops = 0;
  nhap = 1;
  while ( ! feof(file))
    { line += 1;
      if (fscanf(file,"h %lld %f\n",&nreads,&rate) != 2)
        { fprintf(stderr,"%s: Expecting proper h-line\n",Prog_Name);
          exit (1);
        }
      for (i = 0; i < nreads; i++)
        { if (fscanf(file,"S %d %*d %*d %*lld %*lld\n",&hap) != 1)
            { fprintf(stderr,"%s: Expecting proper S-line\n",Prog_Name);
              exit (1);
            }
          if (hap != nhap)
            { fprintf(stderr,"%s: Read is not in haplotype %d but rather %d\n",Prog_Name,nhap,hap);
              exit (1);
            }
          if (fscanf(file,"O %d ",&nops) != 1)
            { fprintf(stderr,"%s: Expecting proper O-line\n",Prog_Name);
              exit (1);
            }
          if (fscanf(file,"%*s\n") != 0)
            { fprintf(stderr,"%s: Expecting ops string\n",Prog_Name);
              exit (1);
            }
          if (fscanf(file,"L %d ",&nlens) != 1)
            { fprintf(stderr,"%s: Expecting proper L-line\n",Prog_Name);
              exit (1);
            }
          if (nlens != 2*nops+1)
            { fprintf(stderr,"%s: L list should be %d elements, not %d\n",Prog_Name,2*nops+1,nlens);
              exit (1);
            }
          for (j = 0; j < nlens; j++)
            if (fscanf(file," %*d") != 0)
          if (fscanf(file,"\n") != 0)
            { fprintf(stderr,"%s: L list is too long\n",Prog_Name);
              exit (1);
            }
          totops += nops;
        }
      totrds += nreads;
      nhap   += 1;
    }

  truth = Malloc(sizeof(_Read_Truth),"Allocating Read Truth");
  rates = Malloc(sizeof(float)*nhap,"Allocating Read Truth");
  reads = Malloc(sizeof(Script)*(totrds+1),"Allocating Read Truth");
  ops   = Malloc(totops,"Allocating Read Truth");
  lens  = Malloc(2*totops+nreads,"Allocating Read Truth");

  truth->rate   = rates;
  truth->edit   = reads;
  reads->ops    = ops;
  reads->length = lens;

  rewind(file);

  totrds = 0;
  totops = 0;
  for (h = 0; h < nhap; h++)
    { if (fscanf(file,"h %lld %f\n",&nreads,rates+h) != 2)
        { fprintf(stderr,"%s: Expecting proper h-line\n",Prog_Name);
          exit (1);
        }
      for (i = 0; i < nreads; i++)
        { Script *r;

          r = reads+totrds;
          r->ops = ops + totops;
          r->length = lens + (2*totops + totrds);

          fscanf(file,"S %d %d %d %lld %lld\n",&r->hap,&r->scaf,&r->orient,&r->rbeg,&r->rend);
          fscanf(file,"O %d ",&nops);
          fscanf(file,"%s\n",r->ops);
          fscanf(file,"L %d ",&nlens);
          for (j = 0; j < nlens; j++)
            fscanf(file," %d",r->length+j);
          fscanf(file,"\n");

          totops += nops;
          totrds += 1;
        }
    }
  reads[totrds].ops = ops + totops;
  reads[totrds].length = lens + (2*totops + totrds);

  return ((Read_Truth *) truth);
}


#ifdef XXX

char *Haplotype_read(Haplotype *hap, int64 beg, int64 end)
{
  l = hap->block[hap->gene->sfnum] - 1;
  f = hap->block[0]
  b = f + (1.*beg/l[1]->cum) * (l-f);
  e = f + (1.*end/l[1]->cum) * (l-f);
  while (b->cum > beg)
    b -= 1;
  while (beg < b[1].cum)
    b += 1;
  while (e->cum > end)
    e -= 1;
  while (end < e[1].cum)
    e += 1;
  
  if (b < e)
    { len = get_sequence(sbase,b->beg + (beg-b->cum),b[1].cum-beg,seq);
      for (b++; b < e; b++)
        len = get_sequence(sbase,b->beg,b[1].cum-b->cum,seq);
      len = get_sequence(sbase,b->beg,end-b0>cum,seq);
    }
  else
    len = get_sequence(sbase,b->beg + (beg-b->cum),end-beg,seq);
}

#endif
  

/*******************************************************************************************
 *
 *  Interpret the ploidy tree and build the given number of haplotypes, returning
 *    an array of pointers to them along with their number.
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

int main(int argc, char *argv[])
{ Genome     *gene;
  int         nhaps;
  FILE       *f;

  (void) complement;

  if (argc != 3)
    { fprintf(stderr,"Usage <genome.fast> <trush>\n");
      exit (1);
    }

  gene = Load_Genome(argv[1]);

  for (nhaps = 1; 1; nhaps++)
    { f = fopen(Catenate(argv[2],Numbered_Suffix(".",nhaps,"."),"hap",""),"r");
      if (f == NULL)
         break;
      fclose(f);
    }
  nhaps -= 1;

  { Haplotype *haps[nhaps];
    int        h;

    for (h = 0; h < nhaps; h++)
      { f = fopen(Catenate(argv[2],".hap",Numbered_Suffix("",h+1,""),""),"r");
        haps[h] = Load_Haplotype(f,gene);
        fclose(f);
      }

    for (h = 0; h < nhaps; h++)
      Free_Haplotype(haps[h]);
  }

  Free_Genome(gene);

  exit (0);
}
