/*******************************************************************************************
 *
 *  Lib_Sim:
 *     C-library for using HIsim ground truth information of the analysis of an assembly
 *
 *  Author:  Gene Myers
 *  Date  :  November 2021
 *
 ********************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <ctype.h>

#undef   DEBUG_HAPLO_MAKE
#undef   DEBUG_TABLE
#undef   DEBUG_SPLAY
#undef   DEBUG_CHAIN
#define  DEBUG_COMPOSE
#undef   DEBUG_SNPS
#define  DEBUG_SCRIPT

#include "gene_core.h"
#include "lib_sim.h"

#define WIDTH  100    //  Line width for DNA sequence output

#define MAX_GAP 20    //  Maximum gap allowed in a pseudo alignment

static char dna[4] = { 'a', 'c', 'g', 't' };

static char *fmer[256];   //  2-bit packed byte -> 4bp string

static void setup_fmer_table()
{ static char _fmer[1280];
  char *t;
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
  size = sizeof(_Genome) + gene->sfnum*(sizeof(int64)+sizeof(uint8 *)) + nbps;
  return ((size-1)/0x100000+1);
}

int Scaffold_Count(Genome *_gene)
{ return (((_Genome *) _gene)->sfnum); }

int64 Genome_Length(Genome *_gene)
{ return (((_Genome *) _gene)->nbase); }

int64 Scaffold_Length(Genome *_gene, int i)
{ return (((_Genome *) _gene)->sflen[i]); }

char *Scaffold_Sequence(Genome *_gene, int i)
{ _Genome *gene = (_Genome *) _gene;
  int64    k, b;
  char    *seq;

  setup_fmer_table();

  int64  len  = gene->sflen[i];
  uint8 *base = gene->scafs[i];

  seq = Malloc(len+1,"Allocating scaffold string\n");
  if (seq == NULL)
    return (NULL);

  k = b = 0;
  for (k = 0; k+3 < len; k += 4)
    strncpy(seq+k,fmer[base[b++]],4);
  if (k < len)
    strncpy(seq+k,fmer[base[b]],len-k);

  return (seq);
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
  } _HapTruth;

void Free_HapTruth(HapTruth *_hap)
{ _HapTruth *hap = (_HapTruth *) _hap;
  free(hap->blocks[0]);
  free(hap->blocks);
  free(hap->snps);
  free(hap);
}

int64 Size_Of_HapTruth(HapTruth *_hap)
{ _HapTruth *hap = (_HapTruth *) _hap;
  int64 size;

  size = sizeof(_HapTruth) + (hap->gene->sfnum+1)*sizeof(Block *)
       + ((hap->blocks[hap->gene->sfnum] - hap->blocks[0])+1) * sizeof(Block)
       + (hap->blocks[hap->gene->sfnum]->snp - hap->blocks[0]->snp) * sizeof(uint32);
  return ((size-1)/0x100000+1);
}

void Print_HapTruth(HapTruth *_hap, FILE *file)
{ _HapTruth *hap = (_HapTruth *) _hap;
  int     i;
  int64   off, beg, end;
  uint32 *s;
  Block  *b, *e;

  off = 0;
  fprintf(file,"\nHaplotype GT:\n");
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

HapTruth *Load_HapTruth(FILE *file, Genome *_gene)
{ _Genome    *gene = (_Genome *) _gene;
  Block      *blocks;
  Block     **bptrs;
  uint32     *snps;
  _HapTruth  *hap;
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
    { fprintf(stderr,"%s: Line %d: A ONE-Code header missing?\n",Prog_Name,line);
      exit (1);
    }
  if (nlen != 3)
    { fprintf(stderr,"%s: Line %d: ONE-Code type name length is not 3\n",Prog_Name,line);
      exit (1);
    }
  if (fscanf(file," %s\n",name) != 1)
    { fprintf(stderr,"%s: Line %d: ONE-Code type name missing?\n",Prog_Name,line);
      exit (1);
    }
  if (strcmp(name,"hap") != 0)
    { fprintf(stderr,"%s: Line %d: ONE-Code type name is not 'hap'\n",Prog_Name,line);
      exit (1);
    }

  nums = 0;
  numb = 1;
  for (i = 0; i < gene->sfnum; i++)
    { line += 1;
      if (fscanf(file,"c %d %lld\n",&nblock,&nbases) != 2)
        { fprintf(stderr,"%s: Line %d: Expecting proper c-line\n",Prog_Name,line);
          exit (1);
        }
      numb += nblock;
      for (b = 0; b < nblock; b++)
        { line += 1;
          if (fscanf(file,"B %lld %lld\n",&beg,&end) != 2)
            { fprintf(stderr,"%s: Line %d: Expecting proper B-line\n",Prog_Name,line);
              exit (1);
            }
          if (fscanf(file," %c",&type) != 1)
            break;
          if (type == 'P')
            { line += 1;
              if (fscanf(file," %d",&nsnp) != 1)
                { fprintf(stderr,"%s: Line %d: Expecting length of P-line\n",Prog_Name,line);
                  exit (1);
                }
              nums += nsnp;
              last = -1;
              for (s = 0; s < nsnp; s++)
                { if (fscanf(file," %d",&loc) != 1)
                    { fprintf(stderr,"%s: Line %d: Expecting a snp location\n",Prog_Name,line);
                      exit (1);
                    }
                  if (loc < 1 || end-beg < loc)
                    { fprintf(stderr,"%s: Line %d: Snp location is not withing its block\n",
                                     Prog_Name,line);
                      exit (1);
                    }
                  if (loc <= last)
                    { fprintf(stderr,"%s: Line %d: Snp location is not in order\n",Prog_Name,line);
                      exit (1);
                    }
                  last = loc;
                }
              line += 1;
              if (fscanf(file,"\nS %d",&nsnp2) != 1)
                { fprintf(stderr,"%s: Line %d: Expecting start of S-line\n",Prog_Name,line);
                  exit (1);
                }
              if (nsnp2 != nsnp)
                { fprintf(stderr,"%s: Line %d: P- and S-line list don't have the same length\n",
                                 Prog_Name,line);
                  exit (1);
                }
              for (s = 0; s < nsnp; s++)
                { if (fscanf(file," %d",&snp) != 1)
                    { fprintf(stderr,"%s: Line %d: Expecting a snp offset\n",Prog_Name,line);
                      exit (1);
                    }
                  if (snp < 1 || snp > 3)
                    { fprintf(stderr,"%s: Line %d: Snp offest is not in [1,3]\n",Prog_Name,line);
                      exit (1);
                    }
                }
              fscanf(file,"\n");
            }
          else
            ungetc(type,file);
        }
    }
 
  blocks = Malloc(sizeof(Block)*(numb+1),"Allocating Haplotype GT");
  bptrs  = Malloc(sizeof(Block *)*(gene->sfnum+1),"Allocating Haplotype GT");
  snps   = Malloc(sizeof(uint32)*nums,"Allocating Haplotype GT");
  hap    = Malloc(sizeof(_HapTruth),"Allocating Haplotype GT");
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
  cblk->beg = cblk[-1].beg + (ntot-cblk[-1].cum);
  cblk->snp = csnp;
  bptrs[gene->sfnum] = cblk;

  hap->snps    = snps;
  hap->blocks  = bptrs;
  hap->gene    = gene;
  hap->max_blk = mblk;

  return ((HapTruth *) hap);
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

static void mutate_part(uint8 *seq, Block *b, uint32 beg, uint32 end)
{ uint32 *snp;
  uint32 w, p, o;
  int    i, len;

  snp  = b->snp;
  len  = b[1].snp - snp;
  seq -= beg;

  for (i = 0; i < len; i++)
    { w = snp[i];
      p = (w >> 2);
      if (p < beg)
        continue;
      if (p >= end)
        break;
      o = (w & 0x3);
      seq[p] = (seq[p]+o) & 0x3;
    }
}

static void mutate_block(uint8 *seq, Block *b)
{ uint32 *snp;
  uint32 w, p, o;
  int    i, len;

  snp  = b->snp;
  len  = b[1].snp - snp;

  for (i = 0; i < len; i++)
    { w = snp[i];
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

Genome *Haplotype_Genome(HapTruth *_hap)
{ _HapTruth *hap = (_HapTruth *) _hap;
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

  db    = Malloc(sizeof(_Genome),"Allocating haplotype genome");
  sflen = Malloc(sizeof(int64)*sfnum,"Allocating haplotype genome");
  scafs = Malloc(sizeof(uint8 *)*sfnum,"Allocating haplotype genome");
  bases = Malloc(sizeof(uint8)*nbps,"Allocating haplotype genome");
  seq   = Malloc(hap->max_blk,"Allocating haplotype genome");

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
          mutate_block(seq,b);
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
          if (r == 6)
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
 *  Routines to input and manage read ground truth
 *
 ********************************************************************************************/

typedef struct
  { int      src;
    int      orient;
    int64    rbeg;
    int64    rend;
    char    *ops;
  } Script;

typedef struct
  { int        hidx;
    HapTruth  *hap;
    int        ctg;
  } Contig;

typedef struct
  { int64    nreads;
    float   *rate;    // [0..nhaps)
    Script  *edit;    // [0..nreads)
    Contig  *srcs;    // [0..nhaps*nscafs)
    int     *length;
  } _ReadTruth;

int64 Size_Of_ReadTruth(ReadTruth *_truth)
{ _ReadTruth *truth = (_ReadTruth *) _truth;
  int64   size, totops;
  int     maxctg;

  totops = truth->edit[truth->nreads].ops - truth->edit[0].ops;
  maxctg = truth->edit[truth->nreads-1].src;

  size = sizeof(_ReadTruth) + ((truth->srcs[maxctg].hap - truth->srcs[0].hap)+1) *sizeof(float)
       + (truth->nreads+1) * sizeof(Script) + (maxctg+1) * sizeof(Contig)
       + totops + (2*totops+truth->nreads) * sizeof(int);
  return ((size-1)/0x100000+1);
}

int Number_Of_Reads(ReadTruth *_truth)
{ return (((_ReadTruth *) _truth)->nreads); }

void Free_ReadTruth(ReadTruth *_truth)
{ _ReadTruth *truth = (_ReadTruth *) _truth;
  free(truth->length);
  free(truth->edit[0].ops);
  free(truth->edit);
  free(truth->rate);
  free(truth);
}

ReadTruth *Load_ReadTruth(FILE *file, int nhaps, HapTruth **haps)
{ _ReadTruth *truth;
  int     line, nlen, chap, nops, nlens, hap;
  int     lctg, ctg, nctg;;
  int    *lens, *lenp;
  char    name[3], *ops;
  int64   nreads, totrds, totops;
  float   rate, *rates;
  Script *reads;
  Contig *srces;
  int     i, j, h, c;

  line = 1;
  if (fscanf(file,"1 %d ",&nlen) != 1)
    { fprintf(stderr,"%s: Line %d: ONE-Code header missing?\n",Prog_Name,line);
      exit (1);
    }
  if (nlen != 3)
    { fprintf(stderr,"%s: Line %d: ONE-Code type name length is not 3\n",Prog_Name,line);
      exit (1);
    }
  if (fscanf(file," %s\n",name) != 1)
    { fprintf(stderr,"%s: Line %d: ONE-Code type name missing?\n",Prog_Name,line);
      exit (1);
    }
  if (strcmp(name,"err") != 0)
    { fprintf(stderr,"%s: Line %d: ONE-Code type name is not 'err'\n",Prog_Name,line);
      exit (1);
    }

  totrds = 0;
  totops = 0;
  chap   = 0;
  nctg   = 0;
  while (1)
    { chap += 1;
      c = fgetc(file);
      if (feof(file))
        { if (chap < nhaps)
            { fprintf(stderr,"%s: Line %d: File is missing %d haplotypes\n",
                             Prog_Name,line,nhaps-chap);
              exit (1);
            }
          break;
        }
      else
        ungetc(c,file);
      line += 1;
      if (fscanf(file,"h %lld %f\n",&nreads,&rate) != 2)
        { fprintf(stderr,"%s: Line %d: Expecting proper h-line\n",Prog_Name,line);
          exit (1);
        }
      if (chap > nhaps)
        { fprintf(stderr,"%s: Line %d: File refers to more than %d haplotypes\n",
                         Prog_Name,line,nhaps);
          exit (1);
        }
      lctg = 0;
      for (i = 0; i < nreads; i++)
        { if (fscanf(file,"S %d %d %*d %*d %*d\n",&hap,&ctg) != 2)
            { fprintf(stderr,"%s: Line %d: Expecting proper S-line %d\n",Prog_Name,line,i);
              exit (1);
            }
          if (hap != chap)
            { fprintf(stderr,"%s: Line %d: Read is not in haplotype %d but rather %d\n",
                             Prog_Name,line,chap,hap);
              exit (1);
            }
          if (ctg != lctg)
            nctg += 1;
          if (fscanf(file,"O %d ",&nops) != 1)
            { fprintf(stderr,"%s: Line %d: Expecting proper O-line\n",Prog_Name,line);
              exit (1);
            }
          if (nops > 0)
            { if (fscanf(file,"%*s\n") != 0)
                { fprintf(stderr,"%s: Line %d: Expecting ops string\n",Prog_Name,line);
                  exit (1);
                }
            }
          if (fscanf(file,"L %d ",&nlens) != 1)
            { fprintf(stderr,"%s: Line %d: Expecting proper L-line %d\n",Prog_Name,line,i);
              exit (1);
            }
          if (nlens != 2*nops+1)
            { fprintf(stderr,"%s: Line %d: L list should be %d elements, not %d\n",
                             Prog_Name,line,2*nops+1,nlens);
              exit (1);
            }
          for (j = 0; j < nlens; j++)
            if (fscanf(file," %*d") != 0)
              { fprintf(stderr,"%s: Line %d: L list isn't long enough (%d)\n",Prog_Name,line,j);
                exit (1);
              }
          if (fscanf(file,"\n") != 0)
            { fprintf(stderr,"%s: Line %d: L list is too long\n",Prog_Name,line);
              exit (1);
            }
          totops += nops;
        }
      nctg   += 1;
      totrds += nreads;
    }

  truth = Malloc(sizeof(_ReadTruth),"Allocating Read Truth");
  rates = Malloc(sizeof(float)*nhaps,"Allocating Read Truth");
  reads = Malloc(sizeof(Script)*(totrds+1),"Allocating Read Truth");
  srces = Malloc(sizeof(Contig)*nctg,"Allocating Read Truth");
  ops   = Malloc(totops,"Allocating Read Truth");
  lens  = Malloc((2*totops+totrds)*sizeof(int),"Allocating Read Truth");

  truth->nreads = totrds;
  truth->rate   = rates;
  truth->edit   = reads;
  truth->srcs   = srces;
  truth->length = lens;
  reads->ops    = ops;

  rewind(file);
  fscanf(file,"1 %*d %*s\n");

  totrds = 0;
  totops = 0;
  nctg   = 0;
  chap   = 0;
  for (h = 0; h < nhaps; h++)
    { fscanf(file,"h %lld %f\n",&nreads,rates+h);
      lctg = 0;
      for (i = 0; i < nreads; i++)
        { Script *r;

          r = reads+totrds;
          r->ops = ops + totops;
          lenp   = lens + (2*totops + totrds);

          fscanf(file,"S %d %d",&hap,&ctg);
          if (ctg != lctg)
            { nctg += 1;
              lctg = ctg;
              srces[nctg].hidx = chap+1;
              srces[nctg].hap  = haps[chap];
              srces[nctg].ctg  = ctg;
            }
          r->src = nctg;
          fscanf(file," %d %lld %lld\n",&r->orient,&r->rbeg,&r->rend);
          fscanf(file,"O %d ",&nops);
          if (nops > 0)
            fscanf(file,"%s\n",r->ops);
          fscanf(file,"L %d ",&nlens);
          for (j = 0; j < nlens; j++)
            fscanf(file," %d",lenp+j);
          fscanf(file,"\n");

          totops += nops;
          totrds += 1;
        }
      nctg += 1;
      chap += 1;
    }
  reads[totrds].ops = ops + totops;

  return ((ReadTruth *) truth);
}


/*******************************************************************************************
 *
 *  Routines to access read source ground truth
 *
 ********************************************************************************************/

Source *Get_Source(ReadTruth *r, int64 i, Source *src)
{ _ReadTruth  *truth = (_ReadTruth *) r;
  Script  *e     = truth->edit + i;
  Contig  *s     = truth->srcs + e->src;

  if (src == NULL)
    { src = Malloc(sizeof(Source),"Allocating read source");
      if (src == NULL)
        exit (1);
    }

  src->hap    = s->hap;
  src->hidx   = s->hidx;
  src->contig = s->ctg;
  src->orient = e->orient;
  src->beg    = e->rbeg;
  src->end    = e->rend;

  return (src);
}

void Free_Read_Source(Source *source)
{ free(source); }


/*******************************************************************************************
 *
 *  Routines to access read haplotype block ground truth
 *
 ********************************************************************************************/

typedef struct
  { int64      beg;    //  the came from [beg,end] of haplotype 'hap'
    int64      end;
    Block     *fst;    //  the read is the cat of blocks [fst,lst] where the 1st and last
    Block     *lst;    //     are possibly a suffix/prefix thereof
    HapTruth  *hap;
    int64      gbeg;   //  almost all blocks in source genome coords are in the interval
    int64      gend;   //     [gbeg,gend], those that are not total govf base pairs
    int64      govf;
  } _Slice;

Slice *Get_Slice(Source *source, Slice *_slice)
{ _HapTruth  *hap;
  int64       beg, end;
  Block      *s, *f;
  Block      *b, *e;
  _Slice     *slice;
  int         N;

  hap = (_HapTruth  *) (source->hap);
  beg = source->beg + hap->blocks[source->contig-1]->cum;
  end = source->end + hap->blocks[source->contig-1]->cum;
  s = hap->blocks[hap->gene->sfnum] - 1;
  f = hap->blocks[0];
  b = f + (beg * (s-f))/s[1].cum;
  e = f + (end * (s-f))/s[1].cum;
  while (beg < b->cum)
    b -= 1;
  while (beg >= b[1].cum)
    b += 1;
  while (end <= e->cum)
    e -= 1;
  while (end > e[1].cum)
    e += 1;

  if (_slice == NULL)
    { slice = Malloc(sizeof(_Slice),"Allocating Haplotype Slice");
      if (slice == NULL)
        exit (1);
    }
  else
    slice = (_Slice *) _slice;

  slice->beg = beg;
  slice->end = end;
  slice->fst = b;
  slice->lst = e;
  slice->hap = hap;

  //  Compute l.i.s. of slice in genome coordinate system

  N = (e-b)+1;

  { int   P[N], M[N+1];
    int64 S[N+1];
    int    i, n, L;
    int    bot, top;

    //  lis computation proper

    L = 0;
    M[0] = -1;
    S[0] = 0;
    for (i = 0; i < N; i++)
      { int   l, h, m;
        int64 x;

        x = b[i].beg;

        l = 1;
        h = L + 1;
        while (l < h)
          { m = (l+h) >> 1;
            if (b[M[m]].beg < x)
              l = m+1;
            else
              h = m;
          }

        if (l > L || S[l-1] + (b[i+1].cum - b[i].cum) > S[l])
          { P[i] = M[l-1];
            M[l] = i;
            S[l] = S[l-1] + (b[i+1].cum - b[i].cum);
            if (l > L)
              L = l;
          }
      }

    //  build inverse of linked list P in M

    i = top = M[L];
    while ((n = P[i]) >= 0)
      { M[n] = i;
        i = n;
      }
    bot = i;

    //  trim prefix and suffix of lis that have low density (gap/length >= 5)

    while (top > bot)
      { int64 len, bp;

        s = b + top;
        f = b + P[top]; 
        if (s == e)
          len = end-e->cum;
        else
          len = s[1].cum - s->cum;
        bp = f->beg + (f[1].cum - f->cum);
        if ((s->beg - bp) / len < 5)
          break;
        top = P[top];
      }

    while (bot < top)
      { int64 len, ep;

        s = b + bot;
        f = b + M[bot]; 
        if (s == b)
          len = b[1].cum - beg;
        else
          len = s[1].cum - s->cum;
        ep = s->beg + (s[1].cum - s->cum);
        if ((f->beg - ep) / len < 5)
          break;
        bot = M[bot];
      }

    //  accumulate total bp not in lis and record in ->govf field

    { int64 excess;

      P[N]   = top;
      P[bot] = -1;
      excess = 0;
      for (i = N; i >= 0; i = P[i])
        { for (n = P[i]+1; n < i; n++)
            excess += b[n+1].cum - b[n].cum;
        }
      if (top != N-1)
        excess -= (e[1].cum - end); 
      if (bot != 0)
        excess -= (beg - b->cum);

      slice->govf = excess;
    }

    //  store range of trimmed lis

    s = b+bot;
    f = b+top;
    if (s == b)
      slice->gbeg = s->beg + (beg-s->cum);
    else
      slice->gbeg = s->beg;
    if (f == e)
      slice->gend = f->beg + (end-f->cum);
    else
      slice->gend = f->beg + (f[1].cum - f->cum);
    for (i = P[top]; i >= 0; i = P[i])
      { int64 x = b[i].beg + (b[i+1].cum - b[i].cum);
        if (x > slice->gend)
          slice->gend = x;
      }
  }

  return ((Slice *) slice);
}

void Free_Slice(Slice *slice)
{ free(slice); };

static void print_snps(Block *b, uint32 beg, uint32 end, FILE *file)
{ uint32 *snp;
  uint32 w, p;
  int    i, len;

  snp  = b->snp;
  len  = b[1].snp - snp;

  for (i = 0; i < len; i++)
    { w = snp[i];
      p = (w >> 2);
      if (p < beg || p >= end)
        continue;
      fprintf(file," %d(%d)",p,(w&0x3));
    }
  fprintf(file,"\n");
}

void Print_Slice(Slice *slice, int show_snps, FILE *file)
{ _Slice *s = (_Slice *) slice;
  Block *b, *e;

  b = s->fst;
  e = s->lst; 
  if (show_snps)
    { fprintf(file,"  LIS: <%lld,%lld:%lld>\n",s->gbeg,s->gend,s->govf);
      if (b < e)
        { fprintf(file,"  [%lld,%lld]:",b->beg + (s->beg-b->cum), b->beg + (b[1].cum-b->cum));
          print_snps(b,s->beg-b->cum,b[1].cum-b->cum,file);
          for (b++; b < e; b++)
            { fprintf(file,"  [%lld,%lld]:",b->beg, b->beg + (b[1].cum-b->cum));
              print_snps(b,0,b[1].cum-b->cum,file);
            }
          fprintf(file,"  [%lld,%lld]",e->beg, e->beg + (s->end-e->cum));
          print_snps(e,0,s->end-b->cum,file);
        }
      else
        { fprintf(file,"[%lld,%lld]",b->beg + (s->beg-b->cum), e->beg + (s->end-e->cum));
          print_snps(e,s->beg-b->cum,s->end-b->cum,file);
        }
    }
  else
    { if (b < e)
        { fprintf(file,"[%lld,%lld]",b->beg + (s->beg-b->cum), b->beg + (b[1].cum-b->cum));
          for (b++; b < e; b++)
            fprintf(file," [%lld,%lld]",b->beg, b->beg + (b[1].cum-b->cum));
          fprintf(file," [%lld,%lld]",e->beg, e->beg + (s->end-e->cum));
        }
      else
        fprintf(file,"[%lld,%lld]",b->beg + (s->beg-b->cum), e->beg + (s->end-e->cum));
      fprintf(file," :: <%lld,%lld:%lld>",s->gbeg,s->gend,s->govf);
    }
}

int Slice_Length(Slice *slice)
{ _Slice *s = (_Slice *) slice;

  return (s->end-s->beg);
}

int Snps_In_Slice(Slice *slice)
{ _Slice *s = (_Slice *) slice;
  uint32  *f, *g;
  uint32   c, p;

  f = s->fst->snp;
  c = s->beg - s->fst->cum;
  while (1)
    { if (f >= s->fst[1].snp)
        break;
      p = (*f >> 2);
      if (p >= c)
        break;
      f += 1;
    }

  g = s->lst->snp;
  c = s->end - s->lst->cum;
  while (1)
    { if (g >= s->lst[1].snp)
        break;
      p = (*g >> 2);
      if (p > c)
        break;
      g += 1;
    }

  return (g-f);
}

uint8 *Slice_Sequence(Slice *slice, uint8 *seq)
{ _Slice *s = (_Slice *) slice;
  uint8   *sbase  = ((_HapTruth *) s->hap)->gene->scafs[0];
  int      len, bln;
  Block   *b, *e;

  if (seq == NULL)
    { len = Slice_Length(slice);
      seq = Malloc(len,"Allocating Sequence");
      if (seq == NULL)
        exit (1);
    }

  b = s->fst;
  e = s->lst; 
  if (b < e)
    { len = get_sequence(sbase,b->beg + (s->beg-b->cum),b[1].cum-s->beg,seq);
      mutate_part(seq,b,s->beg-b->cum,b[1].cum-b->cum);
      for (b++; b < e; b++)
        { bln = b[1].cum-b->cum;
          get_sequence(sbase,b->beg,bln,seq+len);
          mutate_block(seq+len,b);
          len += bln;
        }
      get_sequence(sbase,e->beg,s->end-e->cum,seq+len);
      mutate_part(seq+len,e,0,s->end-e->cum);
    }
  else
    { get_sequence(sbase,b->beg + (s->beg-b->cum),s->end-s->beg,seq);
      mutate_part(seq,b,s->beg-b->cum,s->end-e->cum);
    }

  return (seq);
}

uint8 *True_Sequence(ReadTruth *r, int64 i, uint8 *seq)
{ Source src;
  _Slice slc;

  Slice_Sequence(Get_Slice(Get_Source(r,i,&src),(Slice *) &slc),seq);
  if (src.orient)
    complement(src.end-src.beg,seq);
  return (seq);
}


/*******************************************************************************************
 *
 *  Routines to access read error edits ground truth
 *
 ********************************************************************************************/

typedef struct
  { int   len;
    char *ops;
    int  *vals;
  } _Edit;

Edit *Get_Edit(ReadTruth *reads, int64 i, Edit *edit)
{ _ReadTruth *r = (_ReadTruth *) reads; 
  _Edit *e;

  if (edit == NULL)
    { e = Malloc(sizeof(_Edit),"Allocating Read Script");
      if (e == NULL)
        exit (1);
    }
  else
    e = (_Edit *) edit;
  e->ops  = r->edit[i].ops;
  e->len  = r->edit[i+1].ops - e->ops;
  e->vals = r->length + (2*(e->ops - r->edit[0].ops) + i);
  return ((Edit *) e);
}

void Free_Edit(Edit *edit)
{ free((_Edit *) edit); }

void Print_Edit(Edit *edit, FILE *file)
{ static char base_pair[] = { 'a', 'c', 'g', 't' };

  char *ops = ((_Edit *) edit)->ops;
  int  *val = ((_Edit *) edit)->vals;
  int   len = ((_Edit *) edit)->len;
  int   i, j;

  fprintf(file,"%d",val[0]);
  for (i = 0, j = 1; i < len; i++, j+=2)
    { fprintf(file," %c",ops[i]);
      switch (ops[i])
      { case 'h': case 'z': case 't':
        case 'H': case 'Z': case 'T':
        case 'd': case 'i': case 's':
          fprintf(file,"%d",val[j]);
          break;
        case 'I':
        case 'S':
          fprintf(file,"%c",base_pair[val[j]]);
          break;
        default:
          break;
      }
      fprintf(file," %d",val[j+1]);
    }
}

int Edit_Length(Edit *edit)
{ char *ops = ((_Edit *) edit)->ops;
  int  *val = ((_Edit *) edit)->vals;
  int   len = ((_Edit *) edit)->len;
  int   i, j, olen;

  olen = val[0];
  for (i = 0, j = 1; i < len; i++, j+=2)
    { switch (ops[i])
      { case 'H': case 'Z': case 'T':
        case 'i': case 's':
          olen += val[j];
          break;
        case 'I':
        case 'S':
        case 'x':
          olen += 1;
          break;
        default:
          break;
      }
      olen += val[j+1];
    }

  return (olen);
}

int Errors_In_Edit(Edit *edit)
{ char *ops = ((_Edit *) edit)->ops;
  int  *val = ((_Edit *) edit)->vals;
  int   len = ((_Edit *) edit)->len;
  int   i, j, sum;

  sum = 0;
  for (i = 0, j = 1; i < len; i++, j+=2)
    { switch (ops[i])
      { case 'h': case 'z': case 't':
        case 'H': case 'Z': case 'T':
        case 'd': case 'i': case 's':
          sum += val[j];
          break;
        default:
          sum += 1;
          break;
      }
    }

  return (sum);
}

int Error_Type_In_Edit(Edit *edit, char kind)
{ char *ops = ((_Edit *) edit)->ops;
  int  *val = ((_Edit *) edit)->vals;
  int   len = ((_Edit *) edit)->len;
  int   i, j, sum;

  sum = 0;
  for (i = 0, j = 1; i < len; i++, j+=2)
    { if (ops[i] != kind)
        continue;
      switch (ops[i])
      { case 'h': case 'z': case 't':
        case 'H': case 'Z': case 'T':
        case 'd': case 'i': case 's':
          sum += val[j];
          break;
        default:
          sum += 1;
          break;
      }
    }

  return (sum);
}

static int Period[128] =
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0
  };

uint8  *Edit_Sequence(Edit *edit, uint8 *in, uint8 *out)
{ char *ops = ((_Edit *) edit)->ops;
  int  *val = ((_Edit *) edit)->vals;
  int   len = ((_Edit *) edit)->len;
  int   i, j, p, v, o, olen;

  if (out == NULL)
    { len = Edit_Length(edit);
      out = Malloc(len,"Allocating Sequence");
      if (out == NULL)
        exit (1);
    }

  v = val[0];
  if (v > 0)
    memcpy(out,in,v);
  in  += v;
  olen = v;
  for (i = 0, j = 1; i < len; i++, j+=2)
    { v = val[j];
      o = ops[i];
      switch (o)
      { case 'h': case 'z': case 't':
          in += v;
          break;
        case 'H': case 'Z': case 'T':
          p = Period[o];
          while (v > p)
            { memcpy(out+olen,in-p,p);
              olen += p;
              v  -= p;
            }
          memcpy(out+olen,in-p,v);
          olen += v;
          break;
        case 'I':
          out[olen++] = v;
          break;
        case 'S':
          out[olen++] = v;
          in += 1;
          break;
        case 'D':
          in += 1;
          break;
        case 'd':
        case 's':
        case 'i':
          for (i = 0; i < v; i++)
            out[olen++] = o; 
          break;
        case 'x':
          out[olen++] = 'x';
          in += 1;
          break;
      }

      v = val[j+1];
      if (v > 0)
        { memcpy(out+olen,in,v);
          in   += v;
          olen += v;
        }
    }

  return (out);
}

uint8 *Read_Sequence(ReadTruth *r, int64 i, uint8 *in, uint8 *out)
{ _Edit *edit;
  uint8 *true, *read;

  true = True_Sequence(r,i,in);
  read = Edit_Sequence(Get_Edit(r,i,(Edit *) (&edit)),in,out);
  if (in == NULL)
    free(true);
  return (read);
}

int Print_Fasta(ReadTruth *_reads, FILE *file)
{ _ReadTruth *truth = (_ReadTruth *) _reads;
  Source src;
  _Slice _slice;
  _Edit  _edit;
  Slice  *slice = (Slice *) (&_slice);
  Edit   *edit  = (Edit  *) (&_edit);
  int     tlen, elen;
  int     tmax, emax;
  uint8  *tseq, *eseq;
  char   *rseq;
  int     i, j, nreads;

  tmax = emax = 0;
  tlen = elen = 0;
  tseq = eseq = NULL;
  rseq = NULL;

  nreads = truth->nreads;
  for (i = 0; i < nreads; i++)
    { Get_Source(_reads,i,&src);
      Get_Slice(&src,slice);
      Get_Edit(_reads,i,edit);

      tlen = Slice_Length(slice);
      elen = Edit_Length(edit);

      if (tlen > tmax)
        { tmax = 1.2*tlen+1000;
          tseq = Realloc(tseq,tmax+1,"Reallocating true sequence buffer");
          if (tseq == NULL)
            return (1);
        }

      if (elen > emax)
        { emax = 1.2*elen+1000;
          eseq = Realloc(eseq,emax+1,"Reallocating true sequence buffer");
          if (eseq == NULL)
            return (1);
          rseq = (char *) eseq;
        }

      Slice_Sequence(slice,tseq);

      if (src.orient)
        complement(tlen,tseq);

      Edit_Sequence(edit,tseq,eseq);

      for (j = 0; j < elen; j++)
        rseq[j] = dna[eseq[j]];
      rseq[elen] = '\0';

      fprintf(file,">Sim %d %d %c %lld %lld\n",src.hidx,src.contig,
                                               src.orient?'-':'+',src.beg,src.end);
      for (j = 0; j+WIDTH < elen; j += WIDTH)
        fprintf(file,"%.*s\n",WIDTH,rseq+j);
      if (j < elen)
        fprintf(file,"%s\n",rseq+j);
    }

  return (1);
}


/*******************************************************************************************
 *
 *  Splay tree routines for ordered list: INSERT, DELETE, FIND, NEXT
 *
 ********************************************************************************************/

typedef struct vtx
  { struct vtx *L, *R;
    int64       V;
    int         score;
    int         link;
  } NODE;

#ifdef DEBUG_CHAIN

static void PRINT_LIST(NODE *v)
{ if (v == NULL)
    return;
  PRINT_LIST(v->L);
  printf(" %lld:%d:%d",v->V,v->score,v->link);
  PRINT_LIST(v->R);
}

#endif

#ifdef DEBUG_SPLAY

static void PRINT_TREE(NODE *v, int deep, NODE *space)
{ if (v == NULL)
    return;
  PRINT_TREE(v->R,deep+3,space);
  printf("%*s %lld:%d:%d (%ld)\n",deep,"",v->V,v->score,v->link,v-space);
  PRINT_TREE(v->L,deep+3,space);
}

#endif

static NODE *SPLAY(NODE *v, int64 x)    //  Assumes x is in the tree
{ NODE *u, *n;

  if (v == NULL || x == v->V)
    return (v);
  if (x < v->V)
    { u = v->L;
      if (x == u->V)
        { v->L = u->R;
          u->R = v;
          return (u);
        }
      if (x < u->V)
        { n = SPLAY(u->L,x);
          v->L = u->R;
          u->R = v;
          u->L = n->R;
          n->R = u;
        }
      else
        { n = SPLAY(u->R,x);
          v->L = n->R;
          u->R = n->L;
          n->L = u;
          n->R = v;
        }
    }
  else
    { u = v->R;
      if (x == u->V)
        { v->R = u->L;
          u->L = v;
          return (u);
        }
      if (x > u->V)
        { n = SPLAY(u->R,x);
          v->R = u->L;
          u->L = v;
          u->R = n->L;
          n->L = u;
        }
      else
        { n = SPLAY(u->L,x);
          v->R = n->L;
          u->L = n->R;
          n->R = u;
          n->L = v;
        }
    }
  return (n); 
}

static NODE *FIND(NODE *v, int64 x)   //  Find v s.t. v->V <= x && x < v->next->V
{ NODE *u;

  if (v == NULL || v->V == x)
    return (v);
  if (x < v->V)
    return (FIND(v->L,x));
  else
    { u = FIND(v->R,x);
      if (u == NULL)
        return (v);
      else
        return (u);
    }
}

static NODE *NEXT(NODE *v, NODE *t, NODE *w)
{ if (v == NULL || t->V == v->V)
    { if (v->R != NULL)
        { w = v->R;
          while (w->L != NULL)
            w = w->L;
        }
      return (w);
    }
  if (t->V < v->V)
    return (NEXT(v->L,t,v));
  else
    return (NEXT(v->R,t,w));
}

static NODE *JOIN(NODE *v, NODE *w)
{ NODE *p;

  if (v == NULL)
    return (w);
  for (p = v; p->R != NULL; p = p->R)
    ;
  v = SPLAY(v,p->V);
  v->R = w;
  return (v);
}

static NODE *INSERT(NODE *v, NODE *new)
{ NODE *u, *p;

  if (v == NULL)
    return (new);
  u = FIND(v,new->V);
  if (u != NULL && u->R == NULL)
    u->R = new;
  else
    { if (u == NULL)
        p = v;
      else  // u->R == NULL
        p = u->R;
      while (p->L != NULL)
        p = p->L;
      p->L = new;
    }
  return (SPLAY(v,new->V));
}

static NODE *DELETE(NODE *v, NODE *old)
{ NODE *u, *w;

  u = FIND(v,old->V);
  if (u == NULL || u->V != old->V)
    return (NULL);
  v = SPLAY(v,old->V);
  w = JOIN(v->L,v->R);
  return (w);
}


/*******************************************************************************************
 *
 *  Auxiliary routines to compose blocks, iterate over SNPs, and map through edits
 *
 ********************************************************************************************/

typedef struct     //  End point of a match interval
  { int  c1;
    int  c2;
  } Point;

typedef struct     //  SNP iterator
  { Block  *blk;
    uint32 *cur;
    int64   pos;
  } SNPit;

static void Init_SNP(SNPit *it, Block *blk)
{ it->blk = blk;
  it->cur = blk->snp;
}

  //  Deliver SNPs in haplo-interval [beg,end] one at a time in order return 1 when
  //    there are no more

static int Next_SNP(SNPit *it, int beg, int end)
{ Block   *b;
  uint32  *f, *g;
  int64    p;

  b = it->blk;
  f = it->cur;
  g = b[1].snp;
  while (1)
    { if (f >= g)
        { b += 1;
          g = b[1].snp;
        }
      p = b->cum + (int) (*f >> 2);
      if (p >= end)
        { it->blk = b;
          it->cur = f;
          it->pos = -1;
          return (0);
        }
      if (p >= beg)
        { it->blk = b;
          it->cur = f+1;
          it->pos = p; 
          return (1);
        }
      f += 1;
    }

  return (1);
}

typedef struct          //  Edit script mapper (from perfect to error or vice versa)
  { int   len;
    char *ops;
    int  *val;
    int   i, j;
    int   ilen, olen;
    int   ipre, opre;
    int   iedl, oedl;
  } Map;

static inline void Init_Map(_Edit *edit, Map *map)
{ map->len  = edit->len;
  map->ops  = edit->ops;
  map->val  = edit->vals;
  map->i    = 0;
  map->j    = 1;
  map->ilen = edit->vals[0];
  map->olen = edit->vals[0];
  map->ipre = 0;
  map->opre = 0;
  map->iedl = 0;
  map->oedl = 0;
}

  //   Map pnt from error space to perfect space.
  //      Only called on the beg coord of blocks that are a refinement of the edit
  //      so mapped point is never in the middle of an edit.

static inline int Inverse_Map(Map *map, int pnt)
{
  while (pnt >= map->olen)
    { switch (map->ops[map->i++])
      { case 'h': case 'z': case 't':
          map->ilen += map->val[map->j];
          break;
        case 'H': case 'Z': case 'T':
          map->olen += map->val[map->j];
          break;
        case 'I':
          map->olen += 1;
          break;
        case 'S':
          map->ilen += 1;
          map->olen += 1;
          break;
        case 'D':
          map->ilen += 1;
          break;
      }
      map->j += 1;
      map->olen += map->val[map->j];
      map->ilen += map->val[map->j];
      map->j += 1;
    }
  return (map->ilen - (map->olen - pnt));
}

  //  Map pnt from perfect space to error space being careful to correctly map
  //    points when the edit script deletes a prefix (first = 1) or suffix (first = 0)
  //    of the perfect read interval.

static inline int Boundary_Map(Map *map, int pnt, int first)
{
  while (pnt >= map->ilen && map->i < map->len)
    { map->ipre = map->ilen;
      map->opre = map->olen;
      switch (map->ops[map->i++])
      { case 'h': case 'z': case 't':
          map->ilen += map->val[map->j];
          break;
        case 'H': case 'Z': case 'T':
          map->olen += map->val[map->j];
          break;
        case 'I':
          map->olen += 1;
          break;
        case 'S':
          map->ilen += 1;
          map->olen += 1;
          break;
        case 'D':
          map->ilen += 1;
          break;
      }
      map->j += 1;
      map->iedl = map->ilen;
      map->oedl = map->olen;
      map->olen += map->val[map->j];
      map->ilen += map->val[map->j];
      map->j += 1;
    }
  if (first)
    { if (pnt >= map->iedl)
        return (map->oedl + (pnt-map->iedl));
      else if (map->iedl-pnt > map->oedl-map->opre)
        return (map->opre);
      else
        return (map->oedl - (map->iedl-pnt));
    }
  else
    { if (pnt > map->iedl)
        return (map->oedl + (pnt-map->iedl));
      else if (pnt-map->ipre > map->oedl-map->opre)
        return (map->oedl);
      else
        return (map->opre + (pnt-map->ipre));
    }
}

static void Print_BA(int no, Point *pairs, FILE *file)
{ int i;

  for (i = 0; i < no; i += 2)
    { fprintf(file,"   [%5d,%5d] x",pairs[i].c1,pairs[i+1].c1);
      fprintf(file," [%5d,%5d]\n",pairs[i].c2,pairs[i+1].c2);
    }
}

  //  Convert an edit script into a list of matching blocks in epair returning its length.
  //    Epair is guaranteed to be big enough.  If flip then from error to perfect map, otherwise
  //    perfect to error.

static int Edits_2_Blocks(_Edit *edit, int beg, int end, Point *epair, int flip)
{ char *ops = edit->ops;
  int  *val = edit->vals;
  int   len = edit->len;
  int   i, j, v, olen, ilen;
  int   npts;

#ifdef DEBUG_COMPOSE
  printf("Read Edit: ");
  Print_Edit((Edit *) edit,stdout);
  printf("\n");
#endif

  npts = 0;
  ilen = olen = 0;
  i = j = -1;
  while (1)
    { v = val[j+1];
      if (v > 0)
        { if (ilen >= end)
            break;
          if (beg < ilen+v)
            { if (end < ilen+v)
                v = end-ilen;
              if (beg > ilen)
                { olen += (beg-ilen);
                  v    -= (beg-ilen);
                  ilen  = beg;
                }
              if (v > 0)
                { epair[npts].c1 = ilen;
                  epair[npts].c2 = olen;
                  npts += 1;
                  ilen  += v;
                  olen  += v;
                  epair[npts].c1 = ilen;
                  epair[npts].c2 = olen;
                  npts += 1;
                }
            }
          else
            { ilen += v;
              olen += v;
            }
        }

      i += 1;
      j += 2;
      if (i >= len)
        break;

      switch (ops[i])
      { case 'h': case 'z': case 't':
          ilen += val[j];
          break;
        case 'H': case 'Z': case 'T':
          olen += val[j];
          break;
        case 'I':
          olen += 1;;
          break;
        case 'S':
          olen += 1;
          ilen += 1;
          break;
        case 'D':
          ilen += 1;
          break;
      }
    }

  if (flip)
    for (i = 0; i < npts; i++)
      { v = epair[i].c1;
        epair[i].c1 = epair[i].c2;
        epair[i].c2 = v;
      }

#ifdef DEBUG_COMPOSE
  Print_BA(npts,epair,stdout);
#endif

  return (npts);
}

  //  Compose two sets of matching blocks along the common dimension placing it
  //    in epair and returning its length.  Epair is guaranteed to be large enough.

static int Compose_Blocks(int n1, Point *e1, int n2, Point *e2, Point *epair)
{ int b, e;
  int i, j;
  int bi, bj;
  int ei, ej;
  int npts;

  npts = 0;
  i = j = 0;
  while (i < n1 && j < n2) 
    { bi = e1[i].c2;
      bj = e2[j].c1;
      ei = e1[i+1].c2;
      ej = e2[j+1].c1;
      if (ei <= bj)
        i += 1;
      else if (ej <= bi)
        j += 1;
      else
        { if (bi < bj)
            b = bj;
          else
            b = bi;
          if (ei < ej)
            e = ei;
          else
            e = ej;
          if (e > b)
            { epair[npts].c1 = e1[i].c1 + (b-bi);
              epair[npts].c2 = e2[j].c2 + (b-bj);
              npts += 1;
              epair[npts].c1 = e1[i].c1 + (e-bi);
              epair[npts].c2 = e2[j].c2 + (e-bj);
              npts += 1;
            }
          if (ei == e)
            i += 2;
          if (ej == e)
            j += 2;
        }
    }

#ifdef DEBUG_COMPOSE
  printf("Composition:\n");
  Print_BA(npts,epair,stdout);
#endif

  return (npts);
}

static int Block_Diff(int npts, Point *pair, int b1, int b2, int e1, int e2)
{ int   diff;
  int   l1, l2;
  int   d1, d2;
  int   i;
  
  diff = 0;
  l1 = b1;
  l2 = b2;
  for (i = 0; i < npts; i += 2)
    { d1 = pair[i].c1 - l1;
      d2 = pair[i].c2 - l2;
      if (d1 > d2)
        diff += d1;
      else
        diff += d2;
      l1 = pair[i+1].c1;
      l2 = pair[i+1].c2;
    }
  d1 = e1 - l1;
  d2 = e2 - l2;
  if (d1 > d2)
    diff += d1;
  else
    diff += d2;
  return (diff);
}


/*******************************************************************************************
 *
 *  Routines to test for implied relationship between pairs of reads
 *
 ********************************************************************************************/

typedef struct
  { Relation  type;       //  alignment type
    int       r1, r2;     //  between reads r1 & r2
    int       b1, e1;     //  if type != NO_ALIGNMENT then
    int       b2, e2;     //    [b1,e1] of the first read aligns to [b2,e2] of the second
    Edit     *edit;       //  if requested, implied edit script when reads overlap
    _Edit     _edit;
    int       emax;
  } _Alignment;

typedef struct
  { int64  gbeg;  //  For a block read:
    int64  gend;  //     the block is [gbeg,gend] in genome space, [hbeg,?] in haplotype space
    int64  hbeg;
  } Interval;

static Point *gfrag;

static int PSORT(const void *l, const void *r)
{ int x = *((int *) l);
  int y = *((int *) r);

  if (gfrag[x].c1 < gfrag[y].c1)
    return (-1);
  if (gfrag[x].c1 > gfrag[y].c1)
    return (1);
  return (gfrag[x].c2 - gfrag[y].c2);
}

Alignment *Align(ReadTruth *r, int64 ri, int64 rj, int min_match, int do_edit, Alignment *align)
{ Source      src1,   src2;
  _Slice      slice1, slice2;
  int         N1,     N2;
  int64       ab1,    ab2;

  //  Allocate an alignment if input align is NULL, set read pair

  if (align == NULL)
    { align = (Alignment *) Malloc(sizeof(_Alignment),"Allocating alignment record");
      if (align == NULL)
        return (NULL);
      ((_Alignment *) align)->emax = 0;
      ((_Alignment *) align)->_edit.vals = NULL;
    }
  align->r1   = ri;
  align->r2   = rj;
  align->edit = NULL;

#ifdef DEBUG_CHAIN
  printf("\nComparing read %lld to %lld\n",ri,rj);
#endif

  //  Get source, if same haplotype, check overlap and exit if not, otherwise get
  //    slices too.  In both cases, want max # of matching segments N1 & N2

  Get_Source(r,ri,&src1);
  Get_Source(r,rj,&src2);

  ab1 = src1.beg;
  ab2 = src2.beg;

  if (src1.hap == src2.hap)
    { if (src1.contig != src2.contig)
        { align->type = NO_ALIGNMENT;
          return (align);
        }
      if (src1.end < ab2 + min_match || src2.end < ab1 + min_match)
        { align->type = NO_ALIGNMENT;
          return (align);
        }
      N1 = N2 = 1;
    }
  else
    { Get_Slice(&src1,(Slice *) (&slice1));
      Get_Slice(&src2,(Slice *) (&slice2));
      N1 = slice1.lst - slice1.fst;
      N2 = slice2.lst - slice2.fst;
    }

  //  If the reads are from the same haplotype, its easy, just overlap of hap source intervals

  { Point pairs[2*(N1+N2+2)];   //  Stack allocate 
    int   npts;
    int64 b, e, ovl;
    
    if (src1.hap == src2.hap)
      { align->type = SAME_HAP_OVL;
        npts = 2;
        if (ab1 < ab2)
          b = ab2;
        else
          b = ab1;
        if (src1.end < src2.end)
          e = src1.end;
        else
          e = src2.end;
        align->b1 = pairs[0].c1 = b - ab1;
        align->b2 = pairs[0].c2 = b - ab2;
        align->e1 = pairs[1].c1 = e - ab1;
        align->e2 = pairs[1].c2 = e - ab2;
        goto phase2;
      }

    //  If different haplotypes then run the "quick test" using the compact source interval
    //    spanned by the l.i.s of source blocks and the # of bp.s not in the lis blocks

    if (slice1.gbeg < slice2.gbeg)
      ovl = slice2.gbeg;
    else
      ovl = slice1.gbeg;
    if (slice1.gend < slice2.gend)
      ovl = slice1.gend - ovl;
    else
      ovl = slice2.gend - ovl;
    if (ovl < 0)
      ovl = 0;
    if (slice1.govf < slice2.govf)
      ovl += slice1.govf;
    else
      ovl += slice2.govf;

#ifdef DEBUG_CHAIN
    printf(" <%lld-%lld:%lld> vs <%lld,%lld:%lld> -> %lld\n",
           slice1.gbeg,slice1.gend,slice1.govf,slice2.gbeg,slice2.gend,slice2.govf,ovl);
#endif

    if (ovl < min_match)
      { align->type = NO_ALIGNMENT;
        return (align);
      }

    //  Next, find the overlaps between blocks in source space, and see if the
    //    greedy chain in haplotype space is optimal

    { Interval p1[N1+1];
      Interval p2[N2+1];
      Interval *q1, *q2;
      Interval *e1, *e2;
      Block    *b1, *b2;
      int64     l1, l2, ln;
      int64     f1, f2;
      Interval *f1p, *f2p;
      int       i, chain, nmat;

      q1 = p1+N1;
      q2 = p2+N2;

      f1p = f2p = NULL;
      l1  = l2  = 0;
      f1  = f2  = 0;

      b1 = slice1.fst;
      for (i = 0; i <= N1; i++)
        { p1[i].gbeg = b1[i].beg;
          p1[i].gend = b1[i].beg + (b1[i+1].cum - b1[i].cum);
          p1[i].hbeg = b1[i].cum;
        }
      p1->gbeg = b1->beg + (ab1 - b1->cum);
      p1->hbeg = b1->cum + (ab1 - b1->cum);
      q1->gend = b1[N1].beg + (slice1.end - b1[N1].cum);

      b2 = slice2.fst;
      for (i = 0; i <= N2; i++)
        { p2[i].gbeg = b2[i].beg;
          p2[i].gend = b2[i].beg + (b2[i+1].cum - b2[i].cum);
          p2[i].hbeg = b2[i].cum;
        }
      p2->gbeg = b2->beg + (ab2 - b2->cum);
      p2->hbeg = b2->cum + (ab2 - b2->cum);
      q2->gend = b2[N2].beg + (slice2.end - b2[N2].cum);

      chain = 1;
      f1    = -1;
      ovl   = 0;
      npts  = 2;
      nmat  = 0;
      for (e1 = p1; e1 <= q1; e1++)
        for (e2 = p2; e2 <= q2; e2++)
          { if (e1->gbeg < e2->gbeg)
              b = e2->gbeg;
            else
              b = e1->gbeg;
            if (e1->gend < e2->gend)
              e = e1->gend;
            else
              e = e2->gend;

            if (b < e)
              { nmat += 1;
                if (f1 < 0)
                  { f1 = e1->hbeg + (b - e1->gbeg);
                    f2 = e2->hbeg + (b - e2->gbeg);
                    f1p = e1;
                    f2p = e2;
                  }
                else if (chain)
                  { if (e1 > p1 && e1->gend == e1[-1].gend &&
                        l1 == e1->hbeg && ln >= e1->gend - e1->gbeg)
                      {
#ifdef DEBUG_CHAIN
                        printf("Skip 1: %ld\n",e1-p1);
#endif
                        continue;
                      }
                    if (e2 > p2 && e2->gend == e2[-1].gend &&
                        l2 == e2->hbeg && ln >= e2->gend - e2->gbeg)
                      {
#ifdef DEBUG_CHAIN
                        printf("Skip 2: %ld\n",e2-p2);
#endif
                        continue;
                      }
                    if (e1->hbeg + (b - e1->gbeg) < l1 ||
                        e2->hbeg + (b - e2->gbeg) < l2)
                      chain = 0;
                  }
                l1 = e1->hbeg + (e - e1->gbeg);
                l2 = e2->hbeg + (e - e2->gbeg);
                ln = e-b;
                npts += 2;
                ovl  += e-b;
#ifdef DEBUG_CHAIN
                printf(" %2d:  %ld vs %ld :: %lld-%lld\n",npts,e1-p1,e2-p2,b,e);
#endif
	      }
          }
  
      //  If insufficent aligned segments or an obvious chain then we are done
  
      if (ovl < min_match)
        { align->type = NO_ALIGNMENT;
          return (align);
        }

      //  if the greedy chain worked, then scan it again (in linear time this go)
      //  and record the chain and alignment type

      else if (chain)
        { if (do_edit)
            { npts = 0;
              for (e1 = f1p; e1 <= q1; e1++)
                { for (e2 = f2p; e2 <= q2; e2++)
                    { if (e1->gbeg < e2->gbeg)
                        b = e2->gbeg;
                      else
                        b = e1->gbeg;
                      if (e1->gend < e2->gend)
                        e = e1->gend;
                      else
                        e = e2->gend;
                      if (b < e)
                        { if (npts > 0)
                            { if (e1 > p1 && e1->gend == e1[-1].gend &&
                                  l1 == e1->hbeg && ln >= e1->gend - e1->gbeg)
                                { f2p = e2;
                                  break;
                                }
                              if (e2 > p2 && e2->gend == e2[-1].gend &&
                                  l2 == e2->hbeg && ln >= e2->gend - e2->gbeg)
                                continue;
                            }
                          l1 = e1->hbeg + (e - e1->gbeg);
                          l2 = e2->hbeg + (e - e2->gbeg);
                          ln = e-b;
                          pairs[npts].c1 = (l1-ln) - ab1;
                          pairs[npts].c2 = (l2-ln) - ab2;
                          npts += 1;
                          pairs[npts].c1 = l1-ab1;
                          pairs[npts].c2 = l2-ab2;
                          npts += 1;
                          f2p = e2;
                          if (e == e1->gend)
                            break;
                        }
                    }
                }
            }
#ifdef DEBUG_CHAIN
          printf("CHAIN!  %lld %lld : %lld %lld\n",f1,f2,l1,l2);
#endif

          printf("CHAIN!  %lld %lld : %lld %lld\n",f1,f2,l1,l2);
          printf(" %lld %lld\n",src1.end,src2.end);
          if ((f1 == src1.beg || f2 == src2.beg) &&
              (l1 == src1.end || l2 == src2.end))
            align->type = DIFF_HAP_OVL;
          else
            align->type = DIFF_HAP_LA;

          align->b1 = f1 - ab1;
          align->e1 = l1 - ab1;
          align->b2 = f2 - ab2;
          align->e2 = l2 - ab2;
          goto phase2;
        }

      npts = (nmat<<1);

      //  At the last, we have to find the optimal chain of aligned segments in haplotype
      //    space using the O(nlogn) algorithm.

      { Point    frag[npts];
        int      sort[npts];
        int      score[npts];
        NODE     space[nmat], *free;
        NODE    *n, *o, *list;
        int      p, val;
        int      vmax, pmax, pmin, npts;

        //  Build a list of all the aligned segments/fragments end points
        //    in *read* coordinates

        npts  = 0;
        for (e1 = p1; e1 <= q1; e1++)
          for (e2 = p2; e2 <= q2; e2++)
            { if (e1->gbeg < e2->gbeg)
                b = e2->gbeg;
              else
                b = e1->gbeg;
              if (e1->gend < e2->gend)
                e = e1->gend;
              else
                e = e2->gend;
              if (b < e)
                { frag[npts].c1 = (e1->hbeg+(b-e1->gbeg)) - ab1;
                  frag[npts].c2 = (e2->hbeg+(b-e2->gbeg)) - ab2;
                  sort[npts] = npts;
                  npts += 1;
                  frag[npts].c1 = (e1->hbeg+(e-e1->gbeg)) - ab1;
		  frag[npts].c2 = (e2->hbeg+(e-e2->gbeg)) - ab2;
                  sort[npts] = npts;
                  npts += 1;
#ifdef DEBUG_CHAIN
                  printf(" %2d: %d,%d:%lld\n",npts/2+1,frag[npts-2].c1,frag[npts-2].c2,e-b);
                  printf("     %d,%d\n",frag[npts-1].c1,frag[npts-1].c2);
#endif
                }
            }

        //  Sort their start and end-points

        gfrag = frag;
        qsort(sort,npts,sizeof(int),PSORT);

#ifdef DEBUG_CHAIN
        for (i = 0; i < npts; i++)
          printf(" %6d %6d\n",frag[sort[i]].c1,frag[sort[i]].c2);
#endif

        //  Run the O(nlogn) algorithm using a splay tree to model the current frontier

        for (n = space+1; n < space+nmat; n++)
          n->L = n+1;
        space[nmat-1].L = NULL;

        vmax = pmax = 0;

        list        = space;
        list->V     = 0;
        list->L     = NULL;
        list->R     = NULL;
        list->score = 0;
        list->link  = -1;
        free = space+1;

        for (i = 0; i < npts; i++)
          { p = sort[i];
            if ((p&0x1) == 0)
              { n = FIND(list,frag[p].c2);  // c2 in [n.pos,NEXT(n).pos)
                score[p]   = n->link; 
                score[p+1] = n->score + (frag[p+1].c2 - frag[p].c2);
#ifdef DEBUG_CHAIN
                printf("%2d: In  %d,%d  -> %d (%d)\n",
                       i,frag[p].c1,frag[p].c2,score[p+1],score[p]);
#endif
              }
            else
              { o = FIND(list,frag[p].c2);
                val = score[p];
                if (val > vmax)
                  { vmax = val;
                    pmax = p;
                  }
#ifdef DEBUG_CHAIN
                printf("%2d: Out %d,%d  -> %d (%d)\n",i,frag[p].c1,frag[p].c2,val,score[p-1]);
#endif
                if (val > o->score)
                  { n = free;
                    free = n->L; 
                    n->V     = frag[p].c2;
                    n->L     = NULL;
                    n->R     = NULL;
                    n->score = val;
                    n->link  = p;
                    if (o->V == n->V)
                      list = DELETE(list,o);
                    list = INSERT(list,n);
                    while ((o = NEXT(list,n,NULL)) != NULL)
                      { if (val < o->score)
                          break;
                        list = DELETE(list,o);
                      }
                  }
#ifdef DEBUG_CHAIN
                printf("  List: ");
                PRINT_LIST(list);
                printf("\n");
#endif
              }
          }

#ifdef DEBUG_CHAIN
        for (p = pmax; p >= 0; p = score[p-1])
          printf("  B[%d] = %d\n",p/2,score[p]);
#endif

        //  Still may fail if chain not heavy enough

        if (score[pmax] < min_match)
          { align->type = NO_ALIGNMENT;
            return (align);
          }

        //  Record the best chaing and alignment type

        npts = 2;
        for (p = pmax; score[p-1] >= 0; p = score[p-1])
          npts += 2;
        pmin = p-1;

#ifdef DEBUG_CHAIN
        printf("  npts = %d\n",npts);
#endif

        if (do_edit)
          { int n = npts;

            for (p = pmax; p >= 0; p = score[p-1])
              { n -= 1;
                pairs[n].c1 = frag[p].c1;
                pairs[n].c2 = frag[p].c2;
                n -= 1;
                pairs[n].c1 = frag[p-1].c1;
                pairs[n].c2 = frag[p-1].c2;
              }
          }

        if ((frag[pmin].c1 == 0 || frag[pmin].c2 == 0) &&
            (frag[pmax].c1 == src1.end-ab1 || frag[pmax].c2 == src2.end-ab2))
          align->type = DIFF_HAP_OVL;
        else
          align->type = DIFF_HAP_LA;
        align->b1 = frag[pmin].c1;
        align->e1 = frag[pmax].c1;
        align->b2 = frag[pmin].c2;
        align->e2 = frag[pmax].c2;
        goto phase2;
      }
    }

    //  There is an overlap/local alignment

phase2:
    { _Alignment *elign;
      _Edit      *edit, edit1,  edit2;
      Map        _map1, *map1 = &_map1;
      Map        _map2, *map2 = &_map2;
      char       *ops;
      int        *vals;
      int        b1, b2;
      int        e1, e2;
      int        nsvs, nerr, nsnp;

      //  Project the perfect read coords to error read coords

      Get_Edit(r,ri,(Edit *) &edit1);
      Get_Edit(r,rj,(Edit *) &edit2);

      Init_Map(&edit1,map1);                 //  Map intervals from true to raw read space
      if (align->b1 == 0)
        b1 = 0;
      else
        b1 = Boundary_Map(map1,align->b1,1);
      if (align->e1 == src1.end - ab1)
        e1 = src1.end - ab1;
      else
        e1 = Boundary_Map(map1,align->e1,0);

      Init_Map(&edit2,map2);
      if (align->b2 == 0)
        b2 = 0;
      else
        b2 = Boundary_Map(map2,align->b2,1);
      if (align->e2 == src2.end - ab2)
        e2 = src2.end - ab2;
      else
        e2 = Boundary_Map(map2,align->e2,0);

      //  If want an edit script then:

      nsvs = nerr = nsnp = 0;
      if (do_edit)
        { Point  efull[2*(edit1.len+edit2.len+2) + npts];

          { Point  ep1[2*(edit1.len+1)];
            Point  ep2[2*(edit2.len+1)];
            Point  ehalf[2*(edit1.len+1) + npts];
            int    n1, n2, nf;

            //  Convert edit scripts to alignment blocks, and then compose with those of chain
            //    to get error read to error read matching blocks

            nsvs = Block_Diff(npts,pairs,b1,b2,e1,e2);

#ifdef DEBUG_COMPOSE
            Print_BA(npts,pairs,stdout);
#endif
            n1 = Edits_2_Blocks(&edit1,align->b1,align->e1,ep1,1);
            n2 = Edits_2_Blocks(&edit2,align->b2,align->e2,ep2,0);
            nf = Compose_Blocks(n1,ep1,npts,pairs,ehalf);
            npts = Compose_Blocks(nf,ehalf,n2,ep2,efull);

            nerr = Block_Diff(npts,efull,b1,b2,e1,e2) - nsvs;
          }

          //  Traverse the SNPs (in haplo-space, arrrrgh) in the error read to error read blocks
          //    and intertwine to form final edit script.  Two passes, first to get max size,
          //    and second stacked to fill in.

	  { SNPit _snp1, *snp1 = &_snp1;
            SNPit _snp2, *snp2 = &_snp2;
            int    p, x1, x2;
            int64  hb1, he1;
            int64  hb2, he2;
            int    do_snps;
            int    i1, i2;
            int    a1, a2;
            int    l1, l2;
            int    s1, s2, htt;
            int    del1, del2;
            int    SN, SL;

            do_snps = (src1.hap != src2.hap);
            Init_SNP(snp1,slice1.fst);
            Init_SNP(snp2,slice2.fst);

            if (do_snps)
              { Init_Map(&edit1,map1);
                Init_Map(&edit2,map2);

                for (p = 0; p < npts; p += 2)
                  { x1 = efull[p].c1;
                    i1 = Inverse_Map(map1,x1);
                    x2 = efull[p].c2;
                    i2 = Inverse_Map(map2,x2);

                    hb1 = i1 + ab1;
                    he1 = (efull[p+1].c1 - x1) + hb1;
                    hb2 = i2 + ab2;
                    he2 = (efull[p+1].c2 - x2) + hb2;

                    a1 = Next_SNP(snp1,hb1,he1);
                    a2 = Next_SNP(snp2,hb2,he2);
                    while (1)
                      { if (a1)
                          if (a2)
                            { s1 = snp1->pos - hb1;
                              s2 = snp2->pos - hb2;
                              if (snp1->pos - hb1 < snp2->pos - hb2)
                                a1 = Next_SNP(snp1,hb1,he1);
                              else if (snp1->pos - hb1 > snp2->pos - hb2)
                                a2 = Next_SNP(snp2,hb2,he2);
                              else
                                { a1 = Next_SNP(snp1,hb1,he1);
                                  a2 = Next_SNP(snp2,hb2,he2);
                                }
                            }
                          else
                            a1 = Next_SNP(snp1,hb1,he1);
		          else
                          if (a2)
                            a2 = Next_SNP(snp2,hb2,he2);
                          else
                            break;
                        nsnp += 1;
                      } 
                  }
              }
            SN = nsnp + npts + 2;

            elign = (_Alignment *) align;
            edit  = &(elign->_edit);
            align->edit = (_Edit *) edit;
            if (SN > elign->emax)
              { elign->emax = SN; 
                edit->vals = Realloc(edit->vals,(2*SN+1)*sizeof(int)+SN,
                                      "Reallocating alignment edit");
                edit->ops  = (char *) (edit->vals + (2*SN+1));
              }
            ops  = edit->ops;
            vals = edit->vals;
        
            if (do_snps)
              { Init_SNP(snp1,slice1.fst);
                Init_SNP(snp2,slice2.fst);
                Init_Map(&edit1,map1);
                Init_Map(&edit2,map2);
              }

            SN = 0;
            SL = 0;
            l1 = b1;
            l2 = b2;
            if (efull[0].c1 != b1 || efull[0].c2 != b2)
              { vals[SL++] = 0;
#ifdef DEBUG_SCRIPT
                printf("      0");
#endif
              }
            for (p = 0; p < npts; p += 2)
              { x1 = efull[p].c1;
                x2 = efull[p].c2;

                del1 = (x1 - l1);
                del2 = (x2 - l2);
                if (del1 > del2)
                  { if (del2 > 0)
                      { ops[SN++] = 's';
                        vals[SL++] = del2;
                        vals[SL++] = 0;
#ifdef DEBUG_SCRIPT
                        printf(" s %d\n      0",del2);
#endif
                      }
                    ops[SN++] = 'd';
                    vals[SL++] = del1-del2;
#ifdef DEBUG_SCRIPT
                    printf(" d %d\n",del1-del2);
#endif
                  }
                else if (del2 > del1)
                  { if (del1 > 0)
                      { ops[SN++] = 's';
                        vals[SL++] = del1;
                        vals[SL++] = 0;
#ifdef DEBUG_SCRIPT
                        printf(" s %d\n      0",del1);
#endif
                      }
                    ops[SN++] = 'i';
                    vals[SL++] = del2-del1;
#ifdef DEBUG_SCRIPT
		    printf(" i %d\n",del2-del1);
#endif
                  }
                else if (del1 > 0)
                  { ops[SN++] = 's';
                    vals[SL++] = del1;
#ifdef DEBUG_SCRIPT
                    printf(" s %d\n",del1);
#endif
                  }
                l1 = efull[p+1].c1;
                l2 = efull[p+1].c2;

                if (do_snps)
                  { i1 = Inverse_Map(map1,x1);
                    i2 = Inverse_Map(map2,x2);
 
                    hb1 = i1 + ab1;
                    he1 = (l1 - x1) + hb1;
                    hb2 = i2 + ab2;
                    he2 = (l2 - x2) + hb2;
#ifdef DEBUG_COMPOSE
                    printf("   [%5lld,%5lld] x [%5lld,%5lld]\n",hb1-ab1,he1-ab1,hb2-ab2,he2-ab2);
#endif
                    htt = 0;

                    a1 = Next_SNP(snp1,hb1,he1);
                    a2 = Next_SNP(snp2,hb2,he2);
                    while (1)
                      { if (a1)
                          if (a2)
                            { s1 = snp1->pos - hb1;
                              s2 = snp2->pos - hb2;
                              if (s1 < s2)
                                { a1 = Next_SNP(snp1,hb1,he1);
#ifdef DEBUG_SNPS
                                  printf("         %5d  (1)\n",s1);
#endif
                                }
                              else if (s1 > s2)
                                { a2 = Next_SNP(snp2,hb2,he2);
                                  s1 = s2;
#ifdef DEBUG_SNPS
                                  printf("         %5d  (2)\n",s2);
#endif
                                }
                              else
                                { a1 = Next_SNP(snp1,hb1,he1);
                                  a2 = Next_SNP(snp2,hb2,he2);
#ifdef DEBUG_SNPS
                                  printf("         %5d  (both)\n",s1);
#endif
                                }
                            }
                          else
                            { s1 = snp1->pos - hb1;
                              a1 = Next_SNP(snp1,hb1,he1);
#ifdef DEBUG_SNPS
                              printf("         %5d  (1*)\n",s1);
#endif
                            }
                        else
                          if (a2)
                            { s1 = snp2->pos - hb2;
                              a2 = Next_SNP(snp2,hb2,he2);
#ifdef DEBUG_SNPS
                              printf("         %5d  (2*)\n",s1);
#endif
                            }
                          else
                            break;
                        vals[SL++] = s1-htt;
                        ops[SN++] = 'x';
                        vals[SL++] = 1;
#ifdef DEBUG_SCRIPT
                        printf("  %5d x 1\n",s1-htt); 
#endif
                        htt = s1+1;
                      } 
                    vals[SL++] = (he1-hb1)-htt;
#ifdef DEBUG_SCRIPT
                    printf("  %5lld",(he1-hb1)-htt);
#endif
                  }
                else
                  { vals[SL++] = l1-x1;
#ifdef DEBUG_SCRIPT
                    printf("  %5d",l1-x1);
#endif
                  }
              }

            del1 = (e1 - l1);
            del2 = (e2 - l2);
            if (del1 > del2)
              { if (del2 > 0)
                  { ops[SN++] = 's';
                    vals[SL++] = del2;
                    vals[SL++] = 0;
#ifdef DEBUG_SCRIPT
                    printf(" s %d\n      0",del2);
#endif
                  }
                ops[SN++] = 'd';
                vals[SL++] = del1-del2;
                vals[SL++] = 0;
#ifdef DEBUG_SCRIPT
                printf(" d %d\n      0\n",del1-del2);
#endif
              }
            else if (del2 > del1)
              { if (del1 > 0)
                  { ops[SN++] = 's';
                    vals[SL++] = del1;
                    vals[SL++] = 0;
#ifdef DEBUG_SCRIPT
                    printf(" s %d\n      0",del1);
#endif
                  }
                ops[SN++] = 'i';
                vals[SL++] = del2-del1;
                vals[SL++] = 0;
#ifdef DEBUG_SCRIPT
	        printf(" i %d\n       0\n",del2-del1);
#endif
              }
            else if (del1 > 0)
              { ops[SN++] = 's';
                vals[SL++] = del1;
                vals[SL++] = 0;
#ifdef DEBUG_SCRIPT
                printf("S %d\n      0\n",del1);
#endif
              }
#ifdef DEBUG_SCRIPT
            else
              printf("\n");
#endif

            edit->len = SN;
          }
        }

      align->b1 = b1;
      align->e1 = e1;
      align->b2 = b2;
      align->e2 = e2;
    }
  }

  return (align);
}

void Free_Alignment_Edit(Alignment *_align)
{ _Alignment *align = (_Alignment *) _align;

  if (align->emax > 0)
    { free(align->_edit.vals);
      align->_edit.vals = NULL;
      align->emax = 0;
    }
  align->edit = NULL;
}
  

static char *Atype[4] = { "None", "True", "Hap Overlap", "Hap LA" };

void Print_Alignment(Alignment *align, FILE *file)
{ fprintf(file,"\nAlignment %s",Atype[align->type]);
  if (align->type == NO_ALIGNMENT)
    fprintf(file," %d x %d\n",align->r1,align->r2);
  else
    { fprintf(file," %d[%d..%d] x %d[%d..%d]\n",align->r1,align->b1,align->e1,
                                                align->r2,align->b2,align->e2);
      if (align->edit != NULL)
        { fprintf(file," Edit: ");
          Print_Edit(align->edit,file);
          fprintf(file,"\n");
        }
    }
}


#ifdef TESTING

/*******************************************************************************************
 *
 *  Interpret the ploidy tree and build the given number of haplotypes, returning
 *    an array of pointers to them along with their number.
 *
 ********************************************************************************************/

int main(int argc, char *argv[])
{ Genome     *gene;
  int         nhaps;
  FILE       *f;

  (void) complement;

  Prog_Name = "libtest";

  if (argc != 3)
    { fprintf(stderr,"Usage <genome.fast> <truth>\n");
      exit (1);
    }

  gene = Load_Genome(argv[1]);

  for (nhaps = 1; 1; nhaps++)
    { f = fopen(Catenate(argv[2],Numbered_Suffix(".hap",nhaps,""),"",""),"r");
      if (f == NULL)
         break;
      fclose(f);
    }
  nhaps -= 1;

  { HapTruth  *haps[nhaps];
    ReadTruth *reads;
    int        h;

    for (h = 0; h < nhaps; h++)
      { f = fopen(Catenate(argv[2],".hap",Numbered_Suffix("",h+1,""),""),"r");
        haps[h] = Load_HapTruth(f,gene);
        fclose(f);
      }

    // for (h = 0; h < nhaps; h++)
      // Print_HapTruth(haps[h],stdout);

    f = fopen(Catenate(argv[2],".err","",""),"r");
    reads = Load_ReadTruth(f,nhaps,haps);
    fclose(f);

    { Source     src;
      Edit      *edit;
      Slice     *slice;
      Alignment *align;
      int64      i, j;
      char       query[1000];

      edit  = Get_Edit(reads,0,NULL);
      slice = Get_Slice(Get_Source(reads,0,&src),NULL);

      for (i = 0; i < reads->nreads; i++)
        { Get_Source(reads,i,&src);
          if (src.hidx > 1)
            break;
        }
      printf("Reads for hap 2 start at %lld\n",i);

      while (1)
        { printf("? "); fflush(stdout);
          fgets(query,1000,stdin);
          if (*query == 'q')
            break;
          sscanf(query," %lld %lld\n",&i,&j);

          Get_Source(reads,i,&src);
          Get_Slice(&src,slice);
          printf("\nRead %lld:\n",i);
          printf("  hap %d ctg %d %c %lld %lld\n",
                 src.hidx+1,src.contig,src.orient?'R':'F',src.beg,src.end);
          printf("  %d %d: ",Slice_Length(slice),Snps_In_Slice(slice));
          Print_Slice(slice,1,stdout);

          Get_Source(reads,j,&src);
          Get_Slice(&src,slice);
          printf("\nRead %lld:\n",j);
          printf("  hap %d ctg %d %c %lld %lld\n",
                 src.hidx+1,src.contig,src.orient?'R':'F',src.beg,src.end);
          printf("  %d %d: ",Slice_Length(slice),Snps_In_Slice(slice));
          Print_Slice(slice,1,stdout);

          align = Align(reads,i,j,1000,1,NULL);
          Print_Alignment(align,stdout);
          free(align);
        }

/*
      for (i = 0; i < reads->nreads; i++)
        { Get_Source(reads,i,&src);
          Get_Edit(reads,i,edit);
          Get_Slice(&src,slice);

          printf("Read %lld:\n",i);

          printf("  hap %d ctg %d %c %lld %lld\n",
                 src.hidx+1,src.contig,src.orient?'R':'F',src.beg,src.end);

          printf("  %d %d: ",Slice_Length(slice),Snps_In_Slice(slice));
          Print_Slice(slice,0,stdout);
          printf("\n");

          printf("  %d %d: ",Edit_Length(edit),Errors_In_Edit(edit));
          Print_Edit(edit,stdout);
          printf("\n");
        }
*/
    }

    Free_ReadTruth(reads);
    for (h = 0; h < nhaps; h++)
      Free_HapTruth(haps[h]);
  }

  Free_Genome(gene);

  exit (0);
}

#endif
