/*********************************************************************************************\
 *
 *  Count errors for every 7-mer context where the middle base is the site of the edit.
 *  Also count length errors for 1-, 2-, and 3-micro satellite extensions.
 *
 *  Author:  Gene Myers
 *  Date  :  October, 2020
 *
 *********************************************************************************************/
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <dirent.h>
#include <pthread.h>
#include <math.h>

#undef  DEBUG_SPLIT      //  Original error debug flags
#undef  DEBUG_HEX
#undef  DEBUG_DEL
#undef  DEBUG_SUB
#undef  DEBUG_INS
#undef  DEBUG_MICRO_INS
#undef  DEBUG_MICRO_DEL

#undef  DEBUG_SCAN    //  Profile analysis debug flags
#undef  DEBUG_STATS
#undef  DEBUG_THREADS

#include "libfastk.h"

static char *Usage = " [-v] [-o<out>[.model]] [-T<int(4)>] -g<int>:<int> -e<int> <source>[.ktab+.prof]";

static int   GOOD_LOW, GOOD_HGH;
static int   ERROR, REPEAT;
static int   VERBOSE;
static char *OUT;
static int   NTHREADS;

static int KMER;

#define  COUNT(p)  (*((uint16 *) (p+kbyte)))
#define  MICRO(p)  (*((uint16 *) (p+tbyte)))
#define  HEPTA(p)  (*((uint16 *) (p+sbyte)))

typedef struct
  { int64  all;     //  # of occurences
    int64  del;     //  # deletion errors
    int64  sub[4];  //  # of substitution errors
    int64  ins[4];  //  # of insertion errors
    int64  mic[4];  //  # of occurences where edit point is the end of a micro with insert a
  } Edits;

typedef struct      //  For 1-, 2-, and 3-micro sats count indel errors up to unit length
  { int64 allD;
    int64 allI;
    int64 del[3];
    int64 ins[3];
  } Micro;


#if defined(DEBUG_MICRO_INS) || defined(DEBUG_MICRO_DEL)

static char *MicType[4] = { "", "Homo", "Di", "Tri" };

#endif

/****************************************************************************************
 *
 *  Print & compare utilities
 *
 *****************************************************************************************/

static char dna[4] = { 'a', 'c', 'g', 't' };

static void print_seq(uint8 *seq, int len)
{ int i, b, k;
  int khalf = KMER/2;

  b = -1;
  k = 0;
  for (i = 0; i < khalf-1; i++)
    { if (k == 0)
        { k = 6; b += 1; }
      else
        k -= 2;
      printf("%c",dna[(seq[b] >> k) & 0x3]);
    }
  for (; i < khalf+2; i++)
    { if (k == 0)
        { k = 6; b += 1; }
      else
        k -= 2;
      if (i == khalf)
        printf(".");
      printf("%c",dna[(seq[b] >> k) & 0x3]-32);
    }
  for (; i < len; i++)
    { if (k == 0)
        { k = 6; b += 1; }
      else
        k -= 2;
      printf("%c",dna[(seq[b] >> k) & 0x3]);
    }
}


/****************************************************************************************
 *
 *  Comparison Routines over bit-compressed k-mers
 *
 *****************************************************************************************/

  //  Return the length of the longest common prefix of a and b

static inline int mypref(uint8 *a, uint8 *b, int n)
{ int   i;
  uint8 x, y;
  
  for (i = 0; i <= n; i += 4)
    { if (*a != *b)
        { x = *a;
          y = *b;
          if ((x & 0xf0) != (y & 0xf0))
            if ((x & 0xc0) != (y & 0xc0))
              return (i);
            else
              return (i + 1);
          else
            if ((x & 0xfc) != (y & 0xfc))
              return (i + 2);
            else
              return (i + 3);
        }
      a += 1;
      b += 1;
    }
  return (n+1);
}

#define CHAR_INDEX(k)      ((k) >> 2)
#define CHAR_SHIFT(k)      (6 - 2*(k&0x3))
#define CHAR_AT(ptr,ki,ks) ((ptr[ki] >> ks) & 0x3)

  //  Assumes la <= lb: compare suffix of la chars of a and b

static inline int mycmp(uint8 *a, int la, uint8 *b, int lb)
{ int ka, sa;
  int kb, sb;
  int j, x, y;

  ka = CHAR_INDEX(KMER-la);
  sa = CHAR_SHIFT(KMER-la);

  kb = CHAR_INDEX(KMER-lb);
  sb = CHAR_SHIFT(KMER-lb);

  for (j = 0; j < la; j++)
    { x = CHAR_AT(a,ka,sa);
      y = CHAR_AT(b,kb,sb);
      if (x != y)
        return (x-y);
      if (sa == 0)
        { sa = 6;
          ka += 1;
        }
      else
        sa -= 2;
      if (sb == 0)
        { sb = 6;
          kb += 1;
        }
      else
        sb -= 2;
    }
  return (0);
}


/****************************************************************************************
 *
 *  Threaded routine to accumulate error statistics
 *
 *****************************************************************************************/

typedef struct
  { Kmer_Stream *T;
    int64        end;
    Edits        heptab[0x4000];
    Micro       *mic3tab[0x40]; 
    Micro       *mic2tab[0x10]; 
    Micro       *mic1tab[0x04]; 
    int          hex[0x1001];
    int          heb[0x1001];
  } Error_Parm;

static void *measure_thread(void *args)
{ Error_Parm  *parm    = (Error_Parm *) args;
  Kmer_Stream *T       = parm->T;
  int64        end     = parm->end;
  Edits       *heptab  =  parm->heptab;
  Micro      **mic3tab = parm->mic3tab;
  Micro      **mic2tab = parm->mic2tab;
  Micro      **mic1tab = parm->mic1tab;
  int         *hex     = parm->hex;
  int         *heb     = parm->heb;

  int    kbyte = T->kbyte;
  int    tbyte = T->tbyte;
  int    sbyte = tbyte+2;
  int    mbyte = sbyte+2;

  int    khalf = KMER/2;
  int    khulf = KMER - khalf;
  int    kmici = khulf+1;
  int    kmicd = khulf-2;

  int    km2i = CHAR_INDEX(khalf-1); // index & mask for chars khalf-1, khalf, & khalf+1
  int    km2s = CHAR_SHIFT(khalf-1);
  int    kh2i = CHAR_INDEX(khalf);
  int    kh2s = CHAR_SHIFT(khalf);
  int    kp2i = CHAR_INDEX(khalf+1);
  int    kp2s = CHAR_SHIFT(khalf+1);

  uint8 *cache, *ctop;
  uint8 *finger[65];

  //  Setup micro detection tables, hex[6-mer] = micro unit-length thereof (0 if not a micro)
  //                                heb[6-mer] = 1st base of unit

  { uint32 x;
    int   *hx = hex+1;
    int   *hb = heb+1;

    bzero(hex,sizeof(int)*0x1001);

    for (x = 0; x < 4096; x += 16)
      { hx[x|0x0] = hx[x|0x5] = hx[x|0xa] = hx[x|0xf] = 1;

        hb[x|0x0] = 0;
        hb[x|0x5] = 1;
        hb[x|0xa] = 2;
        hb[x|0xf] = 3;
      }

    for (x = 0; x < 4096; x += 256)
      { hx[x|0x11] = hx[x|0x22] = hx[x|0x33] = hx[x|0x44] = 2;
        hx[x|0x66] = hx[x|0x77] = hx[x|0x88] = hx[x|0x99] = 2;
        hx[x|0xbb] = hx[x|0xcc] = hx[x|0xdd] = hx[x|0xee] = 2;

        hb[x|0x11] = hb[x|0x22] = hb[x|0x33] = 0;
        hb[x|0x44] = hb[x|0x66] = hb[x|0x77] = 1;
        hb[x|0x88] = hb[x|0x99] = hb[x|0xbb] = 2;
        hb[x|0xcc] = hb[x|0xdd] = hb[x|0xee] = 3;
      }

    for (x = 1; x < 63; x++)
      if (x != 0x15 && x != 0x2a)
        { hx[x << 6 | x] = 3;
          hb[x << 6 | x] = (x >> 4);
        }
  }

  //   Zero statistics tables

  { int     i;
    Micro  *mic1, *mic2, *mic3;

    mic1 = ((Micro *) Malloc(sizeof(Micro)*0x04*(khalf-6),"Homo-Table")) - 2;
    mic2 = ((Micro *) Malloc(sizeof(Micro)*0x10*(khalf-6),"Di-Table"))   - 4;
    mic3 = ((Micro *) Malloc(sizeof(Micro)*0x40*(khalf-6),"Tri-Table"))  - 6;

    for (i = 0; i < 0x04; i++)
      mic1tab[i] = mic1 + (khalf-6)*i;

    for (i = 0; i < 0x10; i++)
      mic2tab[i] = mic2 + (khalf-6)*i;

    for (i = 0; i < 0x40; i++)
      mic3tab[i] = mic3 + (khalf-6)*i;

    bzero(mic3+6,sizeof(Micro)*0x40*(khalf-6));
    bzero(mic2+4,sizeof(Micro)*0x10*(khalf-6));
    bzero(mic1+2,sizeof(Micro)*0x4*(khalf-6));
    bzero(heptab,sizeof(Edits)*0x400);
  }

  cache = Malloc(4097*mbyte,"Allocating entry buffer");  //  Processing block
  ctop  = cache + 4096*mbyte;

  while (T->cidx < end)
    {

    //  Load cache with the next set of entries whose first khalf-1 chars are equal.
    //  Keep track of the 64 sub-segements for each continuation of 3 bases in the
    //   array of finger pointers.

      { uint8 *nptr, *cptr;
        int    f, x;

        cptr = cache;
        Current_Entry(T,cptr);
        f = 0;
        x = (CHAR_AT(cptr,km2i,km2s) << 4)
          | (CHAR_AT(cptr,kh2i,kh2s) << 2)
          |  CHAR_AT(cptr,kp2i,kp2s);
        while (f <= x)
          finger[f++] = cptr;
        nptr = cptr+mbyte;
        for (Next_Kmer_Entry(T); T->csuf != NULL; Next_Kmer_Entry(T))
          { x = mypref(cptr,Current_Entry(T,nptr),khalf+1); 
            if (x < khalf-1)
              break;
            if (x <= khalf+1)
              { x = (CHAR_AT(nptr,km2i,km2s) << 4)
                  | (CHAR_AT(nptr,kh2i,kh2s) << 2)
                  |  CHAR_AT(nptr,kp2i,kp2s);
                while (f <= x)
                  finger[f++] = nptr;
              }
            if (nptr >= ctop)
              { int64 cidx = ctop-cache;
                int64 cmax = ((cidx*14)/(10*mbyte) + 2048)*mbyte; 
                cache = Realloc(cache,cmax+mbyte,"Reallocting entry buffer");
                ctop  = cache + cmax;
                nptr  = cache + cidx;
                for (x = 1; x < f; x++)
                  finger[x] = cache + (finger[x] - finger[0]);
                finger[0] = cache;
              }
            cptr = nptr;
            nptr = cptr+mbyte;
          }
        while (f <= 64)
          finger[f++] = nptr;
      }

#ifdef DEBUG_SPLIT
      { int i;

        printf("\nLoad Cache:\n");
        for (i = 0; i < 64; i++)
          { if (finger[i+1] > finger[i])
              { printf(" %02x: %ld\n",i,(finger[i]-cache)/mbyte);
                printf("  ");
                print_seq(finger[i],KMER);
                printf("\n");
                printf("  ");
                print_seq(finger[i+1]-mbyte,KMER);
                printf("\n");
              }
          }
        printf(" 40: %ld\n",(finger[64]-cache)/mbyte);
      }
#endif

      //  Find polysat tips and mark good heptamers @khalf

      { uint8 *f, *fe;
        int    i, s, t, cn, ma;
        uint32 x, y, z;

        f = cache;
        i = CHAR_INDEX(khalf-6);
        s = CHAR_SHIFT(khalf-6);
        z = CHAR_AT(f,i,s);
        for (t = 1; t < 5; t++)
          { if (s == 0)
              { s = 6; i += 1; }
            else
              s -= 2;
            z = (z << 2) | CHAR_AT(f,i,s);
          }

        for (fe = finger[64]; f < fe; f += mbyte)
          { x = (z << 2) | CHAR_AT(f,km2i,km2s);
            y = (((x++ << 2) & 0xfff) | CHAR_AT(f,kh2i,kh2s)) + 1;
            if (hex[x] > 0 && hex[y] != hex[x])
              { *((uint16 *) (f+tbyte)) = x;
                ma = heb[x];
              }
            else
              { *((uint16 *) (f+tbyte)) = 0;
                ma = -1;
              }
#ifdef DEBUG_HEX
            print_seq(f,KMER);
            printf(" %5d",COUNT(f));
            printf(" [ %03x %1d %03x %1d %c %c ]",x-1,hex[x],y-1,hex[y],dna[heb[x]],
                                                  (hex[x] > 0 && hex[y] != hex[x]) ? 'Y' : 'N');
#endif
            x = (y-1) & 0x3ff;
            i = kh2i;
            s = kh2s;
            for (t = 0; t < 2; t++)
              { if (s == 0)
                  { s = 6; i += 1; }
                else
                  s -= 2;
                x = (x << 2) | CHAR_AT(f,i,s);
              }
            cn = COUNT(f);
            if (GOOD_LOW <= cn && cn <= GOOD_HGH)
              { *((uint16 *) (f+sbyte)) = x+1;
                heptab[x].all += cn;
                if (ma >= 0)
                  heptab[x].mic[ma] += cn;
              }
            else
              *((uint16 *) (f+sbyte)) = 0;
#ifdef DEBUG_HEX
            printf(" < %04x %c >\n",x,(GOOD_LOW <= cn && cn <= GOOD_HGH) ? 'Y' : 'N');
#endif
          }
      }

      //  Deletion: scan Ax.y? versus Ay.? for all x!=y and Ax. is not a micro end

      { uint8 *f, *fe;
        uint8 *g, *ge;
        int    i, c;
        int    en;

#ifdef DEBUG_DEL
        printf("\nDeletion Scan:\n");
#endif
        for (i = 0; i < 16; i++)
          { if (i % 5 == 0)
              continue; 
            c  = (i<<2);
            f  = finger[c];
            fe = finger[c + 4];
#ifdef DEBUG_DEL
            printf(" %02x-%02x (%ld-%ld)  :: ",c,c+4,(f-cache)/mbyte,(fe-cache)/mbyte);
#endif
            c  = ((i&0x3) << 4);
            g  = finger[c];
            ge = finger[c + 16];
#ifdef DEBUG_DEL
            printf("%02x-%02x (%ld-%ld)\n",c,c+16,(g-cache)/mbyte,(ge-cache)/mbyte);
#endif
            while (1)
              { if (f < fe)
                  if (g < ge)
                    c = mycmp(f,khulf-1,g,khulf);
                  else
                    break;
                else
                  break;
#ifdef DEBUG_DEL
                printf(" %2ld ",(f-cache)/mbyte);
                print_seq(f,KMER);
                printf(" %2ld ",(g-cache)/mbyte);
                print_seq(g,KMER);
                printf(" %d\n",c);
#endif
                if (c < 0)
                  f += mbyte;
                else if (c > 0)
                  g += mbyte;
                else
                  { en = COUNT(g);
                    g += mbyte;
                    if (MICRO(f) == 0 && HEPTA(f) > 0)
                      { while (g < ge && mycmp(f,khulf-1,g,khulf) == 0)
                          { en += COUNT(g);
                            g  += mbyte;
                          }
                        if (en <= ERROR)
                          { heptab[HEPTA(f)-1].del += en;
#ifdef DEBUG_DEL
                            printf("Del %c : %d %d : %04x\n",dna[i>>2],COUNT(f),en,HEPTA(f)-1);
#endif
                          }
                      }
                    f += mbyte;
                  }
              }
          }
      }

      //  Substitution: Scan/merge Ax.? for all x

      { uint8 *ep[4];
        uint8 *p[4];
        int    stk[4], stp, in;
        int    i, s, c, n;

#ifdef DEBUG_SUB
        printf("\nSubstitution Scan:\n");
#endif
        n = 0;
        for (i = 0; i < 4; i++)
          { c = (i<<4);
            p[i]  = finger[c];
            ep[i] = finger[c + 16];
            if (p[i] < ep[i])
              { n += 1;
#ifdef DEBUG_SUB
                printf(" %d (%ld-%ld) \n",i,(p[i]-cache)/mbyte,(ep[i]-cache)/mbyte);
                fflush(stdout);
#endif
              }
          }
        if (n >= 2)
         while (1)
          { for (i = 0; i < 4; i++)
              { if (p[i] < ep[i])
                  break;
              }
            if (i > 3)
              break;
            s = i;
            stk[0] = i;
            stp = 1;
            for (i++; i < 4; i++)
              { if (p[i] >= ep[i])
                  continue;
                c = mycmp(p[s],khulf,p[i],khulf);
                if (c > 0)
                  { s = i;
                    stk[0] = i;
                    stp = 1;
                  }
                else if (c == 0)
                  stk[stp++] = i;
              }
#ifdef DEBUG_SUB
            for (s = 0; s < 4; s++)
              if (p[s] < ep[s])
                { printf(" %ld ",(p[s]-cache)/mbyte);
                  print_seq(p[s],KMER);
                }
            for (s = 0; s < stp; s++)
              printf("  %d",stk[s]);
            printf("\n");
#endif
            if (stp == 1)
              { p[stk[0]] += mbyte;
                continue;
              }
            for (s = 0; s < stp; s++)
              if (HEPTA(p[stk[s]]) > 0)
                { for (i = 0; i < stp; i++)
                    if (i != s)
                      { in = COUNT(p[stk[i]]);
                        if (in <= ERROR)
                          { heptab[HEPTA(p[stk[s]])-1].sub[stk[i]] += in;
#ifdef DEBUG_SUB
                            printf("Sub %c for %c : %d %d : %04x\n",
                                   dna[stk[s]],dna[stk[i]],COUNT(p[stk[s]]),in,HEPTA(p[stk[s]])-1);
#endif
                          }
                      }
                }
            for (s = 0; s < stp; s++)
              p[stk[s]] += mbyte;
	  }
      }

      //  Insertion: scan Ax.yz? vs Ax.z? for y!=z and Ax.y? is not a micro extension

      { uint8 *f, *fe;
        uint8 *g, *ge;
        int    i, a, c, h;

#ifdef DEBUG_INS
        printf("\nInsertion Scan:\n");
#endif
        for (i = 0; i < 64; i++)
          { if ((i&0xf) % 5 == 0)
              continue;
            f  = finger[i];
            fe = finger[i+1];
            if (f >= fe)
              continue;
            c  = (i&0x30) | ((i&0x3)<<2);
            g  = finger[c];
            ge = finger[c+4];
            if (g >= ge)
              continue;
#ifdef DEBUG_INS
            printf(" %02x (%ld-%ld)  :: ",i,(f-cache)/mbyte,(fe-cache)/mbyte);
            printf("%02x-%02x (%ld-%ld)\n",c,c+4,(g-cache)/mbyte,(ge-cache)/mbyte);
#endif

            a = (i>>2) & 0x3;
            h = MICRO(g);
            if (h > 0 && heb[h] == a)
              {
#ifdef DEBUG_INS
                printf("  Micro Skip\n");
#endif
                continue;
              }

            while (1)
              { if (f < fe)
                  if (g < ge)
                    c = mycmp(f,khulf-1,g,khulf);
                  else
                    break;
                else
                  break;
#ifdef DEBUG_INS
                printf(" %ld ",(f-cache)/mbyte);
                print_seq(f,KMER);
                printf(" %ld ",(g-cache)/mbyte);
                print_seq(g,KMER);
                printf(" %d\n",c);
#endif
                if (c < 0)
                  f += mbyte;
                else if (c > 0)
                  g += mbyte;
                else
                  { if (COUNT(f) <= ERROR)
                      { do
                          { if (HEPTA(g) > 0)
                              { heptab[HEPTA(g)-1].ins[a] += COUNT(f);
#ifdef DEBUG_INS
                                printf("Ins %c : %d %d : %04x\n",
                                       dna[a],COUNT(g),COUNT(f),HEPTA(g)-1);
#endif
                              }
                            g += mbyte;
                          }
                        while (g < ge && mycmp(f,khulf-1,g,khulf) == 0);
                      }
                    else
                      { do
                          g += mbyte;
                        while (g < ge && mycmp(f,khulf-1,g,khulf) == 0);
                      }
                    f += mbyte;
                  }
              }
          }
      }

      //  Micro-satellite insertions
      //    If khalf-1 prefix A ends in a micro x of unit length m (in 1,2,3) then:
      //       Determine its length k and scan Ac against Ay for all prefixes y of x
      //           where c is not the 1st symbol of x

      { Micro *mic;
        uint8 *f, *fe, *fs;
        uint8 *g[3], *e[3];
        uint32 x, z;
        int i, k;
        int m;
        int s, c;

        f = finger[0];
        c = CHAR_INDEX(khalf-7);
        s = CHAR_SHIFT(khalf-7);
        x = CHAR_AT(f,c,s);
        for (i = 1; i < 6; i++)
          { if (s == 0)
              { s = 6; c += 1; }
            else
              s -= 2;
            x = (x << 2) | CHAR_AT(f,c,s);
          }
        if ((m = hex[x+1]) > 0)
          { c = CHAR_INDEX(khalf-7);
            s = CHAR_SHIFT(khalf-7);
            z = x;
            for (k = khalf-8; k >= 0; k--)
              { if (s == 6)
                  { s = 0; c -= 1; }
                else
                  s += 2;
                z = (z >> 2) | (CHAR_AT(f,c,s) << 10);
                if (hex[z+1] != m)
                  { if (m == 1 && hex[z+1] == 3)
                      { if ((z & 0xf) % 5 == 0)
                          k -= 1;
                      }
                    break;
                  }
              } 

            if (k >= 0 || m == 3 || hex[(z>>2)+1] != m)
              { k = (khalf + 2*m) - (8+k);

                //  m=1: k in [2,khalf-5], m=2: k in [4,khalf-3], m=3: k in [6,khalf-1]

                x = ((x << 2*(3-m)) & 0x3f);
#ifdef DEBUG_MICRO_INS
                printf("\n  Mod %d micro of length %d, ext = %03x\n",m,k,x);
#endif

                for (i = 0; i < m; i++)
                  { s = (1 << (2*(2-i)));
                    c = x & (64 - s);
                    g[i] = finger[c];
                    e[i] = finger[c+s];
#ifdef DEBUG_MICRO_INS
                    printf("%d: %02x-%02x (%ld-%ld)\n",
                           i,c,c+s,(g[i]-cache)/mbyte,(e[i]-cache)/mbyte);
#endif
                  }

                x >>= 2*(3-m);
                if (m == 1)
                  mic = mic1tab[x] + k;
                else if (m == 2)
                  mic = mic2tab[x] + k;
                else
                  mic = mic3tab[x] + k;

	        fe = finger[64];
                fs = g[0];
                while (1)
                  { if (f == fs)
                      f = e[0];
                    if (f >= fe)
                      break;
                    for (i = 0; i < m; i++)
                      { for ( ; g[i] < e[i]; g[i] += mbyte)
                          { c = mycmp(g[i],kmici-(i+1),f,kmici);
                            if (c == 0)
                              {
#ifdef DEBUG_MICRO_INS
                                printf(" %ld ",(f-cache)/mbyte);
                                print_seq(f,KMER);
                                printf(" %ld ",(g[i]-cache)/mbyte);
                                print_seq(g[i],KMER);
                                printf(" %d\n",i);
#endif
                                if (HEPTA(f) > 0 && COUNT(g[i]) <= ERROR)
                                  { mic->ins[i] += COUNT(g[i]);
#ifdef DEBUG_MICRO_INS
                                    printf("%s %d:%d Ins %d : %d %d\n",
                                           MicType[m],x,k,j,COUNT(f),COUNT(g[i]));
#endif
                                  }
                                break;
                              }
                            else if (c > 0)
                              break;
                          }
                      }
                    if (HEPTA(f) > 0)
                      mic->allI += COUNT(f);
                    f += mbyte;
                  }
              }
	  }
      }

      //  Micro-satellite deletions
      //    If khalf+2 prefix Aux (|ux|=3) ends in a micro x of unit length m (in 1,2,3) then:
      //       Determine its length k and scan Auy against Aux for all prefixes y of x making
      //          sure AuxB is a micro tip. 

      { Micro *mic;
        uint8 *f, *fe;
        uint8 *g[3], *e[3];
        uint32 x, z;
        int i, k, t;
        int m;
        int s, c, a;
        int d, r;
        int dnt;

        d  = CHAR_INDEX(khalf+2);
        r  = CHAR_SHIFT(khalf+2);
        for (t = 0; t < 64; t++)
          { if (finger[t] >= finger[t+1])
              continue;
            f = finger[t];
            c = CHAR_INDEX(khalf-4);
            s = CHAR_SHIFT(khalf-4);
            x = CHAR_AT(f,c,s);
            for (i = 1; i < 6; i++)
              { if (s == 0)
                  { s = 6; c += 1; }
                else
                  s -= 2;
                x = (x << 2) | CHAR_AT(f,c,s);
              }
            if ((m = hex[x+1]) > 0)
              { c = CHAR_INDEX(khalf-4);
                s = CHAR_SHIFT(khalf-4);
                z = x;
                for (k = khalf-4; k > 2; k--)
                  { if (s == 6)
                      { s = 0; c -= 1; }
                    else
                      s += 2;
                    z = (z >> 2) | (CHAR_AT(f,c,s) << 10);
                    if (hex[z+1] != m)
                      { if (m == 1 && hex[z+1] == 3)
                          { if ((z & 0xf) % 5 == 0)
                              k -= 1;
                          }
                        break;
                      }
                  } 
                if (k == 2)
                  break;
                k += 4 - 2*m;
                k = khalf - k;

#ifdef DEBUG_MICRO_DEL
                printf("\n  Mod %d micro of length %d, ext = %03x\n",m,k,t);
#endif

                for (i = 0; i < m; i++)
                  { s = (1 << (2*(i+1)));
                    c = t & (64 - s);
                    g[i] = finger[c];
                    e[i] = finger[c+s];
#ifdef DEBUG_MICRO_DEL
                    printf("%d: %02x-%02x (%ld-%ld)\n",
                           i,c,c+s,(g[i]-cache)/mbyte,(e[i]-cache)/mbyte);
#endif
                  }

                a  = heb[x+1];
                x &= (1 << 2*m) - 1;
                if (m == 1)
                  mic = mic1tab[x] + k;
                else if (m == 2)
                  mic = mic2tab[x] + k;
                else
                  mic = mic3tab[x] + k;

                fe = finger[t+1];
                while (1)
                  { while (CHAR_AT(f,d,r) == a)
                      { f += mbyte;
                        if (f >= fe)
                          break;
                      }
                    if (f >= fe)
                      break;
                    for (i = 0; i < m; i++)
                      { dnt = 0;
                        for ( ; g[i] < e[i]; g[i] += mbyte)
                          { c = mycmp(f,kmicd,g[i],kmicd+i+1);
                            if (c < 0)
                              break;
                            if (c == 0)
                              dnt += COUNT(g[i]);
                          }
                        if (dnt > 0)
                          {
#ifdef DEBUG_MICRO_DEL
                            printf(" %ld ",(f-cache)/mbyte);
                            print_seq(f,KMER);
                            printf(" %ld ",(g[i]-cache)/mbyte);
                            print_seq(g,KMER);
                            printf(" %d %d\n",c,i);
#endif
                
                            if (HEPTA(f) > 0 && dnt <= ERROR)
                              { mic->del[i] += dnt;
#ifdef DEBUG_MICRO_DEL
                                printf("%s %d:%d Del %d : %d %d\n",
                                      MicType[m],x,k,j,COUNT(f),COUNT(g));
#endif
                              }
                          }
                      }
                    if (HEPTA(f) > 0)
                      mic->allD += COUNT(f);
                    f += mbyte;
                  }
              }
          }
      }

    }

  free(cache);

  return (NULL);
}


/****************************************************************************************
 *
 *  Analyze Profiles
 *
 *****************************************************************************************/

#define RBUCK      50   //  Histogram buckets for real length are RBUCK bases
#define EBUCK     500   //  Histogram buckets for errors are EBUCK bases where EBUCK % RBUCK = 0
#define RLIMIT 100000
#define ELIMIT    500

typedef struct
  { Profile_Index *P;
    int64          beg;
    int64          end;
    int64         *len;
    int64        **elb;
    int64          rsum;
    int64          rsqr;
    int64          ncnt;
    int64          terr;
    int64          tbps;
  } Read_Parm;

static void *scan_thread(void *args)
{ Read_Parm *parm = (Read_Parm *) args;
  Profile_Index *P = parm->P;
  int64          beg = parm->beg;
  int64          end = parm->end;
  int64         *len = parm->len;
  int64        **elb = parm->elb;

  uint16 *profile;
  int     pmax, plen;
  int     ebeg, eend;
  int     rbeg, rend;
  int     rep, err;
  int64   terr, tbps;
  int64   rsum, rsqr, ncnt;
  int64   i;
  int     p;

  pmax    = 20000;
  profile = Malloc(pmax*sizeof(uint16),"Profile array");

  terr = tbps = 0;
  ncnt = 0;
  rsum = rsqr = 0;
  for (i = beg; i < end; i++)
    { plen = Fetch_Profile(P,i,pmax,profile);
      if (plen > pmax)
        { pmax    = 1.2*plen + 1000;
          profile = Realloc(profile,pmax*sizeof(uint16),"Profile array");
          Fetch_Profile(P,i,pmax,profile);
        }

      rbeg = rend = -1;
      ebeg = eend = -1;
      rep = 0;
      err = 0;
      for (p = 0; p < plen; p++)
        { if (profile[p] >= REPEAT)
            { if (ebeg >= 0)
                { err += (eend - ebeg) / KMER + 1;
                  ebeg = eend = -1;
                }
              if (rbeg < 0)
                rbeg = p;
              rend = p+1;
            }
          else if (rbeg >= 0)
            { if (p >= rend+3*KMER)
                { rep += rend - rbeg;
                  p = rend;
                  rbeg = rend = -1;
                }
            }
          else if (profile[p] <= ERROR)
            { if (ebeg < 0)
                ebeg = p;
              eend = p;
            }
          else if (ebeg >= 0)
            { err += (eend - ebeg) / KMER + 1;
              ebeg = eend = -1;
            }
#ifdef DEBUG_SCAN
          printf("%5d: %5d => [%5d,%5d] [%5d,%5d]  %5d %5d\n",
                 p,profile[p],ebeg,eend,rbeg,rend,err,rep);
#endif
        }
          
      if (rbeg >= 0)
        rep += plen-rbeg;
      if (ebeg >= 0)
        err += (eend - ebeg) / KMER + 1;
#ifdef DEBUG_SCAN
      printf("%5d:   *   => [%5d,%5d] [%5d,%5d]  %5d %5d\n",
             p,ebeg,eend,rbeg,rend,err,rep);
#endif

      if (rep > .5*plen)
        continue;

      terr += err;
      tbps += plen-rep;

      ncnt += 1;
      err   = (10000*err)/(plen-rep);
      plen += KMER-1;
      rsum += plen;
      rsqr += plen*plen;

#ifdef DEBUG_STATS
      printf(" %4.2f%%  %d\n",err/100.,plen);
#endif

      if (plen > RLIMIT)
        { len[RLIMIT/RBUCK] += 1;
          if (err > ELIMIT)
            elb[RLIMIT/EBUCK][ELIMIT] += 1;
          else
            elb[RLIMIT/EBUCK][err] += 1;
        }
      else
        { len[plen/RBUCK] += 1;
          if (err > ELIMIT)
            elb[plen/EBUCK][ELIMIT] += 1;
          else
            elb[plen/EBUCK][err] += 1;
        }
    }

  free(profile);

  parm->ncnt = ncnt;
  parm->rsum = rsum;
  parm->rsqr = rsqr;
  parm->terr = terr;
  parm->tbps = tbps;
  return (NULL);
}

/****************************************************************************************
 *
 *  Main
 *
 *****************************************************************************************/


static uint8 comp[128] =
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 3, 0, 2, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 3, 0, 2, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

static void compress_comp(char *s, int len, uint8 *t)
{ int    i;
  char   c, d, e;
  char  *s0, *s1, *s2, *s3;

  s0 = s;
  s1 = s0-1;
  s2 = s1-1;
  s3 = s2-1;

  c = s1[0];
  d = s2[0];
  e = s3[0];
  s1[0] = s2[0] = s3[0] = 't';

  for (i = len-1; i >= 0; i -= 4)
    *t++ = ((comp[(int) s0[i]] << 6) | (comp[(int) s1[i]] << 4)
         |  (comp[(int) s2[i]] << 2) | comp[(int) s3[i]] );

  s1[0] = c;
  s2[0] = d;
  s3[0] = e;
}

static void examine_table(Kmer_Stream *T, int *full, int *sym)
{
  //  Look at middle 100M counts and see if has counts to 1

  { int hbyte = T->hbyte;

    for (GoTo_Kmer_Index(T,0); T->cidx < T->nels; Next_Kmer_Entry(T))
      if (*((int16 *) (T->csuf+hbyte)) == 1)
        break;

    *full = (T->cidx < T->nels);
  }

  //  Walk to a non-palindromic k-mer and see if its complement is in T

  { int64  sidx;
    char  *seq;
    uint8 *cmp;
    int    kmer;

    kmer = T->kmer;

    sidx = 1;
    GoTo_Kmer_Index(T,sidx);
    seq = Current_Kmer(T,NULL);
    cmp = Current_Entry(T,NULL);
    while (1)
      { compress_comp(seq,kmer,cmp);
        if (GoTo_Kmer_Entry(T,cmp))
          { if (T->cidx != sidx)
              { *sym = 1;
                break;
              }
          }
        else
          { *sym = 0;
            break;
          }
        sidx += 1;
        seq = Current_Kmer(T,seq);
      }
    free(cmp);
    free(seq);
  }
}

int main(int argc, char *argv[])
{ FILE          *Efile;

  Kmer_Stream   *T;
  static Edits   HepTab[0x4000];
  static Micro  *Mic3Tab[0x40];
  static Micro  *Mic2Tab[0x10];
  static Micro  *Mic1Tab[0x4];

  Profile_Index *P;
  int64         *Len;
  int64        **Elb;
  int64          RMean;
  int64          RsDev;
  int64          Nread;
  int            AvErr;

  (void) print_seq;

  { int    i, j, k;
    int    flags[128];
    char  *eptr, *fptr;

    ARG_INIT("HImodel");

    GOOD_LOW = -1;
    ERROR    = -1;
    OUT      = NULL;
    NTHREADS = 4;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("v")
            break;
          case 'e':
            ARG_POSITIVE(ERROR,"Error threshold");
            break;
          case 'g':
            GOOD_LOW = strtol(argv[i]+2,&eptr,10);
            if (eptr > argv[i]+2)
              { if (GOOD_LOW < 1 || GOOD_LOW > 0x7fff)
                  { fprintf(stderr,"%s: Good minimum count %d is out of range\n",
                                   Prog_Name,GOOD_LOW);
                    exit (1);
                  }
                if (*eptr == ':')
                  { GOOD_HGH = strtol(eptr+1,&fptr,10);
                    if (fptr > eptr+1 && *fptr == '\0')
                      { if (GOOD_HGH < 1 || GOOD_HGH > 0x7fff)
                          { fprintf(stderr,"%s: Good maximum count %d is out of range\n",
                                           Prog_Name,GOOD_HGH);
                            exit (1);
                          }
                        if (GOOD_LOW > GOOD_HGH)
                          { fprintf(stderr,"%s: Good count range is invalid\n",Prog_Name);
                            exit (1);
                          }
                        break;
                      }
                  }
              }
            fprintf(stderr,"%s: Syntax of -g option invalid -h<int>:<int>\n",Prog_Name);
            exit (1);
          case 'o':
            OUT = argv[i]+2;
            break;
          case 'T':
            ARG_POSITIVE(NTHREADS,"Number of threads")
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    VERBOSE = flags['v'];

    if (argc != 2 || GOOD_LOW <= 0 || ERROR <= 0)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"      -g: count range (inclusive) of k-mers to be considered good.\n");
        fprintf(stderr,"      -e: counts below which k-mer is considered erroneous.\n");
        fprintf(stderr,"      -o: name of output model file.\n");
        exit (1);
      }

    REPEAT = 3*GOOD_HGH;
    if (OUT == NULL)
      { if (strcmp(argv[1]+(strlen(argv[1])-5),".ktab") == 0)
          OUT = Root(argv[1],"ktab");
        else
          OUT = Root(argv[1],"prof");
      }
    else
      OUT = Root(OUT,"model");
  }

  { Error_Parm *parm;
    pthread_t   threads[NTHREADS];
    int         full, sym;
    char        command[1000];
    char       *troot;
    char        template[15] = "._MODEL.XXXX";
    int         khalf;
    int         t;

    P = Open_Profiles(argv[1]);
    if (P == NULL)
      { fprintf(stderr,"%s: Cannot open profile index %s\n",Prog_Name,argv[1]);
        exit (1);
      }

    T = Open_Kmer_Stream(argv[1]);
    if (T == NULL)
      { fprintf(stderr,"%s: Cannot open table %s\n",Prog_Name,argv[1]);
        exit (1);
      }

    KMER  = T->kmer;
    khalf = KMER/2;

    examine_table(T,&full,&sym);

    if (!full)
      { fprintf(stderr,"%s: Input table has no k-mers with a count of 1 ???\n",Prog_Name);
        exit (1);
      }

    if (!sym)
      { if (VERBOSE)
          { fprintf(stderr,"\n  Making table symmetric\n");
            fflush(stderr);
          }

        Free_Kmer_Stream(T);

        troot = mktemp(template);

        sprintf(command,"Symmex -T8 %s %s",argv[1],troot);

        if (system(command) != 0)
          { fprintf(stderr,"%s: Something went wrong with command:\n    %s\n",Prog_Name,command);
            exit (1);
          }

        T = Open_Kmer_Stream(troot);
      }

    parm = Malloc(sizeof(Error_Parm)*NTHREADS,"ALlocating parms");
    if (parm == NULL)
      exit (1);

    { Kmer_Stream *U;
      int64        p;
      uint8       *ent, wm;
      int          wb, i;
#ifdef DEBUG_SPLIT
      char      *seq;

      seq = Current_Kmer(T,NULL);
#endif
      ent = Current_Entry(T,NULL);
      wb = (khalf-2)>>2;
      switch ((khalf-1)&0x3)
      { case 0: wm = 0xff; break;
        case 1: wm = 0xc0; break;
        case 2: wm = 0xf0; break;
        case 3: wm = 0xfc; break;
      }

      parm[0].T = T;
      First_Kmer_Entry(T);
      for (t = 1; t < NTHREADS; t++)
        { parm[t].T = U = Clone_Kmer_Stream(T);
          p = (T->nels*t)/NTHREADS; 
          GoTo_Kmer_Index(U,p);
#ifdef DEBUG_SPLIT
          printf("\n%d: %0*x\n",t,2*U->ibyte,U->cpre);
          GoTo_Kmer_Index(U,p-1);
          printf(" %lld: %s\n",p-1,Current_Kmer(U,seq));
          GoTo_Kmer_Index(U,p);
          printf(" %lld: %s\n",p,Current_Kmer(U,seq));
#endif  
          Current_Entry(U,ent);                //  Break at prefix boundaries
          ent[wb] &= wm;
          for (i = wb+1; i < T->kbyte; i++)      
            ent[i] = 0;
          GoTo_Kmer_Entry(U,ent);
#ifdef DEBUG_SPLIT
          GoTo_Kmer_Index(U,U->cidx-1);
          printf(" %lld: %s\n",U->cidx,Current_Kmer(U,seq));
          GoTo_Kmer_Index(U,U->cidx+1);
          printf(" %lld: %s\n",U->cidx,Current_Kmer(U,seq));
#endif  
          parm[t-1].end = U->cidx;
        } 
      parm[NTHREADS-1].end = T->nels;
      free(ent);
#ifdef DEBUG_SPLIT
      free(seq);
#endif
    }

#ifdef DEBUG_THREADS
    for (t = 0; t < NTHREADS; t++)
      measure_thread(parm+t);
#else 
    for (t = 1; t < NTHREADS; t++)
      pthread_create(threads+t,NULL,measure_thread,parm+t);
    measure_thread(parm);
    for (t = 1; t < NTHREADS; t++)
      pthread_join(threads[t],NULL);
#endif

    { int     i, k, j;
      Micro **mic;
      Edits  *edt;

      for (i = 0; i < 0x04; i++)
        Mic1Tab[i] = parm[0].mic1tab[i];
      for (i = 0; i < 0x10; i++)
        Mic2Tab[i] = parm[0].mic2tab[i];
      for (i = 0; i < 0x40; i++)
        Mic3Tab[i] = parm[0].mic3tab[i];
      for (i = 0; i < 0x4000; i++)
        HepTab[i]  = parm[0].heptab[i];

      for (t = 1; t < NTHREADS; t++)
        { mic = parm[t].mic1tab;
          for (i = 0; i < 0x04; i++)
            for (k = 2; k < khalf-4; k++)
              { Mic1Tab[i][k].allD += mic[i][k].allD;
                Mic1Tab[i][k].allI += mic[i][k].allI;
                for (j = 0; j < 3; j++)
                  { Mic1Tab[i][k].del[j] += mic[i][k].del[j];
                    Mic1Tab[i][k].ins[j] += mic[i][k].ins[j];
                  }
              }
          mic = parm[t].mic2tab;
          for (i = 0; i < 0x10; i++)
            for (k = 4; k < khalf-2; k++)
              { Mic2Tab[i][k].allD += mic[i][k].allD;
                Mic2Tab[i][k].allI += mic[i][k].allI;
                for (j = 0; j < 3; j++)
                  { Mic2Tab[i][k].del[j] += mic[i][k].del[j];
                    Mic2Tab[i][k].ins[j] += mic[i][k].ins[j];
                  }
              }
          mic = parm[t].mic3tab;
          for (i = 0; i < 0x40; i++)
            for (k = 6; k < khalf; k++)
              { Mic3Tab[i][k].allD += mic[i][k].allD;
                Mic3Tab[i][k].allI += mic[i][k].allI;
                for (j = 0; j < 3; j++)
                  { Mic3Tab[i][k].del[j] += mic[i][k].del[j];
                    Mic3Tab[i][k].ins[j] += mic[i][k].ins[j];
                  }
              }
          edt = parm[t].heptab;
          for (i = 0; i < 0x4000; i++)
            { HepTab[i].all += edt[i].all;
              HepTab[i].del += edt[i].del;
              for (j = 0; j < 4; j++)
                { HepTab[i].sub[j] += edt[i].sub[j];
                  HepTab[i].ins[j] += edt[i].ins[j];
                  HepTab[i].mic[j] += edt[i].mic[j];
                }
            }
        }
    }

    for (t = NTHREADS-1; t >= 1; t--)
      { free(parm[t].mic3tab[0]+6);
        free(parm[t].mic2tab[0]+4);
        free(parm[t].mic1tab[0]+2);
        Free_Kmer_Stream(parm[t].T);
      }

    free(parm);

    Free_Kmer_Stream(T);

    if (!sym)
      { sprintf(command,"Fastrm -f %s",troot);
        system(command);
      }
  }

  //  Write stats tables as rates tables to "Error.tab"

  { typedef struct
      { float  all;
        float  ins;
        float  op[9];
      } E_Rates;

    typedef struct
      { float all;
        float op[6];
      } M_Rates;

    int      i, j, k, n;
    int      mic_size[4] = { 1, 4, 16, 64 };
    double   d;
    int      khalf;
    E_Rates  rate;
    Micro   *mtab;
    M_Rates  dmic;

    khalf = KMER/2;

    Efile = fopen(Catenate(OUT,".model","",""),"w");

    fwrite(&KMER,sizeof(int),1,Efile);

    for (i = 0; i < 0x4000; i++)
      { d = HepTab[i].all;
        if (d == 0)
          { rate.all = 0.;
            rate.ins = 0.;
            for (j = 0; j < 9; j++)
              rate.op[j] = 0.;
          }
        else
          { for (j = 0; j < 4; j++)
              if (d > HepTab[i].mic[j])
                rate.op[j] = HepTab[i].ins[j]/(d-HepTab[i].mic[j]);
              else
                rate.op[j] = 0.;
            for (j = 0; j < 4; j++)
              rate.op[j+4] = HepTab[i].sub[j]/d;
            for (j = 0; j < 4; j++)
              d -= HepTab[i].mic[j];
            if (d > 0)
              rate.op[8] = HepTab[i].del/d;
            else
              rate.op[8] = 0.;
            rate.all = 0.;
            rate.ins = 0.;
            for (j = 0; j < 9; j++)
	      rate.all += rate.op[j];
            for (j = 0; j < 4; j++)
	      rate.ins += rate.op[j];

            if (VERBOSE)
              { printf("\n");
                for (j = 12; j >= 0; j -= 2)
                  printf("%c",dna[(i>>j)&0x3]);
                printf("\n");
                d = HepTab[i].all;
                for (j = 0; j < 4; j++)
                  if (d > HepTab[i].mic[j])
                    printf("  I %c = %10lld / (%10.0f - %10lld) = %.3f\n",
                           dna[j],HepTab[i].ins[j],d,HepTab[i].mic[j],
                           (100.*HepTab[i].ins[j])/(d-HepTab[i].mic[j]));
                  else
                    printf("  I %c = %10lld / (%10.0f - %10lld) = %.3f\n",
                           dna[j],HepTab[i].ins[j],d,HepTab[i].mic[j],0.);
                for (j = 0; j < 4; j++)
                  printf("  S %c = %10lld / %10.0f = %.3f\n",
                         dna[j],HepTab[i].sub[j],d,(100.*HepTab[i].sub[j])/d);
                for (k = 0; k < 4; k++)
                  d -= HepTab[i].mic[k];
                if (d > 0)
                  printf("  D   = %10lld / %10.0f = %.3f\n",
                         HepTab[i].del,d,(100.*HepTab[i].del/d));
                else
                  printf("  D   = %10lld / %10.0f = %.3f\n",
                         HepTab[i].del,d,0.);
              }
          }
        fwrite(&rate,sizeof(E_Rates),1,Efile);
      }

    for (n = 1; n <= 3; n++)
      for (i = 0; i < mic_size[n]; i++)
        { if (n == 1)
            mtab = Mic1Tab[i]+2;
          else if (n == 2)
            mtab = Mic2Tab[i]+4;
          else
            mtab = Mic3Tab[i]+6;
          for (k = 0; k < khalf-6; k++)
            { dmic.all = 0.;
              d = mtab[k].allD;
              if (d == 0)
                { for (j = 0; j < n; j++)
                    dmic.op[j] = 0.;
                }
              else
                { for (j = 0; j < n; j++)
                    { dmic.op[j] = mtab[k].del[j]/d;
                      dmic.all  += dmic.op[j];
                    }
                }
              d = mtab[k].allI;
              if (d == 0)
                { for (j = 0; j < n; j++)
                    dmic.op[j+n] = 0.;
                } 
              else
                { for (j = 0; j < n; j++)
                    { dmic.op[j+n] = mtab[k].ins[j]/d;
                      dmic.all    += dmic.op[j+n];
                    }
                }
              fwrite(&dmic,sizeof(M_Rates),1,Efile);
            }

          if (VERBOSE)
            { printf("\n%d ",n);
              for (k = 2*(n-1); k >= 0; k -= 2)
                printf("%c",dna[(i>>k)&0x3]);
              printf("\n");
              for (k = 0; k < khalf-6; k++)
                { d = mtab[k].allD;
                  if (d > 0)
                    { printf("  D %2d:",k+2*n);
                      for (j = 0; j < n; j++)
                        printf(" %10lld",mtab[k].del[j]);
                      printf(" / %10.0f = ",d);
                      for (j = 0; j < n; j++)
                        printf(" %.3f",(100.*mtab[k].del[j])/d);
                      printf("\n");
                    }
                }
              for (k = 0; k < khalf-6; k++)
                { d = mtab[k].allI;
                  if (d > 0)
                    { printf("  I %2d:",k+2*n);
                      for (j = 0; j < n; j++)
                        printf(" %10lld",mtab[k].ins[j]);
                      printf(" / %10.0f = ",d);
                      for (j = 0; j < n; j++)
                        printf(" %.3f",(100.*mtab[k].ins[j])/d);
                      printf("\n");
                    }
                }
            }
        }
  }

  free(Mic1Tab[0]+2);
  free(Mic2Tab[0]+4);
  free(Mic3Tab[0]+6);

  //  Code for analyzing read profiles to extract length distribution and length-dependent
  //    error rate distributions

  { Read_Parm parm[NTHREADS];
    pthread_t threads[NTHREADS];
    int64 etot, ebps;
    int   t, i, e;

    for (t = 0; t < NTHREADS; t++)
      { parm[t].len  = Malloc(sizeof(int64)*(RLIMIT/RBUCK+1),"Allocating stat tables");
        parm[t].elb  = Malloc(sizeof(int64 *)*(RLIMIT/EBUCK+1),"Allocating stat tables");
        parm[t].elb[0] = Malloc(sizeof(int64)*(ELIMIT+1)*(RLIMIT/EBUCK+1),"Allocating stat tables");

        for (e = 1; e <= RLIMIT/EBUCK; e++)
          parm[t].elb[e] = parm[t].elb[e-1] + (ELIMIT+1);
        bzero(parm[t].len,sizeof(int64)*(RLIMIT/RBUCK+1));
        bzero(parm[t].elb[0],sizeof(int64)*(ELIMIT+1)*(RLIMIT/EBUCK+1));
      }

    Elb = parm[0].elb;
    Len = parm[0].len;

    for (t = 0; t < NTHREADS; t++)
      { parm[t].beg  = (P->nreads * t) / NTHREADS;
        parm[t].end  = (P->nreads * (t+1)) / NTHREADS;
        if (t > 0)
          parm[t].P = Clone_Profiles(P);
        else
          parm[0].P = P;
      }

    KMER = P->kmer;

#ifdef DEBUG_THREADS
    for (t = 0; t < NTHREADS; t++)
      scan_thread(parm+t);
#else
    for (t = 1; t < NTHREADS; t++)
      pthread_create(threads+t,NULL,scan_thread,parm+t);
    scan_thread(parm);
    for (t = 1; t < NTHREADS; t++)
      pthread_join(threads[t],NULL);
#endif

    RMean = parm[0].rsum;
    RsDev = parm[0].rsqr;
    Nread = parm[0].ncnt;
    etot  = parm[0].terr;
    ebps  = parm[0].tbps;
    for (t = 1; t < NTHREADS; t++)
      { for (i = 0; i <= RLIMIT/RBUCK; i++)
          Len[i] += parm[t].len[i];
        for (i = 0; i <= RLIMIT/EBUCK; i++)
          for (e = 0; e <= ELIMIT; e++)
            Elb[i][e] += parm[t].elb[i][e];
        RMean += parm[t].rsum;
        RsDev += parm[t].rsqr;
        Nread += parm[t].ncnt;
        etot  += parm[t].terr;
        ebps  += parm[t].tbps;
      }
    RMean /= Nread;
    RsDev = sqrt (RsDev/Nread - RMean*RMean);
    AvErr = (10000*etot)/ebps;

    for (t = NTHREADS-1; t > 0; t--)
      Free_Profiles(parm[t].P);

    for (t = NTHREADS-1; t > 0; t--)
      { free(parm[t].elb[0]);
        free(parm[t].elb);
        free(parm[t].len);
      }
  }

  { int    i, e;
    int    rmax, emax, itmp;
    int64  sum, etot;
    double cum;

    fwrite(&RMean,sizeof(int64),1,Efile);
    fwrite(&RsDev,sizeof(int64),1,Efile);

    if (Len[RLIMIT/RBUCK] > 0)
      fprintf(stderr,"\nWarning: Data set has reads over %dKbp\n",(RLIMIT*RBUCK)/1000);
    for (rmax = RLIMIT/RBUCK; rmax >= 0; rmax--)
      if (Len[rmax] > 0)
        break;

    fwrite(&rmax,sizeof(int),1,Efile);
    itmp = RBUCK;
    fwrite(&itmp,sizeof(int),1,Efile);
    itmp = EBUCK;
    fwrite(&itmp,sizeof(int),1,Efile);
    itmp = ELIMIT;
    fwrite(&itmp,sizeof(int),1,Efile);
    fwrite(&AvErr,sizeof(int),1,Efile);

    sum = 0;
    for (i = 0; i <= rmax; i++)
      { sum += Len[i];
        cum = (1.*sum)/Nread;
        fwrite(&cum,sizeof(double),1,Efile);
      }

    emax = ((rmax+1)*RBUCK-1)/EBUCK;
    for (i = 0; i <= emax; i++)
      { etot = 0;
        for (e = 0; e <= ELIMIT; e++)
          etot += Elb[i][e];
        if (etot == 0)
          { cum = 0.;
            for (e = 0; e <= ELIMIT; e++)
              fwrite(&cum,sizeof(double),1,Efile);
          }
        else
          { sum = 0;
            for (e = 0; e <= ELIMIT; e++)
              { sum += Elb[i][e];
                cum = (1.*sum)/etot;
                fwrite(&cum,sizeof(double),1,Efile);
              }
          }
      }

    if (VERBOSE)
      { printf("\nHistogram of Real Lengths: Mean = %lld SDev = %lld Averr = %.2f%%\n",
               RMean,RsDev,AvErr/100.);
        sum = 0;
        for (i = 0; i <= rmax; i++)
          if (Len[i] > 0)
            { sum += Len[i];
              printf(" %5d: %10lld %.5f\n",i*RBUCK,Len[i],(1.*sum)/Nread);
            }
        printf("\n");

        printf("\nHistograms of Error Rates\n");
        for (i = 0; i <= emax; i++)
          { etot = 0;
            for (e = 0; e <= ELIMIT; e++)
              etot += Elb[i][e];
            if (etot == 0)
              continue;
            printf("%5d - %5d:\n",i*EBUCK,(i+1)*EBUCK);
            sum = 0;
            for (e = 0; e <= ELIMIT; e++)
              if (Elb[i][e] > 0)
                { sum += Elb[i][e];
                  printf("  %4.2f%%: %10lld %.5f\n",e/100.,Elb[i][e],(1.*sum)/etot);
                }
          }
      }

    fclose(Efile);
  }

  Free_Profiles(P);

  free(OUT);

  Catenate(NULL,NULL,NULL,NULL);
  Numbered_Suffix(NULL,0,NULL);
  free(Prog_Name);

  exit (0);
}
