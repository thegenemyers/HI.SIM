/*******************************************************************************************
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

#include "gene_core.h"
#include "lib_sim.h"

int main(int argc, char *argv[])
{ Genome     *gene;
  int         nhaps;
  FILE       *f;

  Prog_Name = "HIreads";

  if (argc != 3)
    { fprintf(stderr,"Usage <genome>[.fast[aq]] <truth:.err+.hap#>\n");
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
  if (nhaps == 0)
    { fprintf(stderr,"%s: Cannot find haplotype ground truth file %s.hap1\n",Prog_Name,argv[2]);
      exit (1);
    }

  { HapTruth   *haps[nhaps];
    ReadTruth  *reads;
    int         h;

    for (h = 0; h < nhaps; h++)
      { f = fopen(Catenate(argv[2],".hap",Numbered_Suffix("",h+1,""),""),"r");
        haps[h] = Load_HapTruth(f,gene);
        fclose(f);
      }

    f = fopen(Catenate(argv[2],".err","",""),"r");
    if (f == NULL)
      { fprintf(stderr,"%s: Cannot find read ground truth file %s.err\n",Prog_Name,argv[2]);
        exit (1);
      }
    reads = Load_ReadTruth(f,nhaps,haps);
    fclose(f);

    Print_Fasta(reads,stdout);

    Free_ReadTruth(reads);
    for (h = 0; h < nhaps; h++)
      Free_HapTruth(haps[h]);
  }

  Free_Genome(gene);

  exit (0);
}
