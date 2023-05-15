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
{ Genome     *gene, *hapseq;
  HapTruth   *hap;
  int         hapno;
  FILE       *f;

  Prog_Name = "HIhaplo";

  if (argc != 4)
    { fprintf(stderr,"Usage <genome>[.fast[aq]] <truth:.hap#> <int>\n");
      exit (1);
    }

  { char *eptr;

    hapno = strtol(argv[3],&eptr,10);
    if (*eptr != '\0')
      { fprintf(stderr,"%s: '%s' argument is not an integer\n",Prog_Name,argv[3]);
        exit (1);
      }
    if (hapno <= 0)
      { fprintf(stderr,"%s: Haplotype number must be positive (%d)\n",Prog_Name,hapno);
        exit (1);
      }
  }

  gene = Load_Genome(argv[1]);

  f = fopen(Catenate(argv[2],Numbered_Suffix(".hap",hapno,""),"",""),"r");
  if (f == NULL)
    { fprintf(stderr,"%s: Cannot find haplotype ground truth file %s.hap%d\n",
                      Prog_Name,argv[2],hapno);
      exit (1);
    }

  hap = Load_HapTruth(f,gene);
  fclose(f);

  hapseq = Haplotype_Genome(hap);

  Print_Genome(hapseq,stdout);

  Free_Genome(hapseq);
  Free_HapTruth(hap);
  Free_Genome(gene);

  exit (0);
}
