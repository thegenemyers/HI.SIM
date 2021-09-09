#ifndef _LIB_HISIM
#define _LIB_HISIM

#include "gene_core.h"
#include "libfastk.h"

typedef void *Error_Model;

Error_Model *Load_Error_Model(char *name);
float       *Read_Error_Profile(Error_Model *epro, uint8 *read, int len, float *profile);


typedef void *Genome;

Genome *Load_Genome(char *name);
void    Free_Genome(Genome *gene);
int64   Size_Of_Genome(Genome *gene);
void    Print_Genome(Genome *gene, FILE *file);


typedef void *Haplotype;
typedef void *Hrange;

Haplotype *Load_Haplotype(FILE *file, Genome *gene);
void       Free_Haplotype(Haplotype *hap);
int64      Size_Of_Haplotype(Haplotype *hap);
void       Print_Haplotype(Haplotype *hap, FILE *file);

Genome    *Haplotype_Sequence(Haplotype *hap);
char      *Haplotype_Read(Haplotype *hap, int64 beg, int64 end);

typedef void *Read_Truth;

Read_Truth *Load_Read_Truth(FILE *file);

#endif  // _LIB_HISIM
