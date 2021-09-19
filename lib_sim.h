#ifndef _LIB_HISIM
#define _LIB_HISIM

#include "gene_core.h"
#include "libfastk.h"


typedef void *Error_Model;

Error_Model *Load_Error_Model(char *name);
float       *Read_Error_Profile(Error_Model *epro, uint8 *read, int len, float *profile);


typedef void Genome;

Genome *Load_Genome(char *name);
void    Free_Genome(Genome *gene);
int64   Size_Of_Genome(Genome *gene);
void    Print_Genome(Genome *gene, FILE *file);


typedef void Haplotype;

Haplotype *Load_Haplotype(FILE *file, Genome *gene);
void       Free_Haplotype(Haplotype *hap);
int64      Size_Of_Haplotype(Haplotype *hap);
void       Print_Haplotype(Haplotype *hap, FILE *file);
Genome    *Haplotype_Sequence(Haplotype *hap);


typedef struct
  { int64  nreads;
    float *rate;
    void  *hidden[2];
  } Reads;

Reads *Load_Reads(FILE *file, int nhaps, Haplotype **haps);
void   Free_Reads(Reads *reads);


typedef struct
  { Haplotype *hap;
    int        contig;
    int        orient;
    int64      beg;
    int64      end;
  } Source;

Source *Get_Read_Source(Reads *r, int64 i, Source *source);
void    Free_Read_Source(Source *source);


typedef void Slice;

Slice  *Get_Read_Slice(Source *source, Slice *slice);
void    Free_Slice(Slice *slice);
void    Print_Slice(Slice *slice, FILE *file);
int     Slice_Length(Slice *slice);
int     Snps_In_Slice(Slice *slice);
uint8  *Get_Slice_Sequence(Slice *slice, uint8 *seq);
uint8  *Get_True_Sequence(Reads *r, int64 i, uint8 *seq);


typedef void Edit;

Edit   *Get_Read_Edit(Reads *r, int64 i, Edit *edit);
void    Free_Read_Edit(Edit *edit);
void    Print_Read_Edit(Edit *edit, FILE *file);
int     Edit_Length(Edit *edit);
int     Errors_In_Edit(Edit *edit);
uint8  *Edit_Sequence(Edit *edit, uint8 *in, uint8 *out);
uint8  *Get_Read_Sequence(Reads *r, int64 i, uint8 *in, uint8 *out);

#endif  // _LIB_HISIM
