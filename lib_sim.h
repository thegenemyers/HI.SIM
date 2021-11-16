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

int64   Genome_Length(Genome *gene);
int     Scaffold_Count(Genome *gene);
int64   Scaffold_Length(Genome *gene, int i);
char   *Scaffold_Sequence(Genome *gene, int i);


typedef void HapTruth;

HapTruth  *Load_HapTruth(FILE *file, Genome *gene);
void       Free_HapTruth(HapTruth *hap);
int64      Size_Of_HapTruth(HapTruth *hap);
void       Print_HapTruth(HapTruth *hap, FILE *file);
Genome    *Haplotype_Genome(HapTruth *hap);


typedef struct
  { int64  nreads;      //  # of reads in data set
    float *rate;        //  for i in [0,H-1] rate of variation from source for haplotype i+1
    void  *hidden[2];
  } ReadTruth;

ReadTruth *Load_ReadTruth(FILE *file, int nhaps, HapTruth **haps);
int64      Size_Of_ReadTruth(ReadTruth *reads);
void       Free_ReadTruth(ReadTruth *reads);
int        Print_Fasta(ReadTruth *reads, FILE *file);


typedef struct
  { HapTruth  *hap;      //  haplotype read sampled from
    int        hidx;     //  index of haplotype starting at 1
    int        contig;   //  scaffold/contig of haplotype starting at 1
    int        orient;   //  0 => forward, 1 => reverse
    int64      beg;      //  interval [beg,end]
    int64      end;
  } Source;

Source *Get_Source(ReadTruth *r, int64 i, Source *source);
void    Free_Source(Source *source);


typedef void Slice;

Slice  *Get_Slice(Source *source, Slice *slice);
void    Free_Slice(Slice *slice);
void    Print_Slice(Slice *slice, int show_snps, FILE *file);
int     Slice_Length(Slice *slice);
int     Snps_In_Slice(Slice *slice);
uint8  *Slice_Sequence(Slice *slice, uint8 *seq);
uint8  *True_Sequence(ReadTruth *r, int64 i, uint8 *seq);


typedef void Edit;

Edit   *Get_Edit(ReadTruth *r, int64 i, Edit *edit);
void    Free_Edit(Edit *edit);
void    Print_Edit(Edit *edit, FILE *file);
int     Edit_Length(Edit *edit);
int     Errors_In_Edit(Edit *edit);
int     Error_Type_In_Edit(Edit *edit, char type);
uint8  *Edit_Sequence(Edit *edit, uint8 *in, uint8 *out);
uint8  *Read_Sequence(ReadTruth *r, int64 i, uint8 *in, uint8 *out);


typedef enum { NO_ALIGNMENT = 0, SAME_HAP_OVL = 1, DIFF_HAP_OVL = 2, DIFF_HAP_LA = 3 } Relation;

typedef struct
  { Relation  type;       //  alignment type
    int       r1, r2;     //  between reads r1 & r2
    int       b1, e1;     //  if type != NO_ALIGNMENT then
    int       b2, e2;     //    [b1,e1] of the first read aligns to [b2,e2] of the second
    Edit     *edit;       //  if requested, implied edit script when reads overlap
    int       nsvs;       //  # of differences due to structural variations between hap's
    int       nsnp;       //  # of differences due to SNPs between hap's
    int       nerr;       //  # of differences due to errors in the reads
    void     *hiddent[4];
  } Alignment;

Alignment *Align(ReadTruth *r, int64 i, int64 j, int min_match, int do_edit, Alignment *align);
void       Free_Alignment_Edit(Alignment *align);
void       Print_Alignment(Alignment *align, FILE *file);

#endif  // _LIB_HISIM
