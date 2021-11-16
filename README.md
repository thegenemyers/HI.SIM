# A HiFi Shotgun Simulator
  
<font size ="4">**_Author:  Gene Myers_**<br>
**_First:   Sept. 1, 2021_**<br>
**_Current: Sept. 1, 2021_**</font>

- [Commands](#command-line)
  - [Himodel](#himodel)
  - [HIsim](#hisim)
  - [HIfasta](#hifasta)
- [The Error Model](#error-model)

- [Ground Truth C-Library](#C-lib)
            
<a name="command-line"></a>

## Commands

This module contains a program, **HImodel**, that builds an error and read
length model for a highly accurate read data set such as produced by Pacbio machines in CCS mode, and another, **HIsim**, that uses the model to produce a simulated
shotgun read data set for a hypothetical multi-ploid genome with accompanying, easy-to-interpret, ground truth.
These programs were developed with the aim of allowing one
to properly test and evaluate a shotgun assembler
for HiFI data or other highly accurate read data.

<a name="himodel"></a>

```
1. HImodel [-v] [-o<out>[.model] -g<int>:<int> -e<int> <source>[.ktab+prof]
```

HImodel takes a **FastK** k&#8209;mer table, \<source>.ktab, and read profiles, \<source>.prof, of a HiFi shotgun read data set as input and from it produces an error and read length model for such data.
The k&#8209;mer table must be produced by FastK with option <code>-t1</code> so that k&#8209;mers that occur as rarely as once in the data set are in the table.  For HiFI data we further recommend that k be on the order of 40 or more.

For efficiency reasons, HImodel, needs a k&#8209;mer table that is *symmetric*,
i.e. every k&#8209;mer occurs in both orientations in the table.  Native [FastK](http:url)
tables are *canonical*, i.e. every k&#8209;mer occurs once in its minimum orientation, but these can be converted into a larger symmetric version with the program [Symmex](httt:url).  One can give as input to HImodel either a symmetric or a canonical table.  The difference is that in the later case, HImodel must take more time to produce the symmetric table internally.

The model is based on the inference that k&#8209;mers with less than a count given by the -e threshold are likely errors, and those with a count in
the closed interval given by the -g argument are likely valid haploid
or diploid k&#8209;mers.  To determine these thresholds for a given input data
set one typically examines the histogram of the k&#8209;mers with say [Histex](http:url)
or even more carefully with a tool like [GenomeScope.FK](http:url).

An encoding of the error model is placed in a file with the extension
<code>.model</code>.  By default the file takes its root name from the source, i.e. <code>\<source>.model</code>.  One can give a specific
root path with the -o option.  With the -v option set, HImodel, outputs
a representation of the set of tables for the model (see the description below).
            
<a name="hisim"></a>

```
2. HIsim <genome>[.fast[aq]] <model>[.model] [-p<tree>] [-o<name>] [-hef]
         [-c<double(50.)>] [-m<int(10000)>] [-s<int(2000)>] [x<int(4000)>]
         [-vCU] [-w<int(100)>] [-r<int>]
```

The FastK simulator, **HIsim**, takes as input a source genomic sequence, a
read error model and an optional "ploidy tree".
It first produces haplotype sequences from the genomic sequence according
to the ploidy tree if given, otherwise it simply takes the genomic
sequence as the sole haplotype sequence.  It then simulates shotgunning
reads from the haplotype sequences with equal likelihood, introducing
errors into the reads according to the supplied error model.

**Source Sequence.**  The source genome is a fasta file and should *not* be
a haplo-type phased assembly as haplotypes will be generated by the
simulator so that the ground truth of haplotypic variation is known.
What is desired is simply a realistic genomic sequence with the repetitive
and structural features of a real genome (which are not simulated).
So one should use either an unphased reconstruction or the primary haplotype
sequence of a phased reconstruction.
If the input is a sequence of scaffolds with runs of N's denoting gaps
between contigs, then all the gaps are removed so that each scaffold
is a single contiguous sequence.

**Error Model.** The model for introducting errors into reads and for determining their lengths and overall error rates is one produced by **HImodel**.  The simulated shotgun reads will thus have an error profile and read length profile similar to that of the data set the model was extracted from.  In particular, read lengths are generated with the distribution observed
for the data set from which the model was generated.  If the -m or -s options
are given, then these are further offset and scaled so that the generated
read lengths have the given mean and/or standard deviation.  The model also
includes the distribution of error rates seen for individual reads as a function
of their length so that HIsim generates reads with varying overall error rates
presumably reflective of the number of "ccs wraps" for the read.

**Ploidy Tree.** The ploidy tree guides how haplotypes are "evolved" from the source genome.
Conceptually it is a tree with mutation rates on each edge.  The source
genome is placed at the root and then mutated by the rate on each edge
emanating from it to produce a sequence at each child of the root.  These sequences are then recursively mutated by
the rate on each edge emanating from them, and so on, until there is a sequence
at each leaf.  The sequences at the leaves collectively represent the haplotype instances of the source genome generated by the simulator.
For example the tree <code>.2(.1,.2),.3</code> specifies 3 haplotypes
where the first is a .1% mutation of a genome X that is a .2% mutation
of the source genome.  The second is a .2% mutation of X, and the third
is a .3% mutation of the source.  In general, a ploidy tree
<code>\<tree></code> is specified with the simple grammar:

```
    <tree> <- <tree> ( ',' <tree> )*
            | <rate> [ '(' <tree> ')' ]
```

where rate is any fixed precision number representing a *percentage*
mutation rate.

**Haplotype Generation.** Given a rate of mutation &rho;, a genome sequence of length G will have
&rho;G bases of structural variations and another &rho;G bases of individual
SNP mutations.  The structural variations are either block deletions or
block insertions, each equally likely.  The size of the blocks is generated
with a Paretto distribution with parameters such that the minimum (and median) block
size is 1, and the mean is 50bp.  Insertion
blocks are either *random* or *tandem* with equal probability.  A random
insertion takes a block chosen at random from the genome and inserts it
at the selected insertion point, whereas a tandem insertion takes the
block ending at the insertion point and inserts it at the insertion point, effectively duplicating
the block in tandem.  The various rates are all defined constants in the
underlying code and can be changed if desired, but we wanted to keep the
design simple and so did not make these constants explicit user parameters.
Our sense is that this model of haplotype variation is sufficient for test
and evaluation purposes.

**Read Generation.**  Reads are sampled from each haplotype with equal probability until the coverage set by the -c option is reached (or 50X
by default).  The read length distribution is taken from the supplied
error model, \<model>, scaled and offset if necessary so that the mean and
standard deviation match the -m and -s options of the command line.
Reads less than the -x option threshold are ignored.  Reads are sampled from either strand with equal probability.  The error read for the
reads is selected from a distribution supplied by the error model as a function
of the read length to the nearest Kbp.  Given an error rate, errors are then
introduced into the read using the supplied model, with the probability of
individual rates scaled so that the read acheives the selected error rate.s

The reads are output in fasta format either to stdout if the -f option
is *not* set or to <code>\<name>.fasta</code> if it is set, where
<code>\<name></code> is either the path of the -o option if present,
or the root of the genome fasta file.

**Ground Truth / Haplotypes.**
If the -h option is given, then a [ONE-code](http:url) encoding of
the relationship of each haplotype to the source genome is output to
files with the names <code>\<name>.hap#</code> where # is the number of
the haplotype (starting with 1) and <code>\<name></code> is either the
path of the -o option if present, or the root of the genome fasta file.  Each haplotype is described as a
sequence of blocks of the source where each block can have a number
of SNPs within it.  Each block is simply an interval of the source,
and each SNP is an offset within the block and the variant at that
offset.  The variant is specified as an offset of 1, 2, or 3 from the
source base at that location (e.g. G+1 is T, G+2 is A, and G+3 is C).

The ONE-code file consists of a sequence of lines the first of which is
always <code>1 3 hap</code> indicating that this file is a ".hap" ONE-code
file.  The file then describes a sequence of chromosomes/contigs where
the following grammar specifies the syntax of a chromosome:

```
<chrom> <- || 'c' <blocks:int> <length:int> ||  <block>+
<block> <- || 'B' <beg-coord:int> <end-coord:int> ||  [<snp>]
<snps>  <- || 'P' <n:int> <off_1:int> ... <off_n:int> ||
           || 'S' <n:int> <var_1:[1-3]> ... <var_n:[1-3]> ||
```

The tokens between || are individual lines and it is guaranteed that each
line begins with the given character and tokens are separated by a single blank.
The 'c'-line (c for chromosome/contig) gives the # of blocks that follow
and the total length of the chromosome in base pairs.  The 'B'-line (B for
block) gives the start and end coordinates of the block in the source genome.  Here coordinates are "ticks" between bases starting at 0.  The
'P'- and 'S'-lines (P for position, S for SNP) come in pairs each giving
the correlated lists of positions and variants within the block specified
in the immediately preceeding line.  Each list is encoded as an initial integer giving the length of the following sequence of integers.


**Ground Truth / Reads.**
If the -e option is given then the source of each simulated read and the errors
introduced into it, are recorded in a [ONE-code](http:url) file of the name <code>\<name>.err</code> where <code>\<name></code> is either the
path of the -o option if present, or the root of the genome fasta file.

The ONE-code file consists of a sequence of lines the first of which is
always <code>1 3 err</code> indicating that this file is a ".err" ONE-code
file.  The file then lists the origin of each read in a haplotype and the
errors introduced into it in a sequence of triples of 'S', 'L', and 'O' lines
in exactly the same order the reads occur in the read sequence fasta file:

```
<hapreads> <- || 'h' <reads:int> <rate:real> || <reads>+
<reads>    <- || 'S' <contig:int> <orient:[0-1]> <start:int> <end:int> ||
              || 'O' <int:n> <string of length n over [=IDShztHZT]> ||
              || 'L' <int:2n+1> <len_1:int> ... <len_2n+1:int> ||
```

The tokens between || are individual lines and it is guaranteed that each
line begins with the given character and tokens are separated by a single blank.
The -h line (h for haplotype) indicates that the ensuing read descriptions are
for the next haplotype and gives the number of these reads and the rate of
divergences of the haplotype from the original source sequence.
The -S lines (S for source) indicate which contig, orientation, and interval
of the haplotype the read was sampled from.
Together, the -O (O for operation) and -L (L for length) lines encode an edit
script for converting the sequence in the haplotype interval into the read.
If the length of the -O string is n then the -L array has length 2n+1.
The i'th operator symbol is paired with the 2i+1'st length.  The 2i'th length
value gives the number of equal aligned bases between the effect of the i-1'st
and i'th operators.  The operater/length pairs are interpreted as follows:

```
S,#   substitute the next base for # where 0 => A, 1 => C, 2 => G, and 3 => T
I,#   insert the base # where 0 => A, 1 => C, 2 => G, and 3 => T
D,0   delete the next base (the # value of 0 contains no information)
```

The error model treats insertions and deletions at the end of microsatellites of periodicity 1, 2, and 3 as special cases,
as the error rates for these operations are unusually elevated. 
To reflect this, we extend the I and D alignment operators with the additional operators h, z, t, and H, Z, T indicating deletions and insertions at the end of homopolymer, di-nucleotide, and tri-nucleotide sequences, respectively.

```
h,#	  delete the next # symbols (of a homopolymer string)
z,#	  delete the next # symbols (of a dinucleotide string)
t,#	  delete the next # symbols (of a trinucleotide string)
H,#  insert the previous base # times (to extend a homopolymer string)
Z,#  insert the dinucleotide pair at the current point, adding # bases
T,#  insert the trinucleotide triple at the current point, adding # bases
```

**In order to facilitate evaluation and analysis, the reads are guaranteed
to be output in order of haplotype, then contig, and then start position in the
contig/scaffold.**

**Miscellaneous Options.**
If the -v option is set, then HIsim outputs a summary of the output it
produces.

If the -r option is set, then HIsim seeds the random number
generator with the supplied integer.  This is useful if you want HIsim
to produce the same data set with each call.  Ordinarily, the random number
generator is seeded with the process id of the call, implying a different
pseudo-random number sequence with each invocation.

The -C option asks HIsim to consider each contig of the source to be
a circular molecule.  This is useful if one wishes to avoid the lower
coverage that inevitably occurs at contig tips.

The -w option controls how many base pairs per line to output in the
Fasta read file, and the -U option specifies that upper-case letters should
be used (default is lower-case).

<a name="hifasta"></a>

```
3. HIfasta <genome>[.fast[aq]] <truth:.err+.hap#>
```

Given a source genome fasta file \<genome> and the root path \<truth> to the ground truth
generated by HIsim for a data set derived from this genome, `HIfast` will output
a fasta file of the simulated reads onto the standard output.
The ground truth files \<truth>.err and \<truth>.hap# for # = 1, 2, ..., occupy
2-5% of the space of the simulated data set, so a convenient way to record or pass a data
set on is to simply pass the smaller ground truth files.  Of course the tersest but
more time consuming method for recreating a data set is to save the source, the model
produced by HImodel (about 1MB), and the exact call to HIsim including an explicit
setting of the random number generator seed with -r.  The tradeoff is that regenerating
the data set from scratch takes about 3-4 times as long compared to generating it with
`HIfasta`.     
                     
<a name="error-model"></a>

## The Error Model

The error model provides a probability of an error, either insertion, deletion or substitution at each position in a read sequence.
For a given location, the model considers the 7-mer consisting of the
base in question, the 3 bases prior to it, and the 3 bases after it.
For every such 7-mer, HImodel, computes the frequency of each error type
in all the contexts where the 7-mer occurs at the center of a k-mer in
a FastK table.

So the core model consists of a table over all 7-mers, that for each
gives the measured likelihood of a deletion, or an insertion or substitution of a particular base.  However it is well known that for homopolymer runs, the longer the run, the higher the
likelihood of an insertion or deletion of the homopolymer base -- well
above the rate that would be observed for the 7-mer at the end of the
homopolymer.  A similar elevation is seen at the end of di-nucleotide
and tri-nucleotide satellites albeit the effect becomes smaller with
unit length.  So our error model explicitly determines the error rates
of insertions and deletions at the ends of micro-satellite runs of periodicity 1, 2, and 3 for each distinct micro-satellite unit as the
rate varies with the bases involved.  So for each microsatellite unit
(4 for homopolymers, 20 for di-nucleotides, and 60 for tri-nucleoties)
and each microsatellite unit length (up to some maximum limited by the
k-mer size of the analyzed data set), the frequency with which insertions
or deletions of satellite bases is recorded and is a part of the error
model.

So in overview, at any given location in a sequence, if the location is the end of a
microsatellite then those error likelihoods are considered first, and then
those of the general 7-mer model that do not overlap with the micro-satellite cases.  In general, the tables and their use is arranged so
that any possible edit is considered in only one table context, but we
will not burden the reader with the details of how we treat such
redundancies.

The model also records the distribution of read lengths in the data set and
HIsim uses this distribution to generate reads with a distribution of the same
shape, but possibly offset and scaled to have a user-specified mean and standard
deviation.  Lastly, the distribution of the overall error of every read in the data set is recorded as a function of length (to the nearest 1Kbp).  These
distributions are then used by the simulator to either increase or decrease
the base level error model rate to give reads that have higher or lower
overall error as seen in realistic data sets.
                     
<a name="C-lib"></a>

## Ground Truth C-Library

One can regenerate an entire simulated read data set, given the source genome and the haplotype and read ground truth files produced with the -h and -e options of HIsim.  This of course should be true, as the ground truth specifies exactly how each read was obtained from each haplotype, and how each haplotype was obtained from the source genome.  To assist one in leveraging this ground-truth in the testing and evaluation of assemblies, the C-library `lib_sim.[ch]` is included in the package and is described in what follows organized around several principal data types.  The library routines are all re-entrant requiring
one to optionally provide pointers to work data (e.g. `Get_Source`).

```
typedef void Genome;

Genome *Load_Genome(char *name);
void    Free_Genome(Genome *gene);
int64   Size_Of_Genome(Genome *gene);
void    Print_Genome(Genome *gene, FILE *file);

int64   Genome_Length(Genome *gene);
int     Scaffold_Count(Genome *gene);
int64   Scaffold_Length(Genome *gene, int i);
char   *Scaffold_Sequence(Genome *gene, int i);
```

To begin one needs an encoding of the source genome.  `Load_Genome` creates an in-memory `Genome` object given the name of its .fasta or .fastq file.  One can query the memory size of the encoding with `Size_Of_Genome`, print it out as a sequence of contigs with `Print_Genome`, and free the object with `Free_Genome`.

For a given Genome object one can query its total length in base pairs with
`Genome_Length`, get the number of scaffolds in the genome with `Scaffold_Count`, get the length of a given scaffold with `Scaffold_Length`,
and get the dna sequence of a given scaffold with `Scaffold_Sequence`.
Note carefully, that runs of N's denoting gaps between contigs, if present,
are removed from the scaffolds by `Load_Genome`, just as is done by **HIsim**.

```
typedef void HapTruth;

HapTruth  *Load_HapTruth(FILE *file, Genome *gene);
void       Free_HapTruth(HapTruth *hap);
int64      Size_Of_HapTruth(HapTruth *hap);
void       Print_HapTruth(HapTruth *hap, FILE *file);
Genome    *Haplotype_Genome(HapTruth *hap);
```

A `HapTruth` object encodes the ground truth needed to generate a haplotype from a source genome.  `Load_HapTruth` creates such an object from a ground truth file output by HIsim, where the source genome must be supplied as an argument.  The ground truth for each haplotype is written to a separate file by HIsim, so each haplotype must be created with individual calls.   One can query the memory size of the encoding with `Size_Of_HapTruth`, print it out as a series of blocks with embedded SNPs with `Print_HapTruth`, and free the object with `Free_HapTruth`.  Lastly, one can produce an explicit `Genome` object for the haplotype with `Haplotype_Genome`.

```
typedef struct
  { int    nreads;   //  number of reads in the data set
    float *rate;     //  rate[i] = divergence of haplotype i from the source
    void  *hidden[2];
  } ReadTruth;

ReadTruth *Load_ReadTruth(FILE *file, int nhaps, HapTruth **haps);
int64      Size_Of_ReadTruth(ReadTruth *ReadTruth);
void       Free_ReadTruth(ReadTruth *ReadTruth);
int        Print_Fasta(ReadTruth *reads, FILE *file);
```

A `ReadTruth` object encodes all the ground truth about a read data set and is
created by calling `Load_ReadTruth` with the ground truth file, and an array `haps`
of the `nhaps` haplotype ground truth objects describing the haplotypes the reads were sample from.  One can query the memory size of the encoding with `Size_Of_ReadTruth`.  `Free_ReadTruth` frees all the memory associated with such an object.  Lastly, `Print_Fasta` will recreate the
simulated reads in fasta format on the given file.


```
typedef struct           //  Read was sample from:
  { HapTruth  *hap;      //    this haplotype 
    int        hidx;     //    index of haplotype (starting at 1)
    int        contig;   //    this contig of the haplotype (starting at 1)
    int        orient;   //    in the forward (0) / reverse (1) orientation
    int64      beg;      //    the interval [beg,end] of the contig
    int64      end;
  } Source;

Source *Get_Source(ReadTruth *r, int64 i, Source *src);
void    Free_Source(Source *source);
```

One can get the information about how the read was sampled from the haplotype by calling `Get_Source`.  It fills in the object pointed at by `src` if not NULL, otherwise it allocates an object.  In both cases it returns a pointer to the filled in object.  An allocated object can later be freed with `Free_Source`.  The position `beg` and `end` are conceptually
between base pairs starting at 0, e.g. [0,3] is the first 3 bases of a contig.

```
typedef void Slice;

Slice  *Get_Slice(Source *source, Slice *slice);
void    Free_Slice(Slice *slice);
void    Print_Slice(Slice *slice, FILE *file);
int     Slice_Length(Slice *slice);
int     Snps_In_Slice(Slice *slice);
uint8  *Slice_Sequence(Slice *slice, uint8 *seq);
uint8  *True_Sequence(ReadTruth *r, int64 i, uint8 *seq);
```

Each read was sampled from an interval of a haplotype which in turn is a series of SNP mutated blocks/intervals of the source genome.  Therefore each read can
be viewed as a series of SNP mutated intervals of the source genome.
An `Slice` encodes this representation of a read and is produced by calling
`Get_Slice` with a read's source information.  Like `Get_Source` the routine either fills in the slice provided, or if `slice` is NULL then it allocates one, in either case returning a pointer to the filled in object.
An `Slice` object is
freed with `Free_Slice` and printed out with `Print_Slice`.  `Slice_Length` returns the length, in base pairs, of the slice which is also the length of
the read before the introduction of any sequencing errors.  `Snps_In_Slice` returns the number of SNPs introduced into the intervals constituting the slice.
Finally, `Slice_Sequence` reconstructs the sequence of the slice in the argument `seq` if it is not NULL, otherwise allocating adequate memory for the sequence.  It returns a pointer to the constructed sequence.  Please note that this sequence is encoded as a uint8 array of values from 0-3 representing a, c, g, and t, respectively.  `True_Sequence` is analogous to `Slice_Sequence` but does so directly for the i'th read in the data set, bypassing the need to produce a source and slice object.  Note carefully that it further
orients the read according to which strand it was sampled from, whereas `Slice_Sequence` does not.

```
typedef void Edit;

Edit   *Get_Edit(ReadTruth *r, int64 i, Edit *edit);
void    Free_Edit(Edit *edit);
void    Print_Edit(Edit *edit, FILE *file);
int     Edit_Length(Edit *edit);
int     Errors_In_Edit(Edit *edit);
int     Error_Type_In_Edit(Edit *edit, char kind);
uint8  *Edit_Sequence(Edit *edit, uint8 *in, uint8 *out);
uint8  *Read_Sequence(ReadTruth *r, int64 i, uint8 *in, uint8 *out);
```

Each read had a series of errors introduced into it which are recorded in the ground truth as an editing script that edits the true sequence from the haplotype into the final read.  Such an editing script is encoded in a `Edit` object and is produced for the i'th read in a data set with `Get_Edit`.
Like `Get_Source` the routine either fills in the edit provided, or if `edit` is NULL then it allocates one, in either case returning a pointer to the filled in object.
A `Edit` object is free with `Free_Edit` and a a representation of it displayed with `Print_Edit`.  The length of the edited/final read is returned by `Edit_Length`, and the number of base pairs of error introduced by a script by `Errors_In_Edit`.  `Error_Type_In_Edit` returns the number of base pairs of error of type `kind` in the
script where kind is the letter I, S, D, H, Z, T, h, z, or t.
The final, error-laden read can be produced by calling `Edit_Sequence` with the
edit script `edit` for the read and its true sequence `in` (as produced by say `True_Sequence` descript above).  The final result is place in `out` if it
is not NULL, otherwise space is allocated for it and the routine returns a poiter to the result in both cases.
One can mored directly produce a read's sequence by calling `Read_Sequence`
which takes buffers `in` and `out` for the true and final sequences.  If either of these is NULL, then the necessary space is allocated.  Note carefully that if `in` is NULL, then the space for the true sequence is allocated and then immediately freed after production of the final sequence.  It is thus more efficient to supply this buffer in a production code.
Finally, note carefully that `Read_Sequence` 
orients the read according to which strand it was sampled from, whereas `Edit_Sequence`
operates on the supplied input sequence assuming it is in the correct orientation.

```
typedef enum
          { NO_ALIGNMENT = 0, SAME_HAP_OVL = 1, DIFF_HAP_OVL = 2, DIFF_HAP_LA = 3 }
        Relation;

typedef struct
  { Relation  type;       //  alignment type
    int       r1, r2;     //  between reads r1 & r2
    int       b1, e1;     //  if type != NO_ALIGNMENT then
    int       b2, e2;     //    [b1,e1] of the first read aligns to [b2,e2] of the second
    Edit     *edit;       //  if requested, implied edit script when reads overlap
    int       nsvs;       //  # of differences due to structural variations between hap's
    int       nsnp;       //  # of differences due to SNPs between hap's
    int       nerr;       //  # of differences due to errors in the reads
    void     *hidden[4];
  } Alignment;

Alignment *Align(ReadTruth *r, int64 i, int64 j, int min_match, int do_edit, Alignment *align);
void       Free_Alignment_Edit(Alignment *align);
void       Print_Alignment(Alignment *align, FILE *file);
```

`Align` determines if, based on the ground truth alone, there is a "true" overlap between
reads `i` and `j` involving at least `min_match` base pairs excluding SNP variation and
read errors.
As per convention adopted for this library, the routine fills in its answer in the input
parameter `align` if it is non-NULL, and otherwise allocates one, and in either case returns
a pointer to the object filled in.
The `Alignment` record gives the pair of read indices, `r1` and `r2` compared and the `type` of relationship
found between the reads where there are four possibilities: `NO_ALIGNMENT` implying there
is not a true overlap between them, `SAME_HAP_OVL` if they both come from the same haplotype and the intervals they were sampled from intersect, and `DIFF_HAP_OVL` or `DIFF_HAP_LA` if
they come from different haplatypes, but there is enough synteny between the haplotypes in
the intervals from which the reads were sampled, that more than `min_match` bases are
expected to align.  `DIFF_HAP_OVL` is returned if the matched segments form a proper overlap
and `DIF_HAP_LA` (where LA stands for Local Alignment) otherwise.  The intervals, `[b1,e1]`
and `[b2,e2]`, of the
two reads that are expected to align is also determined when there is an alignment
relationship between the reads.

If `do_edit` is non-zero, then `Align` further finds an edit script that converts the
aligned interval of the first read into the aligned interval of the second read and places
it in `edit`.  In addition the number of mismatches due to SV's and SNP's between the haplotypes, and due to error in the reads is recorded in the fields `nsvs`, `nsnp`, and
`nerr`, respectively.  The edit script can be operated on as a proper `Edit` object but
there are two important differences as follows.

The first is that the edit operations are `s:#`, `d:#`, and `i:#`, and `x:1` where s, d, and i are substitutions, deletions, and insertions of `#` bases,
and x specifically denotes a one base substitution induces by a SNP variation
in one or the other haplotype in the event the reads do not come from the same haplotype.
These edit operations are needed in this context as gaps introduced by structural variations
can be large and the actual bases of the reads are not considered as the basis for matching,
but rather the simulation ground truth.  All the routines for `Edit` objects work on these
"alignment" edit scripts as well.

The second difference with alignment edits is that they are allocated as part of the
storage of an `Alignment` object.  Indeed the memory for the edit script is cached with
the Alignment object and recycled if it is given repeated as the input parameter
to `Align` and expanded when necessary.  So one needs to explicitly free this memory, especially if leaving the scope of its containing alignment object if it is stack allocated.  To do so call `Free_Alignment_Edit`.

Finally, `Print_Alignment` simply outputs the type and overlap intervals and if an
edit script is present, then that too.

```
typedef void *Error_Model;

Error_Model *Load_Error_Model(char *name);
float       *Read_Error_Profile(Error_Model *epro, uint8 *read, int len,
                                float *profile);
```

HIsim uses an error model produced by HImodel.  This model can be loaded from a source `.model`
file with the routine `Load_Error_Model`.  The model source file, by the way, is always
exactly 933,672 bytes or about 1MB.  Given the error model, one can then scan a read sequence with `Read_Error_Profile` which returns a real-valued array of the same length
as the read that contains the *relative* probability of an error at each position along the read.  That is the true probabilites are the given numbers scaled by some fixed factor reflective of the actual error rate of the particular read in question.