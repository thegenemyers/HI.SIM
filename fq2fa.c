#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <ctype.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "gene_core.h"

//  Read next line into a buffer and return a pointer to the buffer and set *plen
//    the length of the line.  NB: replaces '\n' with '\0'.

static char *read_line(FILE *input, int *plen)
{ static char *buffer;
  static int   bmax = 0;
  int len;

  if (bmax == 0)
    { bmax = 500;
      buffer = (char *) Malloc(bmax,"Allocating read buffer");
      if (buffer == NULL)
        exit (1);
    }

  if (fgets(buffer,bmax,input) == NULL)
    return (NULL);

  len = strlen(buffer);
  while (buffer[len-1] != '\n')
    { bmax = ((int) (1.4*bmax)) + 100;
      buffer = (char *) Realloc(buffer,bmax,"Reallocating read buffer");
      if (buffer == NULL)
        exit (1);
      if (fgets(buffer+len,bmax-len,input) == NULL)
        { fprintf(stdout,"%s: Last line of file does not end with new-line\n",Prog_Name);
          exit (1);
        }
      len += strlen(buffer+len);
    }
  buffer[--len] = '\0';

  if (plen != NULL)
    *plen = len;
  return (buffer);
}

//  Read .fastq data set from input output fasta for each entry

static void fastq2fasta(FILE *input)
{ char  *line;
  int    nread, maxlen;
  int    len, qlen;
  int64  totbp;

  nread  = 0;
  totbp  = 0;
  maxlen = 0;
  while ((line = read_line(input,NULL)) != NULL)
    { if (line[0] != '@')
        { fprintf(stderr,"%s: Entry header does not start with an @-sign\n",Prog_Name);
          exit (1);
        }

      fprintf(stdout,">%s\n",line+1);

      line = read_line(input,&len);
      if (line == NULL)
        { fprintf(stderr,"%s: Incomplete fastq entry\n",Prog_Name);
          exit (1);
        }

      fprintf(stdout,"%s\n",line);

      nread += 1;
      totbp += len+1;
      if (len > maxlen)
        maxlen = len;

      line = read_line(input,NULL);
      if (line == NULL)
        { fprintf(stderr,"%s: Incomplete fastq entry\n",Prog_Name);
          exit (1);
        }
      if (line[0] != '+')
        { fprintf(stderr,"%s: Divider line does not start with a +-sign\n",Prog_Name);
          exit (1);
        }
      line = read_line(input,&qlen);
      if (line == NULL)
        { fprintf(stderr,"%s: Incomplete fastq entry\n",Prog_Name);
          exit (1);
        }
      if (len != qlen)
        { fprintf(stderr,"%s: QV line does not have the same length as sequence line\n",Prog_Name);
          exit (1);
        }
    }

  fprintf(stderr,"\n   # Reads = ");
  Print_Number((int64) nread,0,stderr);
  fprintf(stderr,"\n   Max len = ");
  Print_Number((int64) maxlen,0,stderr);
  fprintf(stderr,"\n   Total   = ");
  Print_Number(totbp,0,stderr);
  fprintf(stderr,"\n");
}

int main(int argc, char *argv[])
{ (void) argv;

  if (argc != 1)
    { fprintf(stderr,"Usage: fq2fa is a pipe\n");
      exit (1);
    }

  //  Open input and output files

  fastq2fasta(stdin);

  exit (0);
}
