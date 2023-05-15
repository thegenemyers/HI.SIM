DEST_DIR = ~/bin

CFLAGS = -O3 -Wall -Wextra -Wno-unused-result -fno-strict-aliasing

CC = gcc

ALL = HImodel HIsim HIreads HIhaplo libtest fq2fa

all: $(ALL)

libfastk.c : gene_core.c
libfastk.h : gene_core.h

HImodel: HImodel.c libfastk.c libfastk.h
	$(CC) $(CFLAGS) -o HImodel HImodel.c libfastk.c -lpthread -lm

HIsim: HIsim.c gene_core.c gene_core.h
	$(CC) $(CFLAGS) -o HIsim HIsim.c gene_core.c -lpthread -lm

HIreads: HIreads.c lib_sim.c lib_sim.h gene_core.c gene_core.h
	$(CC) $(CFLAGS) -o HIreads HIreads.c lib_sim.c gene_core.c -lpthread -lm

HIhaplo: HIhaplo.c lib_sim.c lib_sim.h gene_core.c gene_core.h
	$(CC) $(CFLAGS) -o HIhaplo HIhaplo.c lib_sim.c gene_core.c -lpthread -lm

libtest: lib_sim.c lib_sim.h gene_core.c gene_core.h
	$(CC) $(CFLAGS) -DTESTING -o libtest lib_sim.c gene_core.c -lpthread -lm

fq2fa: fq2fa.c gene_core.c gene_core.h
	$(CC) $(CFLAGS) -DTESTING -o fq2fa fq2fa.c gene_core.c -lpthread -lm

clean:
	rm -f $(ALL)
	rm -fr *.dSYM
	rm -f HIsim.tar.gz

install:
	cp $(ALL) $(DEST_DIR)

package:
	make clean
	tar -zcf HIsim.tar.gz LICENSE README.md Makefile *.h *.c
