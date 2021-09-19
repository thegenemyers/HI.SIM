DEST_DIR = ~/bin

CFLAGS = -O3 -Wall -Wextra -Wno-unused-result -fno-strict-aliasing

CC = clang

ALL = HImodel HIsim libtest

all: $(ALL)

libfastk.c : gene_core.c
libfastk.h : gene_core.h

HImodel: HImodel.c libfastk.c libfastk.h
	$(CC) $(CFLAGS) -o HImodel HImodel.c libfastk.c -lpthread -lm

HIsim: HIsim.c gene_core.c gene_core.h
	$(CC) $(CFLAGS) -o HIsim HIsim.c gene_core.c -lpthread -lm

libtest: lib_sim.c gene_core.c gene_core.h
	$(CC) $(CFLAGS) -DTESTING -o libtest lib_sim.c gene_core.c -lpthread -lm

clean:
	rm -f $(ALL)
	rm -fr *.dSYM
	rm -f HIsim.tar.gz

install:
	cp $(ALL) $(DEST_DIR)

package:
	make clean
	tar -zcf HIsim.tar.gz LICENSE README.md Makefile *.h *.c
