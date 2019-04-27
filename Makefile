# DIRECTORIES
INSTALLDIR=/usr/local
TESTDIR=test
SOURCEDIR = doc
BUILDDIR = doc/_build

# OPTIONS 
CC = gcc
CFLAGS = -Wall -Ofast # show warnings as erros, optimize with fast-math
INCLUDES = -I/usr/include # include headers
LFLAGS = -L/usr/lib
LIBS = -lgsl -lgslcblas -lm -lflint-arb# gsl, gscblas, math library
CHECKLIBS = -lcheck -lpthread -lrt -lsubunit

# FILES
SRC = $(wildcard src/*c) 
OBJS = ${SRC:.c=.o}
PROG = bin/quaca

# documentation
SPHINXBUILD = sphinx-build


all: $(PROG)

# link .o's to program (qfnum)
$(PROG): $(OBJS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $(PROG) $(OBJS) $(LFLAGS) $(LIBS)  

# make every .c to .o
.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

# TEST
bin/check_quaca:  src/plate.o src/qfhelp.o test/check_quaca.o
	$(CC) src/plate.o src/qfhelp.o test/check_quaca.o -o bin/check_quaca $(LIBS) $(CHECKLIBS) 

test: bin/check_quaca

# make doc
doc:
	$(SPHINXBUILD) -M html $(SOURCEDIR) $(BUILDDIR)

# install binary
install: $(PROG) clean
	echo Installing executable to $(INSTALLDIR)/bin
	mkdir -p $(INSTALLDIR)/bin
	cp -f bin/quaca $(INSTALLDIR)/bin
	chmod 755 $(INSTALLDIR)/bin/quaca

# clean directory
clean:
	rm src/*.o

.PHONY: all clean
