#
# Makefile 
#

# directory to install qfnum to
PREFIX=/usr/local

# options for compiling
CC = gcc
CFLAGS = -Wall -Ofast # show warnings as erros, optimize with fast-math
INCLUDES = -I/usr/include -Isrc/h # include headers
LFLAGS = -L/usr/lib
LIBS = -lgsl -lgslcblas -lm # gsl, gscblas, math library

# files
SRC = $(wildcard src/*c) 
OBJS = ${SRC:.c=.o}
PROG = bin/quaca

# documentation
SPHINXBUILD = sphinx-build
SOURCEDIR = doc
BUILDDIR = doc/_build

all: $(PROG) clean

# link .o's to program (qfnum)
$(PROG): $(OBJS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $(PROG) $(OBJS) $(LFLAGS) $(LIBS)  

# make every .c to .o
.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

# make doc
doc:
	$(SPHINXBUILD) -M html $(SOURCEDIR) $(BUILDDIR)

# install binary
install: $(PROG) clean
	echo Installing executable to $(PREFIX)/bin
	mkdir -p $(PREFIX)/bin
	cp -f bin/qfnum $(PREFIX)/bin
	chmod 755 $(PREFIX)/bin/qfnum

# clean directory
clean:
	rm src/*.o

.PHONY: all clean tester
