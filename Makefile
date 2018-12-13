#
# Makefile 
#

# options
CC = gcc
CFLAGS = -Wall -O3 # show warnings 
INCLUDES = -I/usr/include -Isrc/h # include headers
LFLAGS = -L/usr/lib
LIBS = -lgsl -lgslcblas -lm -larb -lmpfr # libraries (gsl and arb), on arch -larb on debian -lflint-arb

# file
SRC = $(wildcard src/*c) 
OBJS = ${SRC:.c=.o}
PROG = bin/qfnum

all: $(PROG) clean

$(PROG): $(OBJS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $(PROG) $(OBJS) $(LFLAGS) $(LIBS)

.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

# make doc
doc:
	cd doc && make html	

clean:
	rm src/*.o

.PHONY: all clean tester
