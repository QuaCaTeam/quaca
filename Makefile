#
# Makefile 
#

# options
CC = gcc
CFLAGS = -Wall -O2 # show warnings 
INCLUDES = -I/usr/include -Isrc/h # include headers
LFLAGS = -L/usr/lib
LIBS = -lgsl -lgslcblas -lm -lflint-arb -lmpfr # libraries (gsl and arb), on arch -larb on debian -lflint-arb

# file
SRC = $(wildcard src/*c) 
OBJS = ${SRC:.c=.o}
PROG = bin/qfnum

all: $(PROG) clean

$(PROG): $(OBJS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $(PROG) $(OBJS) $(LFLAGS) $(LIBS)

.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

# make docs
docs:
	@doxygen doc/doxconf

clean:
	rm src/*.o

.PHONY: all clean tester
