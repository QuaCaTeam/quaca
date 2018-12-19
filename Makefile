#
# Makefile 
#

# directory to install qfnum to
PREFIX=/usr/local

# options for compiling
CC = gcc
CFLAGS = -Wall -Werror -O2 # show warnings as erros, optimize level 2
INCLUDES = -I/usr/include -Isrc/h # include headers
LFLAGS = -L/usr/lib
LIBS = -lgsl -lgslcblas -lm # gsl, gscblas, math library

# files
SRC = $(wildcard src/*c) 
OBJS = ${SRC:.c=.o}
PROG = bin/qfnum

all: $(PROG) clean

# link .o's to program (qfnum)
$(PROG): $(OBJS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $(PROG) $(OBJS) $(LFLAGS) $(LIBS)  

# make every .c to .o
.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

# make doc
doc:
	cd doc && make html	&& cd ..

# install binary
install: $(PROG)
	echo Installing executable to $(PREFIX)/bin
	mkdir -p $(PREFIX)/bin
	cp -f bin/qfnum $(PREFIX)/bin
	chmod 755 $(PREFIX)/bin/qfnum

# clean directory
clean:
	rm src/*.o

.PHONY: all clean tester
