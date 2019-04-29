# DIRECTORIES
INSTALLDIR = /usr/local
TESTDIR = test
BINDIR = bin
SRCDIR = src
DOCDIR = doc
BUILDDIR = doc/_build

# PROGRAMS
CC = gcc # compiler
SPHINXBUILD = sphinx-build # documentation

# COMPILING OPTIONS 
CFLAGS = -Wall -Werror -Ofast # show warnings as erros, optimize with fast-math
INCLUDES = -I/usr/include -Iinclude # header directories
LFLAGS = -L/usr/lib # library directories
LIBS = -lgsl -lgslcblas -lm -lflint-arb # gsl, gscblas, math, arb (for cylinder)
CHECKLIBS = -lcheck -lpthread -lrt -lsubunit # libs needed for check

# FILES
SRC = $(wildcard $(SRCDIR)/*c) 
OBJS = ${SRC:.c=.o}
TESTOBJS = $(SRCDIR)/plate.o $(SRCDIR)/qfhelp.o $(TESTDIR)/check_quaca.o 
PROG = $(BINDIR)/quaca
TEST = $(BINDIR)/check_quaca


# show help if you just run 'make'
help:
	@echo "Please use 'make \033[1mtarget\033[0m' where \033[1mtarget\033[0m can be"
	@echo " \033[1mquaca\033[0m    to make quaca executable"
	@echo " \033[1mtest\033[0m     to make test executable"
	@echo " \033[1mdoc\033[0m      to make documentation (html)"
	@echo " \033[1minstall\033[0m  to install quaca executable to $(INSTALLDIR)"
	@echo " \033[1mclean\033[0m    to remove all .o files"

# quaca executable
quaca: $(PROG)

# link .o's to program (quaca)
$(PROG): $(OBJS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ $^ $(LFLAGS) $(LIBS)  

# make every .c to .o
.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

# TEST
$(TEST): $(TESTOBJS) 
	$(CC) $^ -o $@ $(LFLAGS) $(LIBS) $(CHECKLIBS) 

# make test
test: $(TEST)

# make documentation 
doc:
	$(SPHINXBUILD) -M html $(DOCDIR) $(BUILDDIR)

# install binary
install: $(PROG)
	@echo Installing executable to $(INSTALLDIR)/bin
	mkdir -p $(INSTALLDIR)/bin
	cp -f $^ $(INSTALLDIR)/bin
	chmod 755 $(INSTALLDIR)/bin/quaca

# clean object files
clean:
	rm -f $(SRCDIR)/*.o
	rm -f $(SRCDIR)/*.o

.PHONY: help quaca test doc clean
