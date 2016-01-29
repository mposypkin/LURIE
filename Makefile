# To make run
# make dep
# make

# C++ compiler
CPP = c++

# BNB solver root
BNBROOT = ../BNB-solver

# include folders
INCLUDE = $(BNBROOT)

# C++ compiler options
CPPOPTS = --std=c++11 -I$(INCLUDE)

# Libraries to inculde 
LIBS = $(BNBROOT)/util/common/util.a $(BNBROOT)/problems/optlib/optlib.a $(BNBROOT)/libjson/libjson.a 

# Linkers flags
LDFLAGS = -pthread 

# tests
TESTS = testenergy.exe 

all: testenergy.exe testmatmodel.exe searchmbh.exe searchgdsc.exe locmincheck.exe genpoints.exe searchmcmbh.exe montecarlo.exe searchall.exe

-include deps.inc

clean: 
	rm -f *.exe *.o *.a *~ *.log deps.inc

dep:
	$(CPP) $(CPPOPTS) -MM -c *.cpp >> deps.inc;\
        true
tests:
	@for i in $(TESTS); do if ./$$i > /dev/null; then echo TEST PASSED; continue; else echo TEST FAILED; fi done


.o.exe:
	$(CPP) $(OPTS) -o $@ $< $(LIBS) $(EXTLIBS) $(LDFLAGS)

.cpp.o:
	$(CPP) $(CPPOPTS) -c $<

.c.o:
	$(CC) $(COPTS) -c $<

.SUFFIXES:
.SUFFIXES: .o .a .cpp .c .exe

