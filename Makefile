# To make run
# make dep
# make

# C++ compiler
CPP = c++

all: testenergy.exe

clean: 
	rm -f *.exe *.o *.a *~ *.log deps.inc

dep:
	$(CPP) $(CPPOPTS) -MM -c *.cpp >> deps.inc;\
        true

.o.exe:
	$(CPP) $(OPTS) -o $@ $< $(LIBS) $(EXTLIBS) $(LDFLAGS)

.cpp.o:
	$(CPP) $(CPPOPTS) -c $<

.c.o:
	$(CC) $(COPTS) -c $<

.SUFFIXES:
.SUFFIXES: .o .a .cpp .c .exe
