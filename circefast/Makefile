CC=gcc
#CFITSIO=../cfitsio
CFITSIO=/scratch1/ydallilar/circefast/cfitsio
INCDIR=$(CFITSIO)/include
LIBDIR=$(CFITSIO)/lib
OPTS=-Wall -Wl,-rpath,$(CFITSIO)/lib -funroll-loops
LIBS=-lcfitsio -lm -lgsl -lgslcblas -lfftw3 -O3 -fopenmp

DEPS=fileio.h engine.h
OBJ=fileio.o engine.o main.o 

CFLAGS=$(LIBS) $(OPTS) -I$(INCDIR) -L$(LIBDIR) 

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) -fPIC 

all: prog shared

prog: $(OBJ)
	gcc -o main $(OBJ) $(CFLAGS)

shared: $(OBJ)
	gcc -shared -lgomp  -o circefast.so $(OBJ) $(CFLAGS)

clean: 
	rm -vf $(OBJ) main circefast.so
