# path names of header files, source files, object files and library files， relative to Makefile
IDIR = include
SDIR = src
ODIR = obj
LDIR = lib

# use g++ as the c++ compiler
CC = g++
# compile options
CFLAGS = -I $(IDIR)
# use fftw3 library and math library
LIBS = -lfftw3 -lm

_DEPS = Para.h Mesh.h Field.h Matrix.h IDM.h Statis.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = main.o Para.o Mesh.o Field.o Matrix.o IDM.o Statis.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: $(SDIR)/%.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

channel_LES: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

run: channel_LES
	./channel_LES
	make clean

dir:
	mkdir ../data/fielddata
	mkdir ../data/probedata
	mkdir ../data/statdata
	mkdir ../data/postdata

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ 





