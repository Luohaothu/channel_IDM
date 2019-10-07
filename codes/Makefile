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

channel_WHN: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

run: channel_WHN
	./channel_WHN
	make clean

dir:
	mkdir -p ../data/fielddata
	mkdir -p ../data/probedata
	mkdir -p ../data/statdata
	mkdir -p ../data/postdata

dirclean:
	rm -r ../data/fielddata
	rm -r ../data/probedata
	rm -r ../data/statdata
	rm -r ../data/postdata
	make dir

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ 





