##### path names #####
# header files, source files, middle files and data files
# relative to Makefile

IDIR = include
SDIR = src
ODIR = obj
DATADIR = ../data


##### compile & link options #####

OUTPUT = idm_whn

CC = icc
LIBS = -lfftw3 -lm -qopenmp
CFLAGS = -I $(IDIR) $(LIBS)


##### file lists #####

_DEPS = Basic.h Para.h Field.h Statis.h Matrix.h Interp.h IDM.h SGS.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJS = Bulk.o Mesh.o Scla.o Vctr.o Field.o Matrix.o Interp.o IDM.o SGS.o Statis.o Para.o main.o
OBJS = $(patsubst %,$(ODIR)/%,$(_OBJS))


#### targets & commands #####
# note: $@ target, $^ all dependents, $< the first dependent

$(OUTPUT): $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

$(ODIR)/%.o: $(SDIR)/%.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

dir:
	mkdir -p $(DATADIR)/fielddata
	mkdir -p $(DATADIR)/probedata
	mkdir -p $(DATADIR)/statdata
	mkdir -p $(DATADIR)/postdata

dirclean:
	rm -r $(DATADIR)/fielddata
	rm -r $(DATADIR)/probedata
	rm -r $(DATADIR)/statdata
	rm -r $(DATADIR)/postdata
	make dir

.PHONY: remake
remake:
	make clean
	make
	mv $(OUTPUT) $(DATADIR)/

.PHONY: clean
clean:
	rm -f $(OUTPUT) *~ $(ODIR)/*.o $(INCDIR)/*~
