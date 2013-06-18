FF	= gfortran
LD	= gfortran
CC	= g++

LDFLAGS	= `root-config --libs`
FFFLAGS = -g -c
CCFLAGS = `root-config --cflags` -g -c

all: calcpapcf

# calcpapcf: calcpapcf.o
# 	$(CC) $^ -o $@ $(LDFLAGS)

calcpapcf: calcpapcf.o  FsiTools.o FsiWeightLednicky4.o
	$(LD) $^ -o $@ $(LDFLAGS)

# calcpapcf: calcpapcf.o FsiWeightLednicky4.o
# 	$(LD) $^ -o $@ $(LDFLAGS)

%.o: %.F
	$(FF) $^ -o $@ $(FFFLAGS)

%.o: %.cxx
	$(CC) $^ -o $@ $(CCFLAGS)

clean:
	rm -f *.o calcpapcf
