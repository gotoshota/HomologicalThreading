FC = gfortran
FFLAGS = -O3 -fopenmp

all: main

main: $(obj)
	FC=$(FC) f2py --f90flags="$(FFLAGS)" -m compute -c compute.f90 -lgomp --quiet

clean:
	rm -f *.so *.o *.mod

.PHONY: all clean

