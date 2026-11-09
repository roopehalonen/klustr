# KLUSTR v1.0, Fortran version.
#   make             builds the module and the test
#   make test        runs the test against ../python/test_reference.txt
#   make FC=... FFLAGS="..."            another compiler or other flags (only gfortran has been tested)
# Single precision: add -DSINGLE to FFLAGS (the constants are still computed in double precision).
FC = gfortran
FFLAGS = -O3

all: test_klustr

klustr.o: klustr.F90
	$(FC) $(FFLAGS) -c klustr.F90

test_klustr: test_klustr.F90 klustr.o
	$(FC) $(FFLAGS) -o $@ test_klustr.F90 klustr.o

test: test_klustr
	./test_klustr ../python/test_reference.txt

clean:
	rm -f *.o *.mod test_klustr

.PHONY: all test clean
