#make with gsl and omp

CFLAGS=-I/usr/include/gsl
LDFLAGS=-lgsl -lgslcblas -fopenmp -lm
CC=gcc

%: %.c
	$(CC) $< -o $@ $(CFLAGS) $(LDFLAGS)

clean:
	rm -f *~ *.o core a.out
