CC=gcc
CFLAGS=-Wall -Werror -I/usr/include/gsl/include/
LDFLAGS=-L/usr/include/gsl/lib/
LDLIBS=-lgsl -lgslcblas -lm

main: main.o

clean: 
	@rm -f *.o

.PHONY: clean
