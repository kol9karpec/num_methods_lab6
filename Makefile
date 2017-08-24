CC=gcc
CFLAGS=-Wall -Werror -I/home/kol9karpec/gsl/include
LDFLAGS=-L/home/kol9karpec/gsl/lib
LDLIBS=-lgsl -lgslcblas -lm

main: main.o

clean: 
	@rm -f *.o

.PHONY: clean
