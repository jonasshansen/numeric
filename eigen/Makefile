CFLAGS = -std=gnu99 -Wall
CFLAGS += $$(gsl-config --cflags)
CC = gcc
LDLIBS = -lm
LDLIBS += $$(gsl-config --libs)
OBJS = main.o jacobi.o

all: out.txt
	cat $<

out.txt: main Makefile
	./main 4 2> $@

main: $(OBJS)

clean:
	rm main *.o out.*
