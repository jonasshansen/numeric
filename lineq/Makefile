CFLAGS = -std=gnu99 -Wall
CFLAGS += $$(gsl-config --cflags)
CC = gcc
LDLIBS = -lm
LDLIBS += $$(gsl-config --libs)
OBJS = main.o qr_gs_decomp.o

all:out.txt
	cat $<
out.txt:main
	./$< > $@
main: $(OBJS)
$(OBJS):qr_gs_decomp.h
clean:
	rm *.o main out*
