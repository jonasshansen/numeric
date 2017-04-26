CFLAGS = -std=gnu99 -Wall 
CFLAGS += $$(gsl-config --cflags)
CC = gcc
LDLIBS = -lm
LDLIBS += $$(gsl-config --libs)
OBJS = main.o linterp.o qspline.o
FORMAT = png
SCRIPT = log.plot


plot.$(FORMAT): out.txt Makefile
	graph -T $(FORMAT) out.txt > plot.$(FORMAT)


out.txt: main
	./$^ > out.txt
main: $(OBJS) qspline.h
	$(CC) $(CFLAGS) $(OBJS) $(LDLIBS) -o $@
main.o: main.c qspline.h


linterp.o: linterp.c

qspline.o: qspline.c

.PHONY: clean
clean:	
	rm -rf $(OBJS) main *.txt *.png
