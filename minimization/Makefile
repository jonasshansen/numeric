CFLAGS = -std=gnu99 -Wall
CFLAGS += $$(gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $$(gsl-config --libs)
OBJS = main.o qnewton.o

all: out.txt
	cat $<

out.txt: main
	./$< > $@

main: $(OBJS)

clean:
	rm main *.o out.*
