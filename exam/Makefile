CFLAGS = -Wall -std=gnu11
LDLIBS = -lm

all: out.txt
	cat $<  # Print the content of out.txt to terminal

out.txt: main
	./$< > $@  # Run main and write its output to out.txt

clean:
	rm main out.*
