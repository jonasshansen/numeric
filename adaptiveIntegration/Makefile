CFLAGS = -Wall -std=gnu11
LDLIBS = -lm

all:out.txt
	cat $<

out.txt:adapt
	./$< > $@

clean:
	rm adapt out.*
