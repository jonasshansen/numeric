CFLAGS = -std=gnu99 -Wall 
CFLAGS += $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

main: main.o
	$(CC) $(CFLAGS) $< $(LDLIBS) -o $@

main.o: main.c
	$(CC) $(CFLAGS) -c $<

.PHONY: clean
clean:	
	rm -rf main.o main
