CFLAGS = -O -march=native -std=gnu99 -Wall
CFLAGS += $$(gsl-config --cflags)
CFLAGS += -I../lineq/
LDLIBS += $$(gsl-config --libs)
OBJS = main.o lsfit.o qr_gs_decomp.o

all: plot.svg

plot.svg: out.txt Makefile
	echo '\
	set term svg background rgb "white";\
	set out "$@";\
	set xlabel "x";\
	set ylabel "y";\
	set title "least-squares fit of F_c(x)=c_0/x+c_1+c_2*x to some data" ;\
	plot "$<" index 0 with errorbars title "data"\
	,"$<" index 1 using 1:2 with lines title "F_c(x)"\
	,"$<" index 1 using 1:3 with lines title "F_{c+dc}(x)"\
	,"$<" index 1 using 1:4 with lines title "F_{c-dc}(x)"\
	'|gnuplot

out.txt: main
	./$< > $@

main: $(OBJS)

lsfit.o: ../lineq/qr_gs_decomp.h

qr_gs_decomp.o: ../lineq/qr_gs_decomp.c ../lineq/qr_gs_decomp.h
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm main *.o out.* plot.svg
