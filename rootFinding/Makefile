CFLAGS += -Wall -std=c11 `gsl-config --cflags`
CFLAGS += -I../lineq/
LDLIBS += `gsl-config --libs`
OBJS = main.o newton.o qr_gs_decomp.o

all:out.txt
	cat $<

out.txt: main
	./$< > $@

main: $(OBJS)

newton.o: ../lineq/qr_gs_decomp.h

qr_gs_decomp.o: ../lineq/qr_gs_decomp.c ../lineq/qr_gs_decomp.h
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm main *.o out.*
