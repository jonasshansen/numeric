CFLAGS += -Wall -std=c11 `gsl-config --cflags`
LDLIBS += `gsl-config --libs`
OBJS = main.o ode_driver.o rkstep12.o

all: out.plot.png

out.plot.png: out.txt
	graph \
		--display-type png --bitmap-size 800x600 \
		--x-label "x" --y-label "y" \
		--top-label "ode: y''=-y" \
		< $< > $@

out.txt: main
	./$< > $@

main: $(OBJS)

clean:
	rm main *.o out.*
