CFLAGS = -std=gnu99 -Wall
TARGET = out.A.txt

check: $(TARGET)
	@echo "Checking target A ..."
	@echo "hello, $$(whoami)" > correct.txt
	@diff --brief correct.txt out.A.txt
	@echo "Target A seems to be ok!"
	@rm -f correct.txt

A: $(TARGET)
	cat $(TARGET)

out.A.txt: mainA
	./mainA > $(TARGET)

mainA: mainA.c # build executable file mainA from source mainA.c using the implicit rule below
# $(LINK.c) $^ $(LDLIBS) -o $@

cleanA:
	rm -rf mainA $(TARGET)
