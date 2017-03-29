CFLAGS = -std=gnu99 -Wall

.PHONY: default
default: hello
	./hello

hello: hello.c
	$(CC) $(CFLAGS) $< -o $@

.PHONY: clean
clean:
	rm hello 
