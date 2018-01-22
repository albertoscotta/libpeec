CC = gcc
CC_OPTS = -fopenmp
LIBS = -lpeec -lm
EXAMPLES = example_0

all: $(EXAMPLES) Makefile

%: %.c
	$(CC) $< -o $@ $(CC_OPTS) $(LIBS)

clean:
	rm -f $(EXAMPLES)

.PHONY: all clean
