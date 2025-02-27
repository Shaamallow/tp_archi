CC = gcc
CFLAGS = -pthread -mavx -std=c11 -Wall
TARGETS = parallel_distance avx_demo

all: $(TARGETS)

parallel_distance: parallel_distance.c
	$(CC) $(CFLAGS) -o parallel_distance parallel_distance.c -lm

avx_demo: avx_demo.c
	$(CC) $(CFLAGS) -o avx_demo avx_demo.c -lm

clean:
	rm -f $(TARGETS)
