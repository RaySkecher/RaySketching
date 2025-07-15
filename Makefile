CC=gcc
CFLAGS=-Wall -Wextra -std=c99 -O3
LDFLAGS=-lm

all: raytracer

raytracer: main.c
	$(CC) $(CFLAGS) -o raytracer main.c $(LDFLAGS)

run: raytracer
	./raytracer

clean:
	rm -f raytracer output.ppm 