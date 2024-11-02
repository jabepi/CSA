CC = g++
CFLAGS=-std=c++20 -O3 -march=native -Wall -m64 -g -DCLOCK=CLOCK_REALTIME -fpermissive
OMPFLAGS=-fopenmp
LDFLAGS=-g
LIBS=-lm -lrt

.PHONY: default
default: cgp3d.x

# === Code

cgp3d.x: cgp3d.cc
	$(CC) $(CFLAGS) $(OMPFLAGS) -o $@ $< $(LDFLAGS) $(LIBS)

.PHONY: clean
clean:
	rm -f *~ cgp3d.x
