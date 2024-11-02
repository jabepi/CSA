CC = g++
CFLAGS=-std=c++20 -O3 -march=native -Wall -m64 -g -DCLOCK=CLOCK_REALTIME -fpermissive
OMPFLAGS=-fopenmp
LDFLAGS=-g
LIBS=-lm -lrt

.PHONY: default
default: cgp3d.x

# === Code

cgp3d.x: cgp3d.o
	$(CC) -o $@ $^ $(LFLAGS) $(LIBS) $(OMPFLAGS)

.cc.o:
	$(CC) $(CFLAGS) $(INCLUDE) $(OMPFLAGS) -c $<

.PHONY: clean
clean:
	rm -f *.o *~ cgp3d.x
