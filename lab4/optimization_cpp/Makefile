include Makefile.in

.PHONY: default
default: cgp3d.x

# === Code

cgp3d.x: cgp3d.o
	$(CC) -o $@ $^ $(LFLAGS) $(LIBS) $(OMPFLAGS)

.c.o:
	$(CC) $(CFLAGS) $(INCLUDE) $(OMPFLAGS) -c $<

.PHONY: clean
clean:
	rm -f *.o *~ cgp3d.x
