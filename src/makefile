OBJ= error.o snoopy.o mainloop.o initflow.o timestep.o common.o output.o interface.o gfft.o transpose.o debug.o shear.o

all: snoopy

snoopy: $(OBJ)
	$(CC) -o $@ $^ $(LDFLAGS)

%.o: %.c gvars.h
	$(CC) -o $@ $< $(CFLAGS) -c

clean:
	rm -f *.o snoopy

fullclean:
	rm -f *.o snoopy gvars.h

bench:
	mv ./gvars.h ./gvars.h.old
	cp def/gvars_bench.h ./gvars.h

benchclean:
	rm -f ./gvars.h
	mv ./gvars.h.old ./gvars.h
