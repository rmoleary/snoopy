OBJ= error.o snoopy.o mainloop.o initflow.o timestep.o common.o output.o interface.o gfft.o transpose.o debug.o shear.o read_config.o libconfig/libconfig.o libconfig/scanner.o libconfig/grammar.o

all: snoopy

snoopy: $(OBJ)
	$(CC) -o $@ $^ $(LDFLAGS)

%.o: %.c gvars.h
	$(CC) -o $@ $< $(CFLAGS) -c

clean:
	rm -f *.o libconfig/*.o snoopy

fullclean:
	rm -f *.o libconfig/*.o snoopy gvars.h

bench:
	mv ./gvars.h ./gvars.h.old
	cp ../cfg/bench/gvars.h ./gvars.h

benchclean:
	rm -f ./gvars.h
	mv ./gvars.h.old ./gvars.h
