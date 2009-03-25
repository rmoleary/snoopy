OBJ= error.o snoopy.o mainloop.o initflow.o timestep.o common.o output.o interface.o gfft.o transpose.o debug.o


all: snoopy

snoopy: $(OBJ)
	$(CC) -o $@ $(CFLAGS) $(FFTPATH) $(LDFLAGS) $^

%.o: %.c gvars.h
	$(CC) -o $@ $< $(CFLAGS) -c
	
config: gvars.h

gvars.h:
	cp def/gvars.h .

clean:
	rm -f *.o snoopy

fullclean:
	rm -f *.o snoopy gvars.h
