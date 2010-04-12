OBJ= error.o snoopy.o mainloop.o initflow.o timestep.o common.o output.o interface.o gfft.o transpose.o debug.o shear.o read_config.o symmetries.o boundary.o forcing.o libconfig/libconfig.o libconfig/scanner.o libconfig/grammar.o particles.o

all: snoopy

snoopy: $(OBJ)
	$(CC) -o $@ $^ $(LDFLAGS)

%.o: %.c common.h gvars.h error.h debug.h gfft.h
	$(CC) -o $@ $< $(CFLAGS) -c

clean:
	rm -f *.o libconfig/*.o snoopy

cleanall:
	rm -f *.o libconfig/*.o snoopy gvars.h

bench:
	mv ./gvars.h ./gvars.h.old
	cp problem/bench/gvars.h ./gvars.h

benchclean:
	rm -f ./gvars.h
	mv ./gvars.h.old ./gvars.h


#--ADDITIONAL DEPENDENCIES--

boundaries.o: symmetries.h
gfft.o: transpose.h
initflow.o: output.h symmetries.h
interface.o: mainloop.h output.h transpose.h
mainloop.o: timestep.h output.h interface.h shear.h transpose.h symmetries.h boundary.h particles.h
output.o: shear.h particles.h
particles.o: shear.h transpose.h
read_config.o: libconfig/libconfig.h
shear.o: transpose.h
snoopy.o: mainloop.h output.h initflow.h read_config.h particles.h
timestep.o: forcing.h particles.h
libconfig/libconfig.o: libconfig/libconfig.h libconfig/grammar.h libconfig/scanner.h libconfig/private.h
libconfig/grammar.o: libconfig/libconfig.h libconfig/private.h


