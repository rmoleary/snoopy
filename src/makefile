OBJ= error.o snoopy.o mainloop.o initflow.o timestep.o timestep_compressible.o common.o output/output.o output/output_dump.o output/output_spectrum.o output/output_timevar.o output/output_vtk.o output/output_common.o interface.o gfft.o transpose.o debug.o shear.o read_config.o symmetries.o boundary.o forcing.o libconfig/libconfig.o libconfig/scanner.o libconfig/grammar.o particles.o ltide.o


all: snoopy

snoopy: $(OBJ)
	$(CC) -o $@ $^ $(LDFLAGS)

%.o: %.c common.h gvars.h snoopy.h error.h debug.h gfft.h
	$(CC) -o $@ $< $(CFLAGS) -c

clean:
	rm -f *.o libconfig/*.o output/*.o snoopy

cleanall:
	rm -f *.o libconfig/*.o snoopy gvars.h

bench:
	mv ./gvars.h ./gvars.h.old
	cp problem/bench/gvars.h ./gvars.h

benchclean:
	rm -f ./gvars.h
	mv ./gvars.h.old ./gvars.h


#--ADDITIONAL DEPENDENCIES--

boundaries.o: symmetries.h gvars.h
gfft.o: transpose.h gvars.h
initflow.o: output/output_dump.h symmetries.h gvars.h
interface.o: mainloop.h output/output.h transpose.h gvars.h
mainloop.o: timestep.h output/output.h interface.h shear.h transpose.h symmetries.h boundary.h particles.h gvars.h
output/output.o: shear.h particles.h gvars.h
output/output_vtk.o: shear.h particles.h gvars.h
output/output_timevar.o: shear.h particles.h gvars.h
particles.o: shear.h transpose.h gvars.h
read_config.o: libconfig/libconfig.h gvars.h
shear.o: transpose.h gvars.h
snoopy.o: mainloop.h output/output.h initflow.h read_config.h particles.h gvars.h
timestep.o: forcing.h particles.h gvars.h
libconfig/libconfig.o: libconfig/libconfig.h libconfig/grammar.h libconfig/scanner.h libconfig/private.h gvars.h
libconfig/grammar.o: libconfig/libconfig.h libconfig/private.h gvars.h


