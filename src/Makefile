OBJ= error.o Cflow.o mainloop.o initflow.o timestep.o common.o output.o interface.o vtk_writer.o gfft.o transpose.o debug.o


all: Cflow

Cflow: $(OBJ)
	$(CC) -o $@ $(CFLAGS) -L$(FFTPATH) $(LDFLAGS) $^

%.o: %.c gvars.h
	$(CC) -o $@ $< $(CFLAGS) -c
	
clean:
	rm -f *.o *.mod
