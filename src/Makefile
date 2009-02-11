


all: Cflow

Cflow: error.o Cflow.o mainloop.o initflow.o timestep.o common.o output.o interface.o vtk_writer.o
	$(CC) -o $@ $(CFLAGS) -L$(FFTPATH) $(FFTLIBRARIES) $^

%.o: %.c gvars.h
	$(CC) -o $@ $< $(CFLAGS) -c
	
clean:
	rm -f *.o *.mod
