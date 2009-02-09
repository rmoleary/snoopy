


all: Cflow

Cflow: error.o Cflow.o mainloop.o initflow.o timestep.o common.o output.o interface.o
	$(CC) -o Cflow $(CFLAGS) error.o Cflow.o mainloop.o initflow.o timestep.o common.o output.o interface.o -L$(FFTPATH) $(FFTLIBRARIES)

Cflow.o: gvars.h common.h Cflow.c
	$(CC) Cflow.c $(CFLAGS) -c
    
error.o: error.c
	$(CC) error.c $(CFLAGS) -c

timestep.o: gvars.h common.h timestep.c
	$(CC) timestep.c $(CFLAGS) -c

mainloop.o: gvars.h common.h mainloop.c
	$(CC) mainloop.c $(CFLAGS) -c

initflow.o: gvars.h common.h initflow.c
	$(CC) initflow.c $(CFLAGS) -c

output.o: gvars.h common.h output.c
	$(CC) output.c $(CFLAGS) -c

common.o: gvars.h common.c
	$(CC) common.c $(CFLAGS) -c

interface.o: gvars.h interface.c
	$(CC) interface.c $(CFLAGS) -c
	
clean:
	rm -f *.o *.mod
