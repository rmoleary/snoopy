OPENMP=no
MPI=no
FFTW3_MPI=no
DEBUG=no

CLUSTER="Hyades Cluster"
CC=/opt/intel91039033/bin/icc		#No MPI yet...
FFTPATH=-L/home/raid/chaos/gl293/usr/lib
CFLAGS=-O3 -static
OPENMP_FLAG=-openmp
ifeq ($(DEBUG),yes)
	CFLAGS=-g -DDEBUG
	LDFLAGS+=-g
endif
