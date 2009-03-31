OPENMP=no
MPI=no
FFTW3_MPI=no
DEBUG=no

CLUSTER="HPCF Cambridge"
ifeq ($(MPI),yes)
	CC=mpicc
else
	CC=icc
endif
FFTPATH=-L/home/gl293/usr/lib
CFLAGS=-O3 -I/home/gl293/usr/include
OPENMP_FLAG=-openmp
ifeq ($(DEBUG),yes)
	CFLAGS=-g -DDEBUG
	LDFLAGS+=-g
endif
