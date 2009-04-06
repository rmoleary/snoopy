OPENMP=yes
MPI=no
FFTW3_MPI=no
DEBUG=no

CLUSTER="MacOS 10.5 (Leopard)"
ifeq ($(MPI),yes)
	CC=/usr/local/bin/mpicc
else
	CC=/usr/local/bin/gcc
endif
FFTPATH=-L/usr/local/lib
CFLAGS=-O3 -ffast-math -fomit-frame-pointer
OPENMP_FLAG=-fopenmp
ifeq ($(DEBUG),yes)
	CFLAGS=-g -DDEBUG
	LDFLAGS+=-g
endif
