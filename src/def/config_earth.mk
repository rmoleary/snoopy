OPENMP=yes
MPI=no
FFTW3_MPI=no
DEBUG=no

CLUSTER="SE Linux x86 (64 bits)"
CC=icc
FFTPATH=-L/home/raid/chaos/lj272/usr/lib
CFLAGS=-O3
OPENMP_FLAG=-openmp
ifeq ($(DEBUG),yes)
	CFLAGS=-g -DDEBUG
	LDFLAGS+=-g
endif
