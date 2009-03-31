OPENMP=no
MPI=no
FFTW3_MPI=no
DEBUG=no

CLUSTER="SE Linux x86 (32 bits)"
	CC=icc						#No MPI yet...
	FFTPATH=-L/home/raid/chaos/gl293/usr/lib
	CFLAGS=-O3 -c99
	OPENMP_FLAG=-openmp
ifeq ($(DEBUG),yes)
	CFLAGS=-g -DDEBUG
	LDFLAGS+=-g
endif
