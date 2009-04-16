OPENMP=yes
MPI=no
FFTW3_MPI=no
DEBUG=no

CLUSTER="SE Linux x86 (64 bits)"
CC=icc
FFTPATH=-L/home/raid/chaos/gl293/usr/lib
CFLAGS=-O3 -xHOST -ipo -no-prec-div -static
OPENMP_FLAG=-openmp
ifeq ($(DEBUG),yes)
	CFLAGS=-g -DDEBUG
	LDFLAGS+=-g
endif
