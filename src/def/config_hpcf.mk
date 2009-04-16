OPENMP=no
MPI=yes
FFTW3_MPI=yes
DEBUG=no

CLUSTER="HPCF Cambridge"
ifeq ($(MPI),yes)
	CC=mpicc
else
	CC=icc
endif
FFTPATH=-L/home/gl293/usr/lib
CFLAGS=-xHOST -O3 -ipo -no-prec-div -I/home/gl293/usr/include
OPENMP_FLAG=-openmp
ifeq ($(DEBUG),yes)
	CFLAGS=-g -DDEBUG -I/home/gl293/usr/include
	LDFLAGS+=-g
endif
