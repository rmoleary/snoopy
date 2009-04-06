OPENMP=no
MPI=yes
FFTW3_MPI=no
DEBUG=no

CLUSTER="VARGAS (IBM AIX)"
ifeq ($(MPI),yes)
	CC=mpcc_r
else
	CC=xlc_r
endif
CFLAGS=-O3 -I/usr/local/pub/FFTW/3.2/include
LDFLAGS=-lm
FFTPATH=-L/usr/local/pub/FFTW/3.2/lib
OPENMP_FLAG= -qsmp=omp
ifeq ($(DEBUG),yes)
	CFLAGS=-g -I/usr/local/pub/FFTW/3.2/include -DDEBUG -qnooptimize -qcheck=all -qheapdebug
	LDFLAGS=-g -qnooptimize -qcheck=all -qheapdebug -lm
endif