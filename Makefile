OPENMP=no
MPI=no
FFTW3_MPI=no
DEBUG=yes

MACHINE= $(shell uname -s)
HOSTNAME= $(shell hostname)

#######################################################
## Machine dependant variables
#######################################################
# Default variables:
CLUSTER="Unknown Cluster, using default"
CC=cc
FFTPATH=/usr/local/lib

#Compilation variables for MACOS X (geo laptop)
ifeq ($(MACHINE),Darwin)
	CLUSTER="MacOS 10.5 (Leopard)"
	ifeq ($(MPI),yes)
		CC=/usr/local/bin/mpicc
	else
		CC=/usr/local/bin/gcc
	endif
	FFTPATH=/usr/local/lib
	#FFTPATH=/Users/glesur/test/lib
	CFLAGS=-O3 -Wall
	OPENMP_FLAG=-fopenmp
endif

#Astro2 (DAMTP) on 32 bits
ifeq ($(HOSTNAME),astro2.damtp.cam.ac.uk)
	CLUSTER="SE Linux x86 (32 bits)"
	CC=icc						#No MPI yet...
	FFTPATH=/home/raid/chaos/gl293/usr/lib
	CFLAGS=-O3 -c99
	OPENMP_FLAG=-openmp
endif

# Hyades
ifeq ($(HOSTNAME),master.hyades.private.damtp.cam.ac.uk)
	CLUSTER="Hyades Cluster"
	CC=/opt/intel91039033/bin/icc		#No MPI yet...
	FFTPATH=/home/raid/chaos/gl293/usr/lib
	CFLAGS=-O3 -static
	OPENMP_FLAG=-openmp
endif

#test-vostro and similar configurations
ifeq ($(HOSTNAME),test-vostro.damtp.cam.ac.uk)
	CLUSTER="SE Linux x86 (64 bits)"
	CC=icc64						#No MPI yet...
	FFTPATH=/home/raid/chaos/gl293/vostro/usr/lib
	CFLAGS=-O3 -c99
	OPENMP_FLAG=-openmp
endif

###############################################################
## General compilation variables
###############################################################

LDFLAGS=-lfftw3

ifeq ($(DEBUG),yes)
	CFLAGS=-g -DDEBUG
	LDFLAGS+=-g
endif
ifeq ($(OPENMP),yes)
	CFLAGS+=$(OPENMP_FLAG) -DOPENMP_SUPPORT
	LDFLAGS+=-lfftw3_threads
endif
ifeq ($(FFTW3_MPI),yes)
	CFLAGS+=-DFFTW3_MPI_SUPPORT
	LDFLAGS+=-lfftw3_mpi
endif
ifeq ($(MPI),yes)
	CFLAGS+=-DMPI_SUPPORT
endif

export CC
export FFTPATH
export CFLAGS
export LDFLAGS

###############################################################
## Compilation rules
###############################################################


all: data
	@(cd src && $(MAKE))
	cp src/Cflow .
	@echo "***********************************************************"
	@echo "Make has compiled for: " $(CLUSTER)
ifeq ($(OPENMP),yes)
	@echo "OpenMP support enabled"
else
	@echo "OpenMP support disabled"
endif
ifeq ($(MPI),yes)
	@echo "MPI support enabled"
else
	@echo "MPI support disabled"
endif
	@echo "***********************************************************"

data:
	mkdir data

clean:
	@(cd src && $(MAKE) $@)
	rm Cflow



