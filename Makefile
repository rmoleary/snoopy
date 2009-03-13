OPENMP=yes
MPI=yes
FFTW3_MPI=no
DEBUG=yes

MACHINE= $(shell uname -s)
HOSTNAME= $(shell hostname)
DOMAINNAME= $(shell dnsdomainname)

#######################################################
## Machine dependant variables
#######################################################
# Default variables:
CLUSTER="Unknown Cluster, using default"
CC=cc
FFTPATH=-L/usr/local/lib
LDFLAGS=
CFLAGS=

#Compilation variables for MACOS X (geo laptop)
ifeq ($(MACHINE),Darwin)
	CLUSTER="MacOS 10.5 (Leopard)"
	ifeq ($(MPI),yes)
		CC=/usr/local/bin/mpicc
	else
		CC=/usr/local/bin/gcc
	endif
	FFTPATH=-L/usr/local/lib
	#FFTPATH=/Users/glesur/test/lib
	CFLAGS=-O3 -Wall
	OPENMP_FLAG=-fopenmp
ifeq ($(DEBUG),yes)
	CFLAGS=-g -DDEBUG
	LDFLAGS+=-g
endif
endif

#Astro2 (DAMTP) on 32 bits
ifeq ($(HOSTNAME),astro2.damtp.cam.ac.uk)
	CLUSTER="SE Linux x86 (32 bits)"
	CC=icc						#No MPI yet...
	FFTPATH=-L/home/raid/chaos/gl293/usr/lib
	CFLAGS=-O3 -c99
	OPENMP_FLAG=-openmp
ifeq ($(DEBUG),yes)
	CFLAGS=-g -DDEBUG
	LDFLAGS+=-g
endif
endif

# Hyades
ifeq ($(HOSTNAME),master.hyades.private.damtp.cam.ac.uk)
	CLUSTER="Hyades Cluster"
	CC=/opt/intel91039033/bin/icc		#No MPI yet...
	FFTPATH=-L/home/raid/chaos/gl293/usr/lib
	CFLAGS=-O3 -static
	OPENMP_FLAG=-openmp
ifeq ($(DEBUG),yes)
	CFLAGS=-g -DDEBUG
	LDFLAGS+=-g
endif
endif

#test-vostro and similar configurations
ifeq ($(HOSTNAME),test-vostro.damtp.cam.ac.uk)
	CLUSTER="SE Linux x86 (64 bits)"
	CC=icc64						#No MPI yet...
	FFTPATH=-L/home/raid/chaos/gl293/vostro/usr/lib
	CFLAGS=-O3 -c99
	OPENMP_FLAG=-openmp
ifeq ($(DEBUG),yes)
	CFLAGS=-g -DDEBUG
	LDFLAGS+=-g
endif
endif

ifeq ($(DOMAINNAME),moab.cluster)
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
endif

ifeq ($(MACHINE),AIX)
	CLUSTER="VARGAS (IBM AIX)"
	ifeq ($(MPI),yes)
		CC=mpcc_r
	else
		CC=xlc_r
	endif
	CFLAGS=-O3
	FFTPATH=
	OPENMP_FLAG= -qsmp=omp
ifeq ($(DEBUG),yes)
	CFLAGS=-g -DDEBUG -qnooptimize -qcheck=all -qheapdebug
	LDFLAGS=-g -qnooptimize -qcheck=all -qheapdebug
endif
endif


###############################################################
## General compilation variables
###############################################################




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

LDFLAGS+=-lfftw3
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



