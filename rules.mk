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


all: data snoopy config.mk
	@(cd src && $(MAKE))
	cp src/snoopy .
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
	rm snoopy

fullclean:
	@(cd src && $(MAKE) $@)
	rm -rf snoopy config.mk timevar data dump.dmp