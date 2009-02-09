# Uncomment the required compiler options

#Compilation variables for MACOS X (geo laptop)
export CC=/usr/local/bin/gcc
export FFTPATH=/usr/local/lib
export CFLAGS=-O3 -Wall -fopenmp
export FFTLIBRARIES=-lfftw3_threads -lfftw3

#Compilation variables for Astro2
#export CC=icc
#export FFTPATH=/home/raid/chaos/gl293/usr/lib
#export CFLAGS=-O3 -c99 -openmp
#export FFTLIBRARIES=-lfftw3_threads -lfftw3

#Compilation variables for hyades
#export CC=/opt/intel91039033/bin/icc
#export FFTPATH=/home/gl293/usr/lib
#export CFLAGS=-O3 -openmp -static
#export FFTLIBRARIES=-lfftw3_threads -lfftw3


all: data
	@(cd src && $(MAKE))
	cp src/Cflow .

data:
	mkdir data

clean:
	@(cd src && $(MAKE) $@)
	rm Cflow



