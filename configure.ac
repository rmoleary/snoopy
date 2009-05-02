AC_INIT([snoopy], [2.0])
AC_PROG_CC

# Save debug flags
AC_SUBST(CFLAGS_DEBUG)
AC_SUBST(LDFLAGS_DEBUG)
CFLAGS_DEBUG=$CFLAGS
LDFLAGS_DEBUG=$LDFLAGS

#---------------------------------------------------
# Check standard libraries

AC_CHECK_HEADERS([stdlib.h])
AC_CHECK_HEADERS([math.h])
AC_CHECK_HEADERS([complex.h],[ok=yes],[
				  echo "******************************************"
	              echo "Error! Complex type support not found."
				  echo "******************************************"
                  exit -1
                  ])

AC_CHECK_HEADERS([fftw3.h],[ok=yes],[
				  echo "******************************************"
	              echo "Error! FFTW library not found."
				  echo "******************************************"
                  exit -1
                  ])


AC_CHECK_LIB(fftw3,fftw_execute, [], )


#----------------------------------------------------
# MPI Feature (--enable-mpi)

AC_SUBST(MPI_ACTIVE)

MPI_ACTIVE=notavailable

AC_ARG_ENABLE(mpi, [  --enable-mpi     enable mpi parallelisation], ok=yes, ok=no)
if test "$ok" = "yes"; then
	ACX_MPI(MPI_ACTIVE=yes,[
				  echo "******************************************"
	              echo "Error! MPI libraries not found."
				  echo "******************************************"
                  exit -1
                  ])
fi

#----------------------------------------------------
#FFTW_MPI Feature (--enable-fftw-mpi)

AC_SUBST(FFTW_MPI_ACTIVE)

FFTW_MPI_ACTIVE=notavailable

AC_ARG_ENABLE(fftw-mpi, [  --enable-fftw-mpi     enable fftw3 experimental mpi support], ok=yes, ok=no)
if test "$ok" = "yes"; then
	AC_SEARCH_LIBS(fftw_mpi_init, fftw3_mpi fftw3 fftw3_threads,[
				  echo "******************************************"
	              echo "Error! FFTW MPI libraries not found."
				  echo "******************************************"
                  exit -1
                  ])
fi

#----------------------------------------------------
# OpenMP feature (--enable-openmp)

AC_SUBST(OPENMP_CFLAGS)
AC_SUBST(OPENMP_ACTIVE)

OPENMP_ACTIVE=notavailable

AC_ARG_ENABLE(openmp, [  --enable-openmp     enable OpenMP support], ok=yes, ok=no)
if test "$ok" = "yes"; then
	#check if openMP is supported by CC
	AX_OPENMP(ok=yes,[
				  echo "******************************************"
	              echo "Error! OpenMP not supported by compiler."
				  echo "******************************************"
                  exit -1
                  ])

	#check if Threads or OpenMP is supported by fftw
	
	AC_SEARCH_LIBS(fftw_init_threads, fftw3_threads fftw3, OPENMP_ACTIVE=yes,[
				  echo "************************************************"
	              echo "Error! Threads or OpenMP not supported by FFTW. "
				  echo "************************************************"
                  exit -1
                  ])
				  
fi

#----------------------------------------------------
# Debug options (--enable-debug)

AC_SUBST(DEBUG_ACTIVE)

DEBUG_ACTIVE=no

AC_ARG_ENABLE(debug, [  --enable-debug     enable Debug options], ok=yes, ok=no)
if test "$ok" = "yes"; then
	DEBUG_ACTIVE=yes
fi


#----------------------------------------------------
#Check best optimisation
AX_CC_MAXOPT 

#Create a gvars.h if none exist
if test ! -e src/gvars.h; then
	echo "creating default src/gvars.h"
	cp src/def/gvars_default.h src/gvars.h
fi

AC_CONFIG_FILES([Makefile])
AC_OUTPUT