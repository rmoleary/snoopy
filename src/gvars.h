//#define		DEBUG

// Caution: when using MPI, NX and NY must be multiples of NPROC
#define		NX				128
#define		NY				64
#define		NZ				64

#define		NTHREADS		2

#ifndef NPROC
#define		NPROC			16
#endif

#define		LX				2.0
#define		LY				1.0
#define		LZ				1.0

#define		CFL				1.5
#define		REYNOLDS		1000.0
#define		REYNOLDS_TH		1000.0

#define		BOUSSINESQ
#define		N2				(-1.0)

#define		OMEGA			(0.0)
//#define		WITH_SHEAR
#define		SHEAR			1.0

#define		PER_AMPLITUDE	100.0

#define		FORCING_TIME	4e-2

#define		T_INITIAL		0.0
#define		T_FINAL			7.0

#define		TOUTPUT_TIME	0.1
#define		TOUTPUT_FLOW	0.1
#define		TOUTPUT_DUMP	1.0

//#define		RESTART

//#define		FORTRAN_OUTPUT_ORDER
#define		VTK_OUTPUT

#define		INTERFACE_CHECK	4
//#define		INTERFACE_OUTPUT_FILE

#define ANTIALIASING

// FFT_PLANNING can be either FFTW_ESTIMATE, FFTW_MEASURE, FFTW_PATIENT or FFTW_EXHAUSTIVE (see fftw3 doc)
//#define		FFT_PLANNING	FFTW_ESTIMATE
#define		FFT_PLANNING	FFTW_MEASURE

// Ignore NPROC if MPI_SUPPORT is disabled
#ifndef	MPI_SUPPORT
#undef		NPROC
#define		NPROC			1
#endif

#define		NTOTAL			NX * NY * NZ

#define		NX_COMPLEX		NX
#define		NY_COMPLEX		NY
#define		NZ_COMPLEX		(NZ / 2 + 1)

#define		NTOTAL_COMPLEX	(NX_COMPLEX * NY_COMPLEX * NZ_COMPLEX / NPROC)

#define		PRECISION		double

#define		IDX3D			(k + j * NZ_COMPLEX + NZ_COMPLEX * NY_COMPLEX * i)

