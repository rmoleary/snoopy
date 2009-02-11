#define		DEBUG

#define		NX				64
#define		NY				64
#define		NZ				64

#define		NTHREADS		2

#define		LX				1.0
#define		LY				1.0
#define		LZ				1.0

#define		CFL				1.5
#define		REYNOLDS		4000
#define		REYNOLDS_TH		400

#define		BOUSSINESQ
#define		N2				(-0.01)

#define		OMEGA			(2.0 / 3.0)
#define		WITH_SHEAR
#define		SHEAR			1.0

#define		PER_AMPLITUDE	1.0

#define		FORCING_TIME	4e-2

#define		T_INITIAL		0.0
#define		T_FINAL			10000.0

#define		TOUTPUT_TIME	0.01
#define		TOUTPUT_FLOW	0.01
#define		TOUTPUT_DUMP	1.0

//#define		RESTART

#define		FORTRAN_OUTPUT_ORDER

#define		INTERFACE_CHECK	4
//#define		INTERFACE_OUTPUT_FILE

// ANTIALIASING option need to be coded in init_code routine
#define		ANTIALIASING

// USE the INPLACE option to save some of the memory used by temporary arrays.
//#define		INPLACE

// FFT_PLANNING can be either FFTW_ESTIMATE, FFTW_MEASURE, FFTW_PATIENT or FFTW_EXHAUSTIVE (see fftw3 doc)
//#define		FFT_PLANNING	FFTW_ESTIMATE
#define		FFT_PLANNING	FFTW_MEASURE

#define		NTOTAL			NX * NY * NZ

#define		NX_COMPLEX		NX
#define		NY_COMPLEX		NY
#define		NZ_COMPLEX		(NZ / 2 + 1)

#define		NTOTAL_COMPLEX			NX_COMPLEX * NY_COMPLEX * NZ_COMPLEX

#define		PRECISION		double

#define		IDX3D			(k + j * NZ_COMPLEX + NZ_COMPLEX * NY_COMPLEX * i)

