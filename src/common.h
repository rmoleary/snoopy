#include <math.h>
#include <complex.h>
#include <fftw3.h>

#ifdef MPI_SUPPORT
#include <mpi.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#include "gvars.h"
#include "error.h"


#ifdef MPI_SUPPORT
#define		MPI_Printf			if (rank==0) printf
#else
#define		MPI_Printf			printf
#endif

// Structures

struct Field {
	PRECISION complex *vx;
	PRECISION complex *vy;
	PRECISION complex *vz;
#ifdef BOUSSINESQ
	PRECISION complex *th;
#endif
#ifdef MHD
	PRECISION complex *bx;
	PRECISION complex *by;
	PRECISION complex *bz;
#endif
};

// All these variables may be used in the code as they are initialized by common.c
// Wave number pointers
extern PRECISION	*kx,	*ky,	*kz,	*kxt,	*k2t,	*ik2t;
extern PRECISION	kxmax,	kymax,  kzmax,	kmax;


// Mask for dealiasing
extern PRECISION   *mask;

extern PRECISION	*wr1,	*wr2,	*wr3;
extern PRECISION	*wr4,	*wr5,	*wr6;
extern PRECISION	*wr7,	*wr8,	*wr9;
extern PRECISION	*wr10;

extern struct Field				fld;

extern PRECISION complex		*w1,	*w2,	*w3;
extern PRECISION complex		*w4,	*w5,	*w6;
extern PRECISION complex		*w7,	*w8,	*w9;

// Physics variables 
extern PRECISION	nu;

#ifdef BOUSSINESQ
extern PRECISION	nu_th;
#endif

#ifdef MHD
extern PRECISION	eta;
#endif

// MPI
extern int rank;

// OpenMP
extern int	nthreads;

// Functions provided by the common routine

void init_common ( void );
void finish_common ( void );
double get_c_time(void);

// Useful only if MPI is active. Can be called without though...
void reduce(double *var, const int op);

PRECISION randm_normal(void);
PRECISION randm (void);

void projector( PRECISION complex qx[],
			    PRECISION complex qy[],
			    PRECISION complex qz[]);
				
PRECISION energy(const PRECISION complex q[]);

void symmetrize(PRECISION complex wi[]);
				
				