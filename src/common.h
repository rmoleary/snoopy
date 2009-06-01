#include "snoopy.h"

// All these variables may be used in the code as they are initialized by common.c
// Wave number pointers
extern double	*kx,	*ky,	*kz,	*kxt,	*k2t,	*ik2t;
extern double	kxmax,	kymax,  kzmax,	kmax;


// Mask for dealiasing
extern double   *mask;

extern double	*wr1,	*wr2,	*wr3;
extern double	*wr4,	*wr5,	*wr6;
extern double	*wr7,	*wr8,	*wr9;
extern double	*wr10;

extern struct Field				fld;

extern double complex		*w1,	*w2,	*w3;
extern double complex		*w4,	*w5,	*w6;
extern double complex		*w7,	*w8,	*w9;

extern double complex		*pressure;

// Parameters
extern struct Parameters			param;

// Physics variables 
extern double	nu;

#ifdef BOUSSINESQ
extern double	nu_th;
#endif

#ifdef MHD
extern double	eta;
#endif

// MPI
#ifdef MPI_SUPPORT
extern int		NPROC;									/**< NPROC is a variable when MPI is on. Otherwise, it is preprocessor macro in gvars.h */
#endif
extern int rank;

// OpenMP
extern int	nthreads;

// Functions provided by the common routine

void init_common ( void );
void finish_common ( void );
double get_c_time(void);

// Useful only if MPI is active. Can be called without though...
void reduce(double *var, const int op);

double randm_normal(void);
double randm (void);

void projector( double complex qx[],
			    double complex qy[],
			    double complex qz[]);
				
double energy(const double complex q[]);

void enforce_symm(struct Field fldi);
void symmetrize(double complex wi[]);
				
				