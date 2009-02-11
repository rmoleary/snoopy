#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include "gvars.h"
#include "error.h"


// Structures

struct Field {
	PRECISION complex *vx;
	PRECISION complex *vy;
	PRECISION complex *vz;
#ifdef BOUSSINESQ
	PRECISION complex *th;
#endif
};


// All these variables may be used in the code as they are initialized by common.c
// Wave number pointers
extern PRECISION	*kx,	*ky,	*kz,	*kxt,	*k2t,	*ik2t;
extern PRECISION	kxmax,	kymax,  kzmax,	kmax;

extern fftw_plan	r2cfft,	c2rfft, fft_1d_forward, fft_1d_backward;


// Mask for dealiasing
extern PRECISION   *mask;

extern PRECISION	*wr1,	*wr2,	*wr3;
extern PRECISION	*wr4,	*wr5,	*wr6;
extern PRECISION	*wr7,	*wr8,	*wr9;

extern struct Field				fld;

extern PRECISION complex		*w1,	*w2,	*w3;
extern PRECISION complex		*w4,	*w5,	*w6;
extern PRECISION complex		*w7,	*w8,	*w9;
extern PRECISION complex		*w10;

extern PRECISION complex		*w1d, *w2d;

// Physics variables 
extern PRECISION	nu;

#ifdef BOUSSINESQ
extern PRECISION	nu_th;
#endif

// Functions provided by the common routine

void init_common ( void );
void finish_common ( void );

PRECISION randm (void);

void projector( PRECISION complex qx[],
			    PRECISION complex qy[],
			    PRECISION complex qz[]);
				
PRECISION energy(const PRECISION complex q[]);
				
				