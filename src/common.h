/*	
	This file is part of the Snoopy code.

    Snoopy code is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Snoopy code is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Snoopy code.  If not, see <http://www.gnu.org/licenses/>.
*/



#include "snoopy.h"

#define		CHECK_NAN(XIN)		c_nan(XIN, __func__, __LINE__,__FILE__)
 
// All these variables may be used in the code as they are initialized by common.c
// Wave number pointers
extern double	*kx,	*ky,	*kz,	*kxt,	*k2t,	*ik2t;
extern double	kxmax,	kymax,  kzmax,	kmax;


// Mask for dealiasing
extern double   *mask;

extern double	*wr1,	*wr2,	*wr3;
extern double	*wr4,	*wr5,	*wr6;
extern double	*wr7,	*wr8,	*wr9;
extern double	*wr10,	*wr11,	*wr12;
extern double	*wr13,	*wr14,	*wr15;

extern double complex		*w1,	*w2,	*w3;
extern double complex		*w4,	*w5,	*w6;
extern double complex		*w7,	*w8,	*w9;
extern double complex		*w10,	*w11,	*w12;
extern double complex		*w13,	*w14,	*w15;


// Parameters
extern struct Parameters			param;

// Physics variables 
extern double	nu;

#ifdef BOUSSINESQ
extern double	nu_th;
#ifdef N2PROFILE
extern double *N2_profile;
#endif
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
void allocate_field(struct Field *fldi);
void deallocate_field(struct Field *fldi);
double get_c_time(void);

// Useful only if MPI is active. Can be called without though...
void reduce(double *var, const int op);

float big_endian(float in_number);

double randm_normal(void);
double randm (void);

void projector( double complex qx[],
			    double complex qy[],
			    double complex qz[]);
				
#ifdef ANELASTIC
void projector_anelastic( double complex qx[],
			    double complex qy[],
			    double complex qz[]);
#endif
				
double energy(const double complex q[]);


				