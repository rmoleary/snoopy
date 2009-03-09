#include <math.h>
#include <complex.h>

#include "gvars.h"
#include "common.h"
#include "gfft.h"

void timestep( struct Field dfldo,
			   struct Field fldi,
			   const double t,
			   const double dt ) {
			   
	int i;
	PRECISION complex q0;
	// This is the timesteping algorithm, solving the physics.

/******************************************
** Velocity Self Advection ****************
*******************************************/

#pragma omp parallel private(i) num_threads ( NTHREADS )
{
		/* Compute the convolution */
	#pragma omp for schedule(static) nowait	
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i] =  fldi.vx[i];
		w2[i] =  fldi.vy[i];
		w3[i] =  fldi.vz[i];
	}
}

	gfft_c2r_t(w1);
	gfft_c2r_t(w2);
	gfft_c2r_t(w3);
	
#pragma omp parallel private(i) num_threads ( NTHREADS )
{
		/* Compute the convolution for the advection process */
	#pragma omp for schedule(static) nowait
	for( i = 0 ; i < 2*NTOTAL_COMPLEX ; i++) {
		wr4[i] = wr1[i] * wr1[i] / ((double) NTOTAL*NTOTAL);
		wr5[i] = wr2[i] * wr2[i] / ((double) NTOTAL*NTOTAL);
		wr6[i] = wr3[i] * wr3[i] / ((double) NTOTAL*NTOTAL);
		wr7[i] = wr1[i] * wr2[i] / ((double) NTOTAL*NTOTAL);
		wr8[i] = wr1[i] * wr3[i] / ((double) NTOTAL*NTOTAL);
		wr9[i] = wr2[i] * wr3[i] / ((double) NTOTAL*NTOTAL);
		
	}
}
	gfft_r2c_t(wr4);
	gfft_r2c_t(wr5);
	gfft_r2c_t(wr6);
	gfft_r2c_t(wr7);
	gfft_r2c_t(wr8);
	gfft_r2c_t(wr9);
		
#pragma omp parallel private(i) num_threads ( NTHREADS )
{
	#pragma omp for schedule(static ) nowait
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		dfldo.vx[i] = - I * mask[i] * (
					kxt[i] * w4[i] + ky[i] * w7[i] + kz[i] * w8[i] );
		dfldo.vy[i] = - I * mask[i] * (
					kxt[i] * w7[i] + ky[i] * w5[i] + kz[i] * w9[i] );
		dfldo.vz[i] = - I * mask[i] * (
					kxt[i] * w8[i] + ky[i] * w9[i] + kz[i] * w6[i] );
	}
}

/**********************************************
** BOUSSINESQ TERMS (if needed) ***************
***********************************************/

#ifdef BOUSSINESQ
#pragma omp parallel private(i) num_threads ( NTHREADS )
{
	#pragma omp for schedule(static) nowait	
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w4[i] = fldi.th[i];
	}
	
}

	gfft_c2r_t(w4);
		
#pragma omp parallel private(i) num_threads ( NTHREADS )
{
		/* Compute the convolution */
	#pragma omp for schedule(static) nowait	
	for( i = 0 ; i < 2*NTOTAL_COMPLEX ; i++) {		
		wr5[i] = wr1[i] * wr4[i] / ((double) NTOTAL*NTOTAL);
		wr6[i] = wr2[i] * wr4[i] / ((double) NTOTAL*NTOTAL);
		wr7[i] = wr3[i] * wr4[i] / ((double) NTOTAL*NTOTAL);
	}
}

	gfft_r2c_t(wr5);
	gfft_r2c_t(wr6);
	gfft_r2c_t(wr7);
	
#pragma omp parallel private(i) num_threads ( NTHREADS )
{
	#pragma omp for schedule(static ) nowait
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		dfldo.vx[i] -= N2 * fldi.th[i];
						
		dfldo.th[i] = - I * mask[i] * (
			kxt[i] * w5[i] + ky[i] * w6[i] + kz[i] * w7[i])
			+ fldi.vx[i];
	}
}
#endif

/************************************
** SOURCE TERMS  ********************
************************************/

#pragma omp parallel private(i) num_threads ( NTHREADS )
{
	#pragma omp for schedule(static )
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		dfldo.vx[i] += 2.0 * OMEGA * fldi.vy[i];
#ifdef WITH_SHEAR
		dfldo.vy[i] += (SHEAR - 2.0 * OMEGA) * fldi.vx[i];
#else
		dfldo.vy[i] += (- 2.0 * OMEGA) * fldi.vx[i];
#endif
	}
			
			
/************************************
** PRESSURE TERMS *******************
************************************/
	#pragma omp for schedule(static ) nowait
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
#ifdef WITH_SHEAR
		q0= SHEAR * ky[i] * fldi.vx[i] + kxt[i] * dfldo.vx[i] + ky[i] * dfldo.vy[i] + kz[i] * dfldo.vz[i];
#else
		q0= kxt[i] * dfldo.vx[i] + ky[i] * dfldo.vy[i] + kz[i] * dfldo.vz[i];
#endif
		dfldo.vx[i] += -kxt[i]* q0 * ik2t[i];
		dfldo.vy[i] += -ky[i] * q0 * ik2t[i];
		dfldo.vz[i] += -kz[i] * q0 * ik2t[i];
	}
}
		
	return;
}
		
/************************************
** Implicit steps called by mainloop
*************************************/
 
void implicitstep(
			   struct Field fldi,
			   const double t,
			   const double dt ) {
			   
	PRECISION q0;
	int i;
#pragma omp parallel private(i,q0) num_threads ( NTHREADS )
{
		/* Compute the convolution */
	#pragma omp for schedule(static ) nowait
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		q0 = exp( - nu * dt* k2t[i] );
		fldi.vx[i] = fldi.vx[i] * q0;
		fldi.vy[i] = fldi.vy[i] * q0;
		fldi.vz[i] = fldi.vz[i] * q0;
		
#ifdef BOUSSINESQ
		q0 = exp( - nu_th * dt* k2t[i] );
		fldi.th[i] = fldi.th[i] * q0;
#endif
	}
}
	
	return;
}
	
			   
			   