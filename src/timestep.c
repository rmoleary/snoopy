#include <math.h>
#include <complex.h>
#include <fftw3.h>

#include "gvars.h"
#include "common.h"

void timestep( struct Field dfldo,
			   struct Field fldi,
			   const double t,
			   const double dt ) {
			   
	int i;
	// This is the timesteping algorithm, solving the physics.

#pragma omp parallel private(i) num_threads ( NTHREADS )
{
		/* Compute the convolution */
	#pragma omp for schedule(static) nowait	
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i] =  fldi.vx[i];
		w2[i] =  fldi.vy[i];
		w3[i] =  fldi.vz[i];
#ifdef BOUSSINESQ
		w6[i] = fldi.th[i];
#endif
	}
}
	fftw_execute_dft_c2r( c2rfft, w1, wr1);
	fftw_execute_dft_c2r( c2rfft, w2, wr2);
	fftw_execute_dft_c2r( c2rfft, w3, wr3);
	
	
#ifdef BOUSSINESQ
	fftw_execute_dft_c2r( c2rfft, w6, wr6);
#endif
	
#pragma omp parallel private(i) num_threads ( NTHREADS )
{
		/* Compute the convolution for the advection process */
	#pragma omp for schedule(static) nowait
	for( i = 0 ; i < NTOTAL ; i++) {
		wr4[i] = wr1[i] * wr1[i] / ((double) NTOTAL*NTOTAL);
		wr5[i] = wr2[i] * wr2[i] / ((double) NTOTAL*NTOTAL);
		wr6[i] = wr3[i] * wr3[i] / ((double) NTOTAL*NTOTAL);
		wr7[i] = wr1[i] * wr2[i] / ((double) NTOTAL*NTOTAL);
		wr8[i] = wr1[i] * wr3[i] / ((double) NTOTAL*NTOTAL);
		wr9[i] = wr2[i] * wr3[i] / ((double) NTOTAL*NTOTAL);
		
	}
}
	fftw_execute_dft_r2c( r2cfft, wr4, w4);
	fftw_execute_dft_r2c( r2cfft, wr5, w5);	
	fftw_execute_dft_r2c( r2cfft, wr6, w6);
	fftw_execute_dft_r2c( r2cfft, wr7, w7);
	fftw_execute_dft_r2c( r2cfft, wr8, w8);
	fftw_execute_dft_r2c( r2cfft, wr9, w9);	
	
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

#ifdef BOUSSINESQ
#pragma omp parallel private(i) num_threads ( NTHREADS )
{
	#pragma omp for schedule(static) nowait	
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w4[i] = fldi.th[i];
	}
	
}
	fftw_execute_dft_c2r( c2rfft, w4, wr4);
		
#pragma omp parallel private(i) num_threads ( NTHREADS )
{
		/* Compute the convolution */
	#pragma omp for schedule(static) nowait	
	for( i = 0 ; i < NTOTAL ; i++) {		
		wr5[i] = wr1[i] * wr4[i] / ((double) NTOTAL*NTOTAL);
		wr6[i] = wr2[i] * wr4[i] / ((double) NTOTAL*NTOTAL);
		wr7[i] = wr3[i] * wr4[i] / ((double) NTOTAL*NTOTAL);
	}
}

	fftw_execute_dft_r2c( r2cfft, wr5, w5);
	fftw_execute_dft_r2c( r2cfft, wr6, w6);
	fftw_execute_dft_r2c( r2cfft, wr7, w7);
	
#pragma omp parallel private(i) num_threads ( NTHREADS )
{
	#pragma omp for schedule(static ) nowait
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		dfldo.vx[i] += N2 * fldi.th[i];
						
		dfldo.th[i] = - I * mask[i] * (
			kxt[i] * w5[i] + ky[i] * w6[i] + kz[i] * w7[i])
			- fldi.vx[i];
	}
}
#endif
		
	return;
}
		
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
		fldi.wz[i] = fldi.wz[i] * q0;
		
#ifdef BOUSSINESQ
		q0 = exp( - nu_th * dt* k2t[i] );
		fldi.th[i] = fldi.th[i] * q0;
#endif
	}
}
	
	return;
}
	
			   
			   