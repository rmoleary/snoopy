#include <math.h>
#include <complex.h>

#include "common.h"
#include "gfft.h"
#include "debug.h"

#ifdef FORCING
/*************************************************
** Forcing (if required) *************************
**************************************************/
void forcing(struct Field fldi,
			 double dt) {
			 
// Force random velocity field
	const double kf = 2.0 * M_PI * 5.0;
	const double deltakf = kf * 0.2;
	const double amplitude_forcing = 5.0;
	int found_k = 0;
	int i;
// Find a valid k vector
	while(found_k==0) {
		i = (int) (randm() * NTOTAL_COMPLEX);
		if( (k2t[i]>(kf-deltakf)*(kf-deltakf)) && (k2t[i]<(kf+deltakf)*(kf+deltakf))) {
			found_k=1;
		}
	}
	
// Init V vector
	fldi.vx[i] += amplitude_forcing * randm_normal() * NTOTAL * cexp( I * 2.0*M_PI*randm() ) *dt;
	fldi.vy[i] += amplitude_forcing * randm_normal() * NTOTAL * cexp( I * 2.0*M_PI*randm() ) *dt;
	fldi.vz[i] += amplitude_forcing * randm_normal() * NTOTAL * cexp( I * 2.0*M_PI*randm() ) *dt;
	
	projector(fldi.vx,fldi.vy,fldi.vz);
	
	return;
}
#endif


/**************************************************
*** Timestep, called by runge-kutta loop   ********
***************************************************/

void timestep( struct Field dfldo,
			   struct Field fldi,
			   const double t,
			   const double dt ) {
			   
	int i;
	double complex q0;
	double S;
	// This is the timesteping algorithm, solving the physics.

	// Find the shear at time t
#ifdef WITH_SHEAR
#ifdef TIME_DEPENDANT_SHEAR
	S = param.shear * param.omega_shear * cos(param.omega_shear * t);
#else
	S = param.shear;
#endif
#endif

/******************************************
** Velocity Self Advection ****************
*******************************************/

		/* Compute the convolution */
#ifdef _OPENMP
	#pragma omp parallel for private(i) schedule(static)	
#endif
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i] =  fldi.vx[i];
		w2[i] =  fldi.vy[i];
		w3[i] =  fldi.vz[i];
	}

	gfft_c2r_t(w1);
	gfft_c2r_t(w2);
	gfft_c2r_t(w3);
	
		/* Compute the convolution for the advection process */
#ifdef _OPENMP
	#pragma omp parallel for private(i) schedule(static)	
#endif
	for( i = 0 ; i < 2*NTOTAL_COMPLEX ; i++) {
		wr4[i] = wr1[i] * wr1[i] / ((double) NTOTAL*NTOTAL);
		wr5[i] = wr2[i] * wr2[i] / ((double) NTOTAL*NTOTAL);
		wr6[i] = wr3[i] * wr3[i] / ((double) NTOTAL*NTOTAL);
		wr7[i] = wr1[i] * wr2[i] / ((double) NTOTAL*NTOTAL);
		wr8[i] = wr1[i] * wr3[i] / ((double) NTOTAL*NTOTAL);
		wr9[i] = wr2[i] * wr3[i] / ((double) NTOTAL*NTOTAL);
	}
	
	gfft_r2c_t(wr4);
	gfft_r2c_t(wr5);
	gfft_r2c_t(wr6);
	gfft_r2c_t(wr7);
	gfft_r2c_t(wr8);
	gfft_r2c_t(wr9);

#ifdef _OPENMP
	#pragma omp parallel for private(i) schedule(static)	
#endif
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		dfldo.vx[i] = - I * mask[i] * (
					kxt[i] * w4[i] + ky[i] * w7[i] + kz[i] * w8[i] );
		dfldo.vy[i] = - I * mask[i] * (
					kxt[i] * w7[i] + ky[i] * w5[i] + kz[i] * w9[i] );
		dfldo.vz[i] = - I * mask[i] * (
					kxt[i] * w8[i] + ky[i] * w9[i] + kz[i] * w6[i] );
		
	}


/**********************************************
** BOUSSINESQ TERMS (if needed) ***************
***********************************************/

#ifdef BOUSSINESQ
#ifdef _OPENMP
	#pragma omp parallel for private(i) schedule(static)	
#endif	
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w4[i] = fldi.th[i];
	}
	
	gfft_c2r_t(w4);
		
#ifdef _OPENMP
	#pragma omp parallel for private(i) schedule(static)	
#endif	
	for( i = 0 ; i < 2*NTOTAL_COMPLEX ; i++) {		
		wr5[i] = wr1[i] * wr4[i] / ((double) NTOTAL*NTOTAL);
		wr6[i] = wr2[i] * wr4[i] / ((double) NTOTAL*NTOTAL);
		wr7[i] = wr3[i] * wr4[i] / ((double) NTOTAL*NTOTAL);
	}

	gfft_r2c_t(wr5);
	gfft_r2c_t(wr6);
	gfft_r2c_t(wr7);
	
#ifdef _OPENMP
	#pragma omp parallel for private(i) schedule(static)	
#endif
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
#ifdef VERTSTRAT
		dfldo.vz[i] -= param.N2 * fldi.th[i];
						
		dfldo.th[i] = - I * mask[i] * (
			kxt[i] * w5[i] + ky[i] * w6[i] + kz[i] * w7[i])
			+ fldi.vz[i];
#else
		dfldo.vx[i] -= param.N2 * fldi.th[i];
						
		dfldo.th[i] = - I * mask[i] * (
			kxt[i] * w5[i] + ky[i] * w6[i] + kz[i] * w7[i])
			+ fldi.vx[i];
#endif
	}
	
	
#endif

/*********************************************
**** MHD Terms (if needed)   *****************
*********************************************/
#ifdef MHD

// Start with the induction equation
#ifdef _OPENMP
	#pragma omp parallel for private(i) schedule(static)	
#endif
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w4[i] =  fldi.bx[i];
		w5[i] =  fldi.by[i];
		w6[i] =  fldi.bz[i];
	}


	gfft_c2r_t(w4);
	gfft_c2r_t(w5);
	gfft_c2r_t(w6);
	
	// (vx,vy,vz) is in w1-w3 and (bx,by,bz) is in (w4-w6). It is now time to compute the emfs in w7-w9...
#ifdef _OPENMP
	#pragma omp parallel for private(i) schedule(static)	
#endif
	for( i = 0 ; i < 2*NTOTAL_COMPLEX ; i++) {
		wr7[i] = (wr2[i] * wr6[i] - wr3[i] * wr5[i]) / ((double) NTOTAL*NTOTAL);
		wr8[i] = (wr3[i] * wr4[i] - wr1[i] * wr6[i]) / ((double) NTOTAL*NTOTAL);
		wr9[i] = (wr1[i] * wr5[i] - wr2[i] * wr4[i]) / ((double) NTOTAL*NTOTAL);
	}

	// Compute the curl of the emf to add in the induction equation.
	
	gfft_r2c_t(wr7);
	gfft_r2c_t(wr8);
	gfft_r2c_t(wr9);
	
#ifdef _OPENMP
	#pragma omp parallel for private(i) schedule(static)	
#endif
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		dfldo.bx[i] = I * mask[i] * (ky[i] * w9[i] - kz[i] * w8[i]);
		dfldo.by[i] = I * mask[i] * (kz[i] * w7[i] - kxt[i]* w9[i]);
		dfldo.bz[i] = I * mask[i] * (kxt[i]* w8[i] - ky[i] * w7[i]);
	}


// Let's do the Lorentz Force
// We already have (bx,by,bz) in w4-w6. No need to compute them again...

#ifdef _OPENMP
	#pragma omp parallel for private(i) schedule(static)	
#endif
	for( i = 0 ; i < 2*NTOTAL_COMPLEX ; i++) {
		wr1[i] = wr4[i] * wr4[i] / ((double) NTOTAL*NTOTAL);
		wr2[i] = wr5[i] * wr5[i] / ((double) NTOTAL*NTOTAL);
		wr3[i] = wr6[i] * wr6[i] / ((double) NTOTAL*NTOTAL);
		wr7[i] = wr4[i] * wr5[i] / ((double) NTOTAL*NTOTAL);
		wr8[i] = wr4[i] * wr6[i] / ((double) NTOTAL*NTOTAL);
		wr9[i] = wr5[i] * wr6[i] / ((double) NTOTAL*NTOTAL);
	}


	gfft_r2c_t(wr1);
	gfft_r2c_t(wr2);
	gfft_r2c_t(wr3);
	gfft_r2c_t(wr7);
	gfft_r2c_t(wr8);
	gfft_r2c_t(wr9);

#ifdef _OPENMP
	#pragma omp parallel for private(i) schedule(static)
#endif
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		dfldo.vx[i] += I * mask[i] * (kxt[i] * w1[i] + ky[i] * w7[i] + kz[i] * w8[i]);
		dfldo.vy[i] += I * mask[i] * (kxt[i] * w7[i] + ky[i] * w2[i] + kz[i] * w9[i]);
		dfldo.vz[i] += I * mask[i] * (kxt[i] * w8[i] + ky[i] * w9[i] + kz[i] * w3[i]);
	}
	
#endif


/************************************
** SOURCE TERMS  ********************
************************************/

#ifdef _OPENMP
	#pragma omp parallel for private(i) schedule(static)
#endif
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
#ifdef WITH_ROTATION
		dfldo.vx[i] += 2.0 * param.omega * fldi.vy[i];
		dfldo.vy[i] -= 2.0 * param.omega * fldi.vx[i];
#endif
#ifdef WITH_SHEAR
		dfldo.vy[i] += S  * fldi.vx[i];
#ifdef MHD
		dfldo.by[i] -= S * fldi.bx[i];
#endif		
#endif
	}
			
/************************************
** PRESSURE TERMS *******************
************************************/
#ifdef _OPENMP
	#pragma omp parallel for private(i,q0) schedule(static)
#endif
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
#ifdef WITH_SHEAR
		q0= S * ky[i] * fldi.vx[i] + kxt[i] * dfldo.vx[i] + ky[i] * dfldo.vy[i] + kz[i] * dfldo.vz[i];
#else
		q0= kxt[i] * dfldo.vx[i] + ky[i] * dfldo.vy[i] + kz[i] * dfldo.vz[i];
#endif
		dfldo.vx[i] += -kxt[i]* q0 * ik2t[i];
		dfldo.vy[i] += -ky[i] * q0 * ik2t[i];
		dfldo.vz[i] += -kz[i] * q0 * ik2t[i];
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
			   
	double q0;
	int i;
#ifdef _OPENMP
	#pragma omp parallel for private(i,q0) schedule(static)
#endif
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		q0 = exp( - nu * dt* k2t[i] );
		fldi.vx[i] = fldi.vx[i] * q0;
		fldi.vy[i] = fldi.vy[i] * q0;
		fldi.vz[i] = fldi.vz[i] * q0;
		
#ifdef BOUSSINESQ
		q0 = exp( - nu_th * dt* k2t[i] );
		fldi.th[i] = fldi.th[i] * q0;
#endif
#ifdef MHD
		q0 = exp( - eta * dt* k2t[i] );
		fldi.bx[i] = fldi.bx[i] * q0;
		fldi.by[i] = fldi.by[i] * q0;
		fldi.bz[i] = fldi.bz[i] * q0;
#endif

	}

#ifdef FORCING
	forcing(fldi, dt);
#endif

	return;
}
	
	


			   