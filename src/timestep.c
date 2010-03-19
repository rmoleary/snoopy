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


#include <math.h>
#include <complex.h>

#include "common.h"
#include "gfft.h"
#include "debug.h"
#include "forcing.h"


/**************************************************
*** Timestep, called by runge-kutta loop   ********
***************************************************/

void timestep( struct Field dfldo,
			   struct Field fldi,
			   double complex *po,
			   const double t,
			   const double dt) {
			   
	int i;
	double complex q0,q1;
	double S;
	// This is the timesteping algorithm, solving the physics.

	// Find the shear at time t
#ifdef WITH_SHEAR
#ifdef TIME_DEPENDANT_SHEAR
	S = param.shear * cos(param.omega_shear * t);	// This is the real shear: -dvy/dx
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
#ifndef WITH_2D
		wr6[i] = wr3[i] * wr3[i] / ((double) NTOTAL*NTOTAL);
#endif
		wr7[i] = wr1[i] * wr2[i] / ((double) NTOTAL*NTOTAL);
		wr8[i] = wr1[i] * wr3[i] / ((double) NTOTAL*NTOTAL);
		wr9[i] = wr2[i] * wr3[i] / ((double) NTOTAL*NTOTAL);
	}
	
	gfft_r2c_t(wr4);
	gfft_r2c_t(wr5);
#ifndef WITH_2D
	gfft_r2c_t(wr6);
#endif
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
					kxt[i] * w8[i] + ky[i] * w9[i] + kz[i] * w6[i] );	// since kz=0 in 2D, kz*w6 gives 0, even if w6 is some random array
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
#ifndef WITH_2D
		wr7[i] = wr3[i] * wr4[i] / ((double) NTOTAL*NTOTAL);
#endif
#ifdef N2PROFILE
		wr8[i] = N2_profile[i] * wr4[i] / ((double) NTOTAL);
#endif
	}

	gfft_r2c_t(wr5);
	gfft_r2c_t(wr6);
#ifndef WITH_2D
	gfft_r2c_t(wr7);
#endif
#ifdef N2PROFILE
	gfft_r2c_t(wr8);
#endif
	
#ifdef _OPENMP
	#pragma omp parallel for private(i) schedule(static)	
#endif
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
#ifdef VERTSTRAT
#ifdef N2PROFILE
		dfldo.vz[i] -= w8[i] * mask[i];
#else
		dfldo.vz[i] -= param.N2 * fldi.th[i];
#endif
						
		dfldo.th[i] = - I * mask[i] * (
			kxt[i] * w5[i] + ky[i] * w6[i] + kz[i] * w7[i])
			+ fldi.vz[i];
#else
#ifdef N2PROFILE
		dfldo.vx[i] -= w8[i] * mask[i];
#else
		dfldo.vx[i] -= param.N2 * fldi.th[i];
#endif
						
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
	#pragma omp parallel for private(i,q0,q1) schedule(static)
#endif
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
	
#ifdef ANELASTIC
		// We have to use a different prescription here since we don't have div v = 0 but div rho v = 0
#ifdef WITH_SHEAR
		q0 = S * ky[i] * fldi.vx[i] + kxt[i] * dfldo.vx[i] + ky[i] * dfldo.vy[i] + kz[i] * dfldo.vz[i] - I * dfldo.vx[i] / param.anelastic_lambda;
#else
		q0 = kxt[i] * dfldo.vx[i] + ky[i] * dfldo.vy[i] + kz[i] * dfldo.vz[i] - I * dfldo.vx[i] / param.anelastic_lambda;
#endif
		
		// inverse Poisson equation in anelastic
		q1 = 1.0 / (param.anelastic_lambda * param.anelastic_lambda) + 2.0 * I * kxt[i] / param.anelastic_lambda - k2t[i];
		
		if(q1 != 0)
			q0 = q0 / q1;
		else
			q0=0.0;		// That means the Poisson equation has a singularity (only expected if lambda=infinity)
			
		dfldo.vx[i] += (kxt[i] - I / param.anelastic_lambda) * q0;
		dfldo.vy[i] += ky[i] * q0;
		dfldo.vz[i] += kz[i] * q0;
	
#else
		
#ifdef WITH_SHEAR
		q0= S * ky[i] * fldi.vx[i] + kxt[i] * dfldo.vx[i] + ky[i] * dfldo.vy[i] + kz[i] * dfldo.vz[i];
#else
		q0= kxt[i] * dfldo.vx[i] + ky[i] * dfldo.vy[i] + kz[i] * dfldo.vz[i];
#endif
		if(po != NULL) {
			po[i] = - I * ik2t[i] * q0;	// Save the pressure field (if needed)
		}
		dfldo.vx[i] += -kxt[i]* q0 * ik2t[i];
		dfldo.vy[i] += -ky[i] * q0 * ik2t[i];
		dfldo.vz[i] += -kz[i] * q0 * ik2t[i];
#endif
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
			   
	double q0, q1;
	int i,j,k;
	
#ifdef SGS		// Subgrid model
		// This is the Chollet-Lesieur Model (1981)
		// We have nu(k)=nu_i(k)*(E(kc)/kc)^(1/2)
		// nu_i(k)=0.267+9.21*exp(-3.03 kc/k)
		
		// Compute E(kc)
	double kc, dk;
	
	kc = 2.0 * M_PI * 50;
	dk = 2.5;
	
	q0 = 0.0;
	
	for( i = 0; i < NX_COMPLEX/NPROC; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			for( k = 0; k < NZ_COMPLEX; k++) {
				if( (k2t[ IDX3D ] < (kc+dk) * (kc+dk)) & (k2t[ IDX3D ] > (kc-dk) * (kc-dk) )) {
#ifdef WITH_2D
					if( j == 0)
#else
					if( k == 0)
#endif
						q0 = q0 + creal( fldi.vx[ IDX3D ] * conj( fldi.vx[ IDX3D ] ) +
										 fldi.vy[ IDX3D ] * conj( fldi.vy[ IDX3D ] ) +
										 fldi.vz[ IDX3D ] * conj( fldi.vz[ IDX3D ] ) ) / ((double) NTOTAL*NTOTAL);
					else 
						// k>0, only half of the complex plane is represented.
						q0 = q0 + 2.0 * creal( fldi.vx[ IDX3D ] * conj( fldi.vx[ IDX3D ] ) +
											   fldi.vy[ IDX3D ] * conj( fldi.vy[ IDX3D ] ) +
											   fldi.vz[ IDX3D ] * conj( fldi.vz[ IDX3D ] ) ) / ((double) NTOTAL*NTOTAL);
				}
			}
		}
	}
	
#ifdef MPI_SUPPORT
	reduce(&q0, 1);
#endif

	q0 = q0 / (2.0*dk);
	
	// q0 is E(kc)
	
	q0 = pow(q0/kc, 0.5);
	
#ifdef _OPENMP
	#pragma omp parallel for private(i) schedule(static)
#endif
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i] = (double complex) q0 * ( 0.267 + 9.21 * exp(-3.03 * kc * pow(ik2t[i],0.5) ) );	// Original Chollet-Lesieur
//		w1[i] = 0.1 * (1.0 + 5.0*pow(k2t[i]/(kc*kc), 4.0)) * q0;		// Ponty el al 2003
		pressure[i] = w1[i];
	}
#ifdef DEBUG
	MPI_Printf("w1:\n");
	D_show_field(w1);
	
	printf("nu_t=%g\n",pow(q0/kc, 0.5));
	
#endif
	
	//printf("nu_t=%g\n",pow(q0/kc, 0.5));
#endif


#ifdef _OPENMP
	#pragma omp parallel for private(i,q0) schedule(static)
#endif
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
#ifdef SGS
		q0 = exp( - w1[i] * dt* k2t[i] );
#else
		q0 = exp( - nu * dt* k2t[i] );
#endif
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
	
	


			   