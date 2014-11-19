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
#include "particles.h"
#include "ltide.h"

#ifndef COMPRESSIBLE

/***************************************************************/
/**
	Compute the right hand side of the INCOMPRESSIBLE dynamical equation
	
	@param dfldo: (output) right hand side of the dynamical equation
	@param fldi: (input) current status of the flow
	@param t: current time of the simulation
	@param tremap: current remap time (only when shear is on)
	@param dt: current timestep size
*/
/***************************************************************/
void timestep( struct Field dfldo,
			   struct Field fldi,
			   const double t,
			   const double tremap,
			   const double dt) {
			   
	int i;
	double complex q0,q1;
	double qr0;
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



#ifdef ELSASSER_FORMULATION
/******************************************
** ELSASSER variable formulation **********
** To be used 
*******************************************/

// Solve the MHD equations using Elsasser fields
#ifdef _OPENMP
	#pragma omp parallel for private(i) schedule(static)	
#endif
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i] =  fldi.vx[i]+fldi.bx[i];
		w2[i] =  fldi.vy[i]+fldi.by[i];
		w3[i] =  fldi.vz[i]+fldi.bz[i];
		
		w4[i] =  fldi.vx[i]-fldi.bx[i];
		w5[i] =  fldi.vy[i]-fldi.by[i];
		w6[i] =  fldi.vz[i]-fldi.bz[i];
	}
	
	// These fields should have no divergence.
	// When shear is on, however, divergence is conserved up to the timeintegrator precision.
	// Let's clean it.
	projector(w1,w2,w3);
	projector(w4,w5,w6);

	gfft_c2r_t(w1);
	gfft_c2r_t(w2);
	gfft_c2r_t(w3);
	
	gfft_c2r_t(w4);
	gfft_c2r_t(w5);
	gfft_c2r_t(w6);
	
// Compute the Elsasser tensor

#ifdef _OPENMP
	#pragma omp parallel for private(i) schedule(static)	
#endif
	for( i = 0 ; i < 2*NTOTAL_COMPLEX ; i++) {
		wr7[i]  = wr1[i] * wr4[i] / ((double) NTOTAL*NTOTAL);
		wr8[i]  = wr1[i] * wr5[i] / ((double) NTOTAL*NTOTAL);
		wr9[i]  = wr1[i] * wr6[i] / ((double) NTOTAL*NTOTAL);
		wr10[i] = wr2[i] * wr4[i] / ((double) NTOTAL*NTOTAL);
		wr11[i] = wr2[i] * wr5[i] / ((double) NTOTAL*NTOTAL);
		wr12[i] = wr2[i] * wr6[i] / ((double) NTOTAL*NTOTAL);
		wr13[i] = wr3[i] * wr4[i] / ((double) NTOTAL*NTOTAL);
		wr14[i] = wr3[i] * wr5[i] / ((double) NTOTAL*NTOTAL);
		wr15[i] = wr3[i] * wr6[i] / ((double) NTOTAL*NTOTAL);
	}
	
	gfft_r2c_t(wr7);
	gfft_r2c_t(wr8);
	gfft_r2c_t(wr9);
	gfft_r2c_t(wr10);
	gfft_r2c_t(wr11);
	gfft_r2c_t(wr12);
	gfft_r2c_t(wr13);
	gfft_r2c_t(wr14);
	gfft_r2c_t(wr15);
	
// Compute the volution of the Elssaser fields (u= ik.
#ifdef _OPENMP
	#pragma omp parallel for private(i) schedule(static)	
#endif
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		dfldo.vx[i] = - I * 0.5 * mask[i] * (
						kxt[i] * ( 2.0 * w7[i]  ) + ky[i] * ( w8[i] + w10[i]) + kz[i] * (w9[i]  + w13[i]) );
		dfldo.vy[i] = - I * 0.5 * mask[i] * (
						kxt[i] * (w10[i] + w8[i]) + ky[i] * ( 2.0   * w11[i]) + kz[i] * (w12[i] + w14[i]) );
		dfldo.vz[i] = - I * 0.5 * mask[i] * (
						kxt[i] * (w13[i] + w9[i]) + ky[i] * (w14[i] + w12[i]) + kz[i] * ( 2.0 * w15[i]  ) );
		
										
		dfldo.bx[i] = - I * 0.5 * mask[i] * (
						                            ky[i] * ( w8[i] - w10[i]) + kz[i] * (w9[i]  - w13[i]) );
		dfldo.by[i] = - I * 0.5 * mask[i] * (
						kxt[i] * (w10[i] - w8[i])                             + kz[i] * (w12[i] - w14[i]) );
		dfldo.bz[i] = - I * 0.5 * mask[i] * (
						kxt[i] * (w13[i] - w9[i]) + ky[i] * (w14[i] - w12[i])  );
		

	}
		
// Compute real(U) in case it is used later.
	for( i = 0 ; i < 2*NTOTAL_COMPLEX ; i++) {
		wr1[i] = 0.5 * (wr1[i] + wr4[i]);
		wr2[i] = 0.5 * (wr2[i] + wr5[i]);
		wr3[i] = 0.5 * (wr3[i] + wr6[i]);
	}

#else
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

	// These fields should have no divergence.
	// When shear is on, however, divergence is conserved up to the timeintegrator precision.
	// Let's clean it.
	projector(w1,w2,w3);
	
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
	
#endif


/**********************************************
** Particles (if needed) **********************
***********************************************/

#ifdef WITH_PARTICLES
	particle_step(dfldo, fldi, wr1, wr2, wr3, t, tremap, dt);
#endif

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
#ifndef ELSASSER_FORMULATION		// If Elssaser is on, MHD are already computed...

// Start with the induction equation
#ifdef _OPENMP
	#pragma omp parallel for private(i) schedule(static)	
#endif
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w4[i] =  fldi.bx[i];
		w5[i] =  fldi.by[i];
		w6[i] =  fldi.bz[i];
	}

	// These fields should have no divergence.
	// When shear is on, however, divergence is conserved up to the timeintegrator precision.
	// Let's clean it.
	projector(w4,w5,w6);
	
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
/*************************************************
/** Braginskii viscosity.               ***********
/*************************************************/

// Only to be used without the ELSASSER formulation (otherwise we have
// to calculate the maxwell stress twice

// Compute the b_i b_j tensor

#ifdef WITH_BRAGINSKII

#ifdef _OPENMP
	#pragma omp parallel for private(i,qr0) schedule(static)	
#endif
	for( i = 0 ; i < 2*NTOTAL_COMPLEX ; i++) {
		
		qr0 = wr4[i] * wr4[i] + wr5[i] * wr5[i] + wr6[i] * wr6[i];	// norm(B)^2
		
		wr10[i] = wr4[i] * wr4[i] / qr0;
		wr11[i] = wr5[i] * wr5[i] / qr0;
		wr12[i] = wr6[i] * wr6[i] / qr0;
		wr13[i] = wr4[i] * wr5[i] / qr0;
		wr14[i] = wr4[i] * wr6[i] / qr0;
		wr15[i] = wr5[i] * wr6[i] / qr0;
	}
	
// Compute the stress tensor
#ifdef _OPENMP
	#pragma omp parallel for private(i) schedule(static)
#endif
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i] = I * kxt[i] * fldi.vx[i];
		w2[i] = I * ky[i] * fldi.vy[i];
		w3[i] = I * kz[i] * fldi.vz[i];
		w4[i] = I * ( kxt[i] * fldi.vy[i] + ky[i] * fldi.vx[i] );
		w5[i] = I * ( kxt[i] * fldi.vz[i] + kz[i] * fldi.vx[i] );
		w6[i] = I * ( ky[i]  * fldi.vz[i] + kz[i] * fldi.vy[i] );
	}
	
	gfft_c2r_t(w1);
	gfft_c2r_t(w2);
	gfft_c2r_t(w3);
	gfft_c2r_t(w4);
	gfft_c2r_t(w5);
	gfft_c2r_t(w6);
	
// Compute the pressure anisotropy

#ifdef _OPENMP
	#pragma omp parallel for private(i) schedule(static)
#endif
	for( i = 0 ; i < 2*NTOTAL_COMPLEX ; i++) {
		wr7[i] = (wr10[i] * wr1[i] 
				+ wr11[i] * wr2[i]
				+ wr12[i] * wr3[i]
				+ wr13[i] * (wr4[i] - param.shear*((double) NTOTAL))		// Take into account the background shear in the stress tensor
				+ wr14[i] * wr5[i]
				+ wr15[i] * wr6[i]
				) / ((double) NTOTAL);
	}
	
	
// Compute the viscous stress tensor
#ifdef _OPENMP
	#pragma omp parallel for private(i) schedule(static)
#endif
	for( i = 0 ; i < 2*NTOTAL_COMPLEX ; i++) {
		wr10[i] = wr10[i] * wr7[i];
		wr11[i] = wr11[i] * wr7[i];
		wr12[i] = wr12[i] * wr7[i];
		wr13[i] = wr13[i] * wr7[i];
		wr14[i] = wr14[i] * wr7[i];
		wr15[i] = wr15[i] * wr7[i];
	}
	
	gfft_r2c_t(wr10);
	gfft_r2c_t(wr11);
	gfft_r2c_t(wr12);
	gfft_r2c_t(wr13);
	gfft_r2c_t(wr14);
	gfft_r2c_t(wr15);

	
#ifdef _OPENMP
	#pragma omp parallel for private(i) schedule(static)
#endif
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		dfldo.vx[i] += 3.0 * mask[i] * I / param.reynolds_B * ( kxt[i] * w10[i] + ky[i] * w13[i] + kz[i] * w14[i] );
		dfldo.vy[i] += 3.0 * mask[i] * I / param.reynolds_B * ( kxt[i] * w13[i] + ky[i] * w11[i] + kz[i] * w15[i] );
		dfldo.vz[i] += 3.0 * mask[i] * I / param.reynolds_B * ( kxt[i] * w14[i] + ky[i] * w15[i] + kz[i] * w12[i] );
	}

#endif

#endif

/*************************************************/
/** HALL EFFECT **********************************/
/*************************************************/

#ifdef WITH_HALL

#ifdef _OPENMP
	#pragma omp parallel for private(i) schedule(static)	
#endif
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i] = I * (ky[i]  * fldi.bz[i] - kz[i]  * fldi.by[i]);
		w2[i] = I * (kz[i]  * fldi.bx[i] - kxt[i] * fldi.bz[i]);
		w3[i] = I * (kxt[i] * fldi.by[i] - ky[i]  * fldi.bx[i]);
		
		w4[i] =  fldi.bx[i];
		w5[i] =  fldi.by[i];
		w6[i] =  fldi.bz[i];
	}

	// These fields should have no divergence.
	// When shear is on, however, divergence is conserved up to the timeintegrator precision.
	// Let's clean it.
	projector(w4,w5,w6);
	
	gfft_c2r_t(w1);
	gfft_c2r_t(w2);
	gfft_c2r_t(w3);
	gfft_c2r_t(w4);
	gfft_c2r_t(w5);
	gfft_c2r_t(w6);
	
	// J is in w1-w3, B is in w4-w6
	// q0 is the Hall parameter
	q0 = 1.0 / param.x_hall;
	
#ifdef _OPENMP
	#pragma omp parallel for private(i) schedule(static)	
#endif
	for( i = 0 ; i < 2*NTOTAL_COMPLEX ; i++) {
		wr7[i] = -q0*(wr2[i] * wr6[i] - wr3[i] * wr5[i]) / ((double) NTOTAL*NTOTAL);
		wr8[i] = -q0*(wr3[i] * wr4[i] - wr1[i] * wr6[i]) / ((double) NTOTAL*NTOTAL);
		wr9[i] = -q0*(wr1[i] * wr5[i] - wr2[i] * wr4[i]) / ((double) NTOTAL*NTOTAL);
	}

	// Compute the curl of the emf to add in the induction equation.
	
	gfft_r2c_t(wr7);
	gfft_r2c_t(wr8);
	gfft_r2c_t(wr9);
	
#ifdef _OPENMP
	#pragma omp parallel for private(i) schedule(static)	
#endif
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		dfldo.bx[i] += I * mask[i] * (ky[i] * w9[i] - kz[i] * w8[i]);
		dfldo.by[i] += I * mask[i] * (kz[i] * w7[i] - kxt[i]* w9[i]);
		dfldo.bz[i] += I * mask[i] * (kxt[i]* w8[i] - ky[i] * w7[i]);
	}

#endif	//WITH_HALL



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
** EXPLICIT LINEAR DISSIPATION ******
*************************************/

#ifdef WITH_EXPLICIT_DISSIPATION
#ifdef _OPENMP
	#pragma omp parallel for private(i) schedule(static)
#endif
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		dfldo.vx[i] += - nu * k2t[i] * fldi.vx[i];
		dfldo.vy[i] += - nu * k2t[i] * fldi.vy[i];
		dfldo.vz[i] += - nu * k2t[i] * fldi.vz[i];
		
#ifdef MHD
		dfldo.bx[i] += - eta * k2t[i] * fldi.bx[i];
		dfldo.by[i] += - eta * k2t[i] * fldi.by[i];
		dfldo.bz[i] += - eta * k2t[i] * fldi.bz[i];
#endif	// MHD

#ifdef BOUSSiNESQ
		dfldo.th[i] += - nu_th * k2t[i] * fldi.th[i];
#endif	// BOUSSINESQ
	}

#endif	// WITH_EXPLICIT_DISSIPATION
	
/************************************
** PRESSURE TERMS *******************
************************************/

#ifdef _OPENMP
	#pragma omp parallel for private(i,q0,q1) schedule(static)
#endif
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
			
#ifdef WITH_SHEAR
		q0= S * ky[i] * fldi.vx[i] + kxt[i] * dfldo.vx[i] + ky[i] * dfldo.vy[i] + kz[i] * dfldo.vz[i];
#else
		q0= kxt[i] * dfldo.vx[i] + ky[i] * dfldo.vy[i] + kz[i] * dfldo.vz[i];
#endif
/* po would contain the pressure field
		if(po != NULL) {
			po[i] = - I * ik2t[i] * q0;	// Save the pressure field (if needed)
		}
*/
		dfldo.vx[i] += -kxt[i]* q0 * ik2t[i];
		dfldo.vy[i] += -ky[i] * q0 * ik2t[i];
		dfldo.vz[i] += -kz[i] * q0 * ik2t[i];
	}
	
#ifdef WITH_LINEAR_TIDE
	ltide_timestep(dfldo, fldi, t, dt);	
#endif

	return;
}

#ifdef SGS
/***************************************************************/
/**
	Subgridscale model 
	This is the Chollet-Lesieur Model (1981)
	We have nu(k)=nu_i(k)*(E(kc)/kc)^(1/2)
	nu_i(k)=0.267+9.21*exp(-3.03 kc/k)
	
	NB: this subgridscale model is applied only to the velocity field,
	even when MHD is active!
	
	This is an implicit model: fldi is modified by this routine
	
	@param fldi: (input and output) current status of the flow
	@param t: current time of the simulation
	@param dt: current timestep size
*/
/***************************************************************/
void sgs_dissipation(struct Field fldi,
			   const double t,
			   const double dt) {
			   
		// Subgrid model
				
		// Compute E(kc)
	double kc, dk;
	double q0, q1;
	int i,j,k;
	
	
	kc = 0.5*kmax;
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
	}
	
	// Apply SGS viscosity to the flow.
	
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {

		q0 = exp( - w1[i] * dt* k2t[i] );

		fldi.vx[i] = fldi.vx[i] * q0;
		fldi.vy[i] = fldi.vy[i] * q0;
		fldi.vz[i] = fldi.vz[i] * q0;
	
	}

}
#endif	// This is for SGS


/***************************************************************/
/**
	Implicit steps of the integrator (essentially linear diffusion terms)
	
	This is an implicit model: fldi is modified by this routine
	
	@param fldi: (input and output) current status of the flow
	@param t: current time of the simulation
	@param dt: current timestep size
*/
/***************************************************************/
void implicitstep(
			   struct Field fldi,
			   const double t,
			   const double dt ) {
			   
	double q0;
	int i;
	
#ifdef SGS
	sgs_dissipation( fldi, t, dt);
#endif

#ifndef WITH_EXPLICIT_DISSIPATION
#ifdef _OPENMP
	#pragma omp parallel for private(i,q0) schedule(static)
#endif
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
#ifndef SGS
		q0 = exp( - nu * dt* k2t[i] );
		fldi.vx[i] = fldi.vx[i] * q0;
		fldi.vy[i] = fldi.vy[i] * q0;
		fldi.vz[i] = fldi.vz[i] * q0;
#endif	// SGS
		
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
#endif	// WITH_EXPLICIT_DISSIPATION

#ifdef FORCING
	forcing(fldi, dt);
#endif

#ifdef WITH_PARTICLES
	particle_implicit_step( fldi, t, dt);
#endif

#ifdef WITH_LINEAR_TIDE
	ltide_implicitstep( fldi, t, dt);
#endif

	return;
}
	


#endif	// This is for COMPRESSIBLE
			   