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

#ifdef COMPRESSIBLE

/***************************************************************/
/**
	Compute the right hand side of the COMPRESSIBLE dynamical equation
	
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

	// NB: in the compressible version, fldi.vx,vy, vz are actually the linear momenta in x,y,z!
	
	int i;
	double S;
	double complex qc;
	
#ifdef WITH_SHEAR
#ifdef TIME_DEPENDANT_SHEAR
	S = param.shear * cos(param.omega_shear * t);	// This is the real shear: -dvy/dx
#else
	S = param.shear;
#endif
#endif

/******************************************
** Momentum Advection ****************
*******************************************/
	
#ifdef _OPENMP
	#pragma omp parallel for private(i) schedule(static)	
#endif
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i] =  fldi.vx[i];
		w2[i] =  fldi.vy[i];
		w3[i] =  fldi.vz[i];
		w4[i] =  fldi.d[i];
	}
	
	gfft_c2r_t(w1);
	gfft_c2r_t(w2);
	gfft_c2r_t(w3);
	gfft_c2r_t(w4);
	
// Compute density

	check_positivity(wr4);
	
	for( i = 0 ; i < 2*NTOTAL_COMPLEX ; i++) {
		wr4[i] = 1.0 / wr4[i];
	}
	
/////////////////////////////////////////////////////////////
/*
#ifdef _OPENMP
	#pragma omp parallel for private(i,q0) schedule(static)	
#endif
	for( i = 0 ; i < 2*NTOTAL_COMPLEX ; i++) {
		wr4[i] = ((double) NTOTAL) exp( wr4[i]/((double) NTOTAL) );
	}
*/
////////////////////////////////////////////////////////////
	
	// Here, we compute the Reynolds stress tensor p_ip_j/rho=v_iv_jrho
#ifdef _OPENMP
	#pragma omp parallel for private(i) schedule(static)	
#endif
	for( i = 0 ; i < 2*NTOTAL_COMPLEX ; i++) {

		wr5[i]  = wr1[i] * wr1[i] * wr4[i] / ((double) NTOTAL);
		wr6[i]  = wr2[i] * wr2[i] * wr4[i] / ((double) NTOTAL);
		wr7[i]  = wr3[i] * wr3[i] * wr4[i] / ((double) NTOTAL);
		wr8[i]  = wr1[i] * wr2[i] * wr4[i] / ((double) NTOTAL);
		wr9[i]  = wr1[i] * wr3[i] * wr4[i] / ((double) NTOTAL);
		wr10[i] = wr2[i] * wr3[i] * wr4[i] / ((double) NTOTAL);
	}
	
	gfft_r2c_t(wr5);
	gfft_r2c_t(wr6);
	gfft_r2c_t(wr7);
	gfft_r2c_t(wr8);
	gfft_r2c_t(wr9);
	gfft_r2c_t(wr10);
	
#ifdef _OPENMP
	#pragma omp parallel for private(i) schedule(static)	
#endif
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		dfldo.vx[i] = - I * mask[i] * (
					kxt[i] * w5[i] + ky[i] * w8[i] + kz[i] * w9[i] );
		dfldo.vy[i] = - I * mask[i] * (
					kxt[i] * w8[i] + ky[i] * w6[i] + kz[i] * w10[i] );
		dfldo.vz[i] = - I * mask[i] * (
					kxt[i] * w9[i] + ky[i] * w10[i] + kz[i] * w7[i] );	
	}
	
/*********************************************
** Other terms will need the velocity field **
** Computation of the velocity field        **
**********************************************/

#ifdef _OPENMP
	#pragma omp parallel for private(i) schedule(static)	
#endif
	for( i = 0 ; i < 2*NTOTAL_COMPLEX ; i++) {
		wr1[i] = wr1[i] * wr4[i];
		wr2[i] = wr2[i] * wr4[i];
		wr3[i] = wr3[i] * wr4[i];
	}

	// wr1-wr3 is the real velocity field
	
/*********************************************
** Viscous term ******************************
**********************************************/
	
	// Fourier transform the velocity field
	
	for( i = 0 ; i < 2*NTOTAL_COMPLEX ; i++) {
		wr4[i] = wr1[i];
		wr5[i] = wr2[i];
		wr6[i] = wr3[i];
	}
	
	gfft_r2c_t(wr4);
	gfft_r2c_t(wr5);
	gfft_r2c_t(wr6);
	
	// Add the viscous term to the linear momentum
	
#ifdef _OPENMP
	#pragma omp parallel for private(i,qc) schedule(static)	
#endif
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		qc = 1.0 / 3.0 * (kxt[i] * w4[i] + ky[i] * w5[i] + kz[i] * w6[i]);
		dfldo.vx[i] -= nu * mask[i] * (k2t[i] * w4[i] + kxt[i] * qc);
		dfldo.vy[i] -= nu * mask[i] * (k2t[i] * w5[i] + ky[i]  * qc);
		dfldo.vz[i] -= nu * mask[i] * (k2t[i] * w6[i] + kz[i]  * qc);
	}
	
	

/*********************************************
**** MHD Terms (if needed)   *****************
*********************************************/
#ifdef MHD

	// We first build up the velocity field from the momentum
	
	
	// do the induction equation
	
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
		wr7[i] = (wr2[i] * wr6[i] - wr3[i] * wr5[i]) / ((double) NTOTAL);
		wr8[i] = (wr3[i] * wr4[i] - wr1[i] * wr6[i]) / ((double) NTOTAL);
		wr9[i] = (wr1[i] * wr5[i] - wr2[i] * wr4[i]) / ((double) NTOTAL);
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

	// The 0.5 factor is here to take into account the magnetic pressure term -B^2/2 delta_ij
#ifdef _OPENMP
	#pragma omp parallel for private(i,qc) schedule(static)
#endif
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		qc = 0.5 * (w1[i] + w2[i] + w3[i]);
		dfldo.vx[i] += I * mask[i] * (kxt[i] * (w1[i]-qc) +     ky[i] * w7[i]      + kz[i] * w8[i]);
		dfldo.vy[i] += I * mask[i] * (kxt[i] * w7[i]      +     ky[i] * (w2[i]-qc) + kz[i] * w9[i]);
		dfldo.vz[i] += I * mask[i] * (kxt[i] * w8[i]      +     ky[i] * w9[i]      + kz[i] * (w3[i]-qc));
	}
	
#endif

/*********************************************
**** Continuity equation     *****************
*********************************************/
	
#ifdef _OPENMP
	#pragma omp parallel for private(i) schedule(static)	
#endif
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		dfldo.d[i] = -I * mask[i] * (kxt[i] * fldi.vx[i] + ky[i] * fldi.vy[i] + kz[i] * fldi.vz[i]);
	}
	
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

/**************************************
** Thermal Pressure *******************
***************************************/

// NB: magnetic pressure was already included in the lorentz force calculation
#ifdef _OPENMP
	#pragma omp parallel for private(i) schedule(static)
#endif
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		dfldo.vx[i] += -I * mask[i] * kxt[i] * param.cs * param.cs * fldi.d[i];
		dfldo.vy[i] += -I * mask[i] * ky[i]  * param.cs * param.cs * fldi.d[i];
		dfldo.vz[i] += -I * mask[i] * kz[i]  * param.cs * param.cs * fldi.d[i];
	}
		
	// Finished
	return;
}

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
			   
	double q0,lambda;
	int exponent;
	int i,j;
	
#ifdef _OPENMP
	#pragma omp parallel for private(i,q0) schedule(static)
#endif

	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {

#ifdef MHD
		q0 = exp( - eta * dt* k2t[i] );
		fldi.bx[i] = fldi.bx[i] * q0;
		fldi.by[i] = fldi.by[i] * q0;
		fldi.bz[i] = fldi.bz[i] * q0;
#endif
	}
	
// Hyperviscosity (diffusive time=grid-scale sound crossing time at the cutoff scale)

/*

	lambda=3;
	exponent=4;
#ifdef _OPENMP
	#pragma omp parallel for private(i,q0) schedule(static)
#endif
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
	
		q0 = exp( - dt* kmax*param.cs*pow(lambda*lambda*k2t[i]/(kmax*kmax),exponent) );
		
		q0 = 1.0;
		fldi.d[i] = fldi.d[i] * q0;
		q0=1.0;
		fldi.vx[i] = fldi.vx[i] * q0;
		fldi.vy[i] = fldi.vy[i] * q0;
		fldi.vz[i] = fldi.vz[i] * q0;
		
#ifdef MHD
		fldi.bx[i] = fldi.bx[i] * q0;
		fldi.by[i] = fldi.by[i] * q0;
		fldi.bz[i] = fldi.bz[i] * q0;
#endif
	}

*/

#ifdef FORCING
	forcing(fldi, dt);
#endif

#ifdef WITH_PARTICLES
	particle_implicit_step( fldi, t, dt);
#endif

	return;
}

#endif
	