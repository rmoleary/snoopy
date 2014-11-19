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

#include <stdlib.h>
#include "common.h"
#include "error.h"
#include "gfft.h"
#include "symmetries.h"
#include "debug.h"

// Mask in reak space
double	*mask_real;


/***********************************/
/** 
	Experimental Bounadry conditions call
	This routine enforces a set of boundary conditions
	on fldi (eg rigid boundary conditions, obstacles)
	
	@param fldi: field stucture on which we apply the BC
************************************/
void boundary_c(struct Field fldi) {
	DEBUG_START_FUNC;
	
	symmetrize_walls_z(fldi);
	
	DEBUG_END_FUNC;
	return;
}

////////////////////////////////////////////////////////////
// The following routines are experimental /////////////////
// They allow the inclusion of a hard cylinder in the flow//
////////////////////////////////////////////////////////////
/******************************************
** init the real mask to introduce some 
** object in the flow
******************************************/
void init_real_mask() {
	double *x,*y,*z;
	int i,j,k;
	
	mask_real = (double *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (mask_real == NULL) ERROR_HANDLER( ERROR_CRITICAL, "no memory for mask_real profile allocation");
	
	/*******************************************************************
	** This part does not need to be modified **************************
	********************************************************************/
	// Allocate coordinate arrays
	x = (double *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (x == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for x allocation");
	
	y = (double *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (y == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for y allocation");
	
	z = (double *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (z == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for z allocation");

	// Initialize the (transposed!) arrays
	for(i = 0 ; i < NX ; i++) {
		for(j = 0 ; j < NY/NPROC ; j++) {
			for(k = 0 ; k < NZ ; k++) {
				x[k + (NZ + 2) * i + (NZ + 2) * NX * j] = - param.lx / 2 + (param.lx * i ) / NX;
				y[k + (NZ + 2) * i + (NZ + 2) * NX * j] = - param.ly / 2 + (param.ly * (j + rank * NY / NPROC)) / NY;
				z[k + (NZ + 2) * i + (NZ + 2) * NX * j] = - param.lz / 2 + (param.lz * k ) / NZ;
			}
		}
	}
	
	// Initialize the extra points (k=NZ and k=NZ+1) to zero to prevent stupid things from happening...
	for(i = 0 ; i < NX ; i++) {
		for(j = 0 ; j < NY/NPROC ; j++) {
			for(k = NZ ; k < NZ + 2 ; k++) {
				x[k + (NZ + 2) * i + (NZ + 2) * NX * j] = 0.0;
				y[k + (NZ + 2) * i + (NZ + 2) * NX * j] = 0.0;
				z[k + (NZ + 2) * i + (NZ + 2) * NX * j] = 0.0;
			}
		}
	}
	
	// Init array to zero
	for(i = 0 ; i < NX ; i++) {
		for(j = 0 ; j < NY/NPROC ; j++) {
			for(k = 0 ; k < NZ + 2 ; k++) {
				mask_real[ k + (NZ + 2) * i + (NZ + 2) * NX * j ] = 1.0;
			}
		}
	}
	// Insert a radius=1 object in the system
	
	for(i = 0 ; i < 2*NTOTAL_COMPLEX ; i++) {
		if((x[i]*x[i]+y[i]*y[i])<1.0) {
			mask_real[i] = 0.0;
		}
	}
	
	free(x);
	free(y);
	free(z);
}

/******************************************
** Enforce the real mask ******************
*******************************************/
void enforce_real_mask(struct Field fldi) {
	int i;
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i] =  fldi.vx[i];
		w2[i] =  fldi.vy[i];
		w3[i] =  fldi.vz[i];
	}
	gfft_c2r_t(w1);
	gfft_c2r_t(w2);
	gfft_c2r_t(w3);
	
	for( i = 0 ; i < 2*NTOTAL_COMPLEX ; i++) {
		wr1[i] = wr1[i] * mask_real[i] / ((double) NTOTAL);
		wr2[i] = wr2[i] * mask_real[i] / ((double) NTOTAL);
		wr3[i] = wr3[i] * mask_real[i] / ((double) NTOTAL);
	}
	
	gfft_r2c_t(wr1);
	gfft_r2c_t(wr2);
	gfft_r2c_t(wr3);
	
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		fldi.vx[i] = w1[i] * mask[i];
		fldi.vy[i] = w2[i] * mask[i];
		fldi.vz[i] = w3[i] * mask[i];
	}
	
	return;
}
