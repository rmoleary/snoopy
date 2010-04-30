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

#ifdef WITH_LINEAR_TIDE

void ltide_timestep(struct Field dfldo,
			   struct Field fldi,
			   const double t,
			   const double dt) {
			   
	int i;
	double q0;
	double A;
	
/************************************************
** This are terms associated with Goodman & Oh **
*************************************************/
	A = - cos( param.omega_shear * t );	// This is A in Goodman's notations
	
#ifdef _OPENMP
	#pragma omp parallel for private(i) schedule(static)
#endif
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		dfldo.tvx[i] = 0.0;
		dfldo.tvy[i] = - 2.0 * A * fldi.vx[i]; 
		dfldo.tvz[i] = 0.0;
		
		// Pressure for the test field
		q0= kxt[i] * dfldo.tvx[i] + ky[i] * dfldo.tvy[i] + kz[i] * dfldo.tvz[i];
		
		dfldo.tvx[i] += -kxt[i]* q0 * ik2t[i];
		dfldo.tvy[i] += -ky[i] * q0 * ik2t[i];
		dfldo.tvz[i] += -kz[i] * q0 * ik2t[i];
		
	}
	
	return;
	
}

	
/************************************
** Implicit steps for tidal disturbances
*************************************/	
void ltide_implicitstep(
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
		
		fldi.tvx[i] = fldi.tvx[i] * q0;
		fldi.tvy[i] = fldi.tvy[i] * q0;
		fldi.tvz[i] = fldi.tvz[i] * q0;
		
	}

	return;
}

#endif
