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



#include "common.h"
#ifdef MPI_SUPPORT
#include "transpose.h"
#endif
#include "debug.h"
#ifdef WITH_SHEAR


double time_shift(double t) {
	double tremap;
#ifdef TIME_DEPENDANT_SHEAR
	tremap = sin(param.omega_shear * t) / param.omega_shear;	// This is formally the angular displacement divded by param.shear= int dt S(t) /<S>
#else
	tremap = fmod(t + param.ly / (2.0 * param.shear * param.lx) , param.ly / (param.shear * param.lx)) - param.ly / (2.0 * param.shear * param.lx);
#endif
	return(tremap);
}

void remap(double complex qi[]) {
	int i, j, k;
	int nx, ny, nxtarget;
	
	DEBUG_START_FUNC;
	
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i]=0.0;
	}
	
#ifdef MPI_SUPPORT
// We have to transpose the array to get the remap properly
	transpose_complex_XY(qi,qi);
	
	for( i = 0; i < NX_COMPLEX; i++) {
		nx = fmod( i + (NX_COMPLEX / 2) ,  NX_COMPLEX ) - NX_COMPLEX / 2 ;
		for( j = 0; j < NY_COMPLEX/NPROC; j++) {
			ny = fmod( j + rank * NY_COMPLEX / NPROC + (NY_COMPLEX / 2) ,  NY_COMPLEX ) - NY_COMPLEX / 2 ;
			
			nxtarget = nx + ny;		// We have a negative shear, hence nx plus ny
			
			if( (nxtarget > -NX_COMPLEX / 2) & (nxtarget < NX_COMPLEX/2)) {
			
				if ( nxtarget <0 ) nxtarget = nxtarget + NX_COMPLEX;
			
				for( k = 0; k < NZ_COMPLEX; k++) {
					w1[k + NZ_COMPLEX * nxtarget + NZ_COMPLEX * NX_COMPLEX * j] = qi[ k + i * NZ_COMPLEX + j * NZ_COMPLEX * NX_COMPLEX];
				}
			}
		}
	}
	
	// transpose back
	transpose_complex_YX(w1,w1);

#else
	for( i = 0; i < NX_COMPLEX; i++) {
		nx = fmod( i + (NX_COMPLEX / 2) ,  NX_COMPLEX ) - NX_COMPLEX / 2 ;
		for( j = 0; j < NY_COMPLEX; j++) {
			ny = fmod( j + (NY_COMPLEX / 2) ,  NY_COMPLEX ) - NY_COMPLEX / 2 ;
			
			nxtarget = nx + ny;		// We have a negative shear, hence nx plus ny
			
			if( (nxtarget > -NX_COMPLEX / 2) & (nxtarget < NX_COMPLEX/2)) {
			
				if ( nxtarget <0 ) nxtarget = nxtarget + NX_COMPLEX;
			
				for( k = 0; k < NZ_COMPLEX; k++) {
					w1[k + NZ_COMPLEX * j + NZ_COMPLEX * NY_COMPLEX * nxtarget] = qi[ IDX3D ];
				
				}
			}
		}
	}
#endif

	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		qi[i] = w1[i] * mask[i];
	}

	DEBUG_END_FUNC;
	
	return;
}

void kvolve(const double tremap) {
	int i, j, k;
#ifdef _OPENMP
	#pragma omp parallel for private(i,j,k) schedule(static)	
#endif	
	for( i = 0; i < NX_COMPLEX/NPROC; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			for( k = 0; k < NZ_COMPLEX; k++) {
				kxt[ IDX3D ] = kx[ IDX3D ] + tremap * param.shear * ky[ IDX3D ];
			
				k2t[ IDX3D ] = kxt[IDX3D] * kxt[IDX3D] +
						       ky[IDX3D] * ky[IDX3D]+
							   kz[IDX3D] * kz[IDX3D];
							  
				if ( k2t[IDX3D] == 0.0 ) ik2t[IDX3D] = 1.0;
				else	ik2t[IDX3D] = 1.0 / k2t[IDX3D];
			}
		}
	}

	return;
}
#endif
