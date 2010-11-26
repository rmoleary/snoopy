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
#include "debug.h"

/**************************************************/
/* This code contains several routines to enforce symmetries
in Fourier space. Since some of these symmetries can be seen
as boundary conditions in real space, boundary_c (found in boundary.c)
calls several of these routines.
/**************************************************/


/*****************************************/
/** Symmetrize the complex space according
**  to the symetries of the real tranform
** @param wi  Field to be symmetrized
*/
/*******************************************/

void symmetrize_complex(double complex wi[]) {
	// Symmetrize an array
	
	int i,j,k,index;
	int idx2d, idx2dconj;
	
	double complex * zplane;
	double complex q0;
	
	DEBUG_START_FUNC;
	
#ifdef MPI_SUPPORT
	// Allocate zplane
	if(rank==0) {
		// 2 cases: if zplane fits in one of the temporary array, we're fine, and we can use this array
		// Otherwise, we have to allocate an array for this specific task
		
		if( NZ_COMPLEX >= NPROC)
			zplane = w2;
		else {
			zplane = (double complex *) fftw_malloc( sizeof(double complex) * NX_COMPLEX * NY_COMPLEX);
			if (zplane == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for zplane allocation");
		}
	}
#endif

#ifndef WITH_2D
	// put kz=0 plane in w1
	k=0;
	index=0;
	for(i=0 ; i< NX_COMPLEX/NPROC ; i++) {
		for(j=0 ; j < NY_COMPLEX ; j++ ) {
			w1[index] = wi[ IDX3D ];
			index++;
		}
	}
#else
// put ky=0 line in w1
	j=0;
	k=0;
	index=0;
	for(i=0 ; i< NX_COMPLEX/NPROC ; i++) {
		w1[index] = wi[ IDX3D ];
		index++;
	}
#endif


#ifdef MPI_SUPPORT
	// construct the full kz=0 plane in zplane on rank=0 process
	MPI_Gather(w1,	   2 * NX_COMPLEX * NY_COMPLEX / NPROC, MPI_DOUBLE,
			   zplane, 2 * NX_COMPLEX * NY_COMPLEX / NPROC, MPI_DOUBLE,
			   0, MPI_COMM_WORLD);
#else
	zplane = w1;
#endif

#ifndef WITH_2D
	// Now let's do the nasty symmetrization bit
	if(rank==0) {
		// kz=0, ky!=0
		for (i = 0 ; i < NX_COMPLEX ; i++) {
			for(j = 1 ; j < NY_COMPLEX/2 ; j++) {
				idx2d     = j +  NY_COMPLEX * i;
				if(i>0) 
					idx2dconj = (NY_COMPLEX - j) + NY_COMPLEX * (NX_COMPLEX - i); // where the complex conjugate should be
				else
					idx2dconj = (NY_COMPLEX - j) ;	// special case when i=0	
					
				q0=0.5*( zplane[ idx2d ] + conj(zplane[ idx2dconj ]) );
				
				zplane[ idx2d ] = q0;
				zplane[ idx2dconj ] = conj(q0);
			}
		}
		// kz=0, ky=0
		
		for (i=1 ; i < NX_COMPLEX / 2 ; i++) {
			idx2d     = NY_COMPLEX * i; // That's actually IDX3D
			idx2dconj = NY_COMPLEX * (NX_COMPLEX - i); // where the complex conjugate should be
			
			q0=0.5*( zplane[ idx2d ] + conj(zplane[ idx2dconj ]) );
				
			zplane[ idx2d ] = q0;
			zplane[ idx2dconj ] = conj(q0);
		}
	}

#else

	if(rank==0) {

		// ky=0
		
		for (i=1 ; i < NX_COMPLEX / 2 ; i++) {
			idx2d     = NY_COMPLEX * i; // That's actually IDX3D
			idx2dconj = NY_COMPLEX * (NX_COMPLEX - i); // where the complex conjugate should be
			
			q0=0.5*( zplane[ idx2d ] + conj(zplane[ idx2dconj ]) );
				
			zplane[ idx2d ] = q0;
			zplane[ idx2dconj ] = conj(q0);
		}
	}
#endif
	
	// Wait until the symmetrization is finished.
#ifdef MPI_SUPPORT	
	MPI_Barrier(MPI_COMM_WORLD);
	
	//Send it back
	MPI_Scatter( zplane, 2 * NX_COMPLEX * NY_COMPLEX / NPROC, MPI_DOUBLE,
				 w1,     2 * NX_COMPLEX * NY_COMPLEX / NPROC, MPI_DOUBLE,
				 0, MPI_COMM_WORLD);
#endif
	// No need to translate that back when no MPI is available since zplane=w1

	// Put it back in the array
#ifndef WITH_2D
	k=0;
	index=0;
	for(i=0 ; i< NX_COMPLEX/NPROC ; i++) {
		for(j=0 ; j < NY_COMPLEX ; j++) {
			wi[ IDX3D ] = w1[index];
			index++;
		}
	}
#else
	j=0;
	k=0;
	index=0;
	for(i=0 ; i< NX_COMPLEX/NPROC ; i++) {
		wi[IDX3D] = w1[index];
		index++;
	}
#endif
	// that's all folks...
	
#ifdef MPI_SUPPORT
	if(rank==0) {
		if(NZ_COMPLEX < NPROC)
			fftw_free(zplane);
	}
#endif
	
	DEBUG_END_FUNC;

	return;
}

/*****************************************
/** Enforce Symmetries of field fld     
** useful if numerical noise produces spurious
** growth of mean velocity field due to
** some linear source terms.
** Might be useful if N^2<0 or kappa^2<0
** @param fld field needed to be symmetrized
*/
/******************************************/

void enforce_complex_symm(struct Field fldi) {
	DEBUG_START_FUNC;
	
	// Enforce symmetries of complex plane
	symmetrize_complex(fldi.vx);
	symmetrize_complex(fldi.vy);
	symmetrize_complex(fldi.vz);
#ifdef BOUSSINESQ
	symmetrize_complex(fldi.th);
#endif
#ifdef MHD
  symmetrize_complex(fldi.bx);
  symmetrize_complex(fldi.by);
  symmetrize_complex(fldi.bz);
#endif

	// Remove mean field (noise is generated by the FFTs performed by nonlinear terms)
	if(rank==0) {
		fldi.vx[ 0 ] = 0.0;
		fldi.vy[ 0 ] = 0.0;
		fldi.vz[ 0 ] = 0.0;
#ifdef BOUSSINESQ
		fldi.th[ 0 ] = 0.0;
#endif
	}
	
	DEBUG_END_FUNC;
	
	return;
}


/*****************************************/
/** Symmetrize the complex space assuming
**  wi is a sine in the z direction
/*******************************************/
void symm_sin_z(double complex wi[]) {
	int i,j,k,itarget,jtarget;
	double complex q0, q1;

#ifdef MPI_SUPPORT
	double complex * zplane;
	int index;
	
	DEBUG_START_FUNC;
	
	if(rank==0) {
		zplane = (double complex *) fftw_malloc( sizeof(double complex) * NX_COMPLEX * NY_COMPLEX);
		if (zplane == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for zplane allocation");
	}
		
	// Loop on all kzs
	for(k = 0 ; k < NZ_COMPLEX ; k++) {
		//copy the k plane into w1
		index=0;
		for(i=0 ; i< NX_COMPLEX/NPROC ; i++) {
			for(j=0 ; j < NY_COMPLEX ; j++ ) {
				w1[index] = wi[ IDX3D ];
				index++;
			}
		}
		// construct the full kz=k plane in zplane on rank=0 process
		MPI_Gather(w1,	   2 * NX_COMPLEX * NY_COMPLEX / NPROC, MPI_DOUBLE,
				   zplane, 2 * NX_COMPLEX * NY_COMPLEX / NPROC, MPI_DOUBLE,
				   0, MPI_COMM_WORLD);
			   
		// loop in rank=0 to symmetrize the thing...
		if(rank==0) {
			for( i = 0; i <= NX_COMPLEX / 2; i++) {
				if(i!=0)
					itarget = NX_COMPLEX - i;
				else
					itarget = i;
				for( j = 0; j < NY_COMPLEX; j++) {
					if(j!=0)
						jtarget = NY_COMPLEX - j;
					else
						jtarget = j;
					
					q0 = zplane[ j + i * NY_COMPLEX];					
					q1 = conj(zplane[ jtarget  + itarget * NY_COMPLEX]);
					q0 = 0.5*(q0-q1);
					zplane[ j + i * NY_COMPLEX] = q0;
					zplane[ jtarget  + itarget * NY_COMPLEX] = -conj(q0);
				}
			}
		}
		
		// everyone wait here
		MPI_Barrier(MPI_COMM_WORLD);
		//Send it back
		MPI_Scatter( zplane, 2 * NX_COMPLEX * NY_COMPLEX / NPROC, MPI_DOUBLE,
					 w1,     2 * NX_COMPLEX * NY_COMPLEX / NPROC, MPI_DOUBLE,
					 0, MPI_COMM_WORLD);
		
		// put it back
		index=0;
		for(i=0 ; i< NX_COMPLEX/NPROC ; i++) {
			for(j=0 ; j < NY_COMPLEX ; j++ ) {
				wi[ IDX3D ] = w1[index];
				index++;
			}
		}
		// end of k-loop
	}
	
	if(rank==0) {
		free(zplane);
	}

#else
	
	for( i = 0; i <= NX_COMPLEX / 2; i++) {
		if(i!=0)
			itarget = NX_COMPLEX - i;
		else
			itarget = i;
			
		//MPI_Printf("i=%d, itarget=%d\n",i,itarget);
		for( j = 0; j < NY_COMPLEX; j++) {
			if(j!=0)
				jtarget = NY_COMPLEX - j;
			else
				jtarget = j;

			for( k = 0; k < NZ_COMPLEX; k++) {
				
				//MPI_Printf("kx1=%g, kx2=%g, ky1=%g, ky2=%g, kz1=%g, kz2=%g\n",kx[ k + j * NZ_COMPLEX + i * NZ_COMPLEX * NY_COMPLEX],kx[ k + j * NZ_COMPLEX + itarget * NZ_COMPLEX * NY_COMPLEX],ky[ k + j * NZ_COMPLEX + i * NZ_COMPLEX * NY_COMPLEX],ky[ k + j * NZ_COMPLEX + itarget * NZ_COMPLEX * NY_COMPLEX],kz[ k + j * NZ_COMPLEX + i * NZ_COMPLEX * NY_COMPLEX],kz[ k + j * NZ_COMPLEX + itarget * NZ_COMPLEX * NY_COMPLEX]);
				q0=wi[ k + j * NZ_COMPLEX + i * NZ_COMPLEX * NY_COMPLEX];
				q1=conj(wi[ k + jtarget * NZ_COMPLEX + itarget * NZ_COMPLEX * NY_COMPLEX]);
				q0 = 0.5*(q0-q1);
				wi[ k + j * NZ_COMPLEX + i * NZ_COMPLEX * NY_COMPLEX]= q0;
				wi[ k + jtarget * NZ_COMPLEX + itarget * NZ_COMPLEX * NY_COMPLEX] = -conj(q0);
			}
		}
	}
#endif

	DEBUG_END_FUNC;
	
	return;
}


/*****************************************/
/** Symmetrize the complex space assuming
**  wi is a cosine in the z direction
/*******************************************/
void symm_cos_z(double complex wi[]) {
	int i,j,k,itarget,jtarget;
	double complex q0, q1;

#ifdef MPI_SUPPORT
	double complex * zplane;
	int index;
	
	DEBUG_START_FUNC;
	
	if(rank==0) {
		zplane = (double complex *) fftw_malloc( sizeof(double complex) * NX_COMPLEX * NY_COMPLEX);
		if (zplane == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for zplane allocation");
	}
		
	// Loop on all kzs
	for(k = 0 ; k < NZ_COMPLEX ; k++) {
		//copy the k plane into w1
		index=0;
		for(i=0 ; i< NX_COMPLEX/NPROC ; i++) {
			for(j=0 ; j < NY_COMPLEX ; j++ ) {
				w1[index] = wi[ IDX3D ];
				index++;
			}
		}
		// construct the full kz=k plane in zplane on rank=0 process
		MPI_Gather(w1,	   2 * NX_COMPLEX * NY_COMPLEX / NPROC, MPI_DOUBLE,
				   zplane, 2 * NX_COMPLEX * NY_COMPLEX / NPROC, MPI_DOUBLE,
				   0, MPI_COMM_WORLD);
			   
		// loop in rank=0 to symmetrize the thing...
		if(rank==0) {
			for( i = 0; i <= NX_COMPLEX / 2; i++) {
				if(i!=0)
					itarget = NX_COMPLEX - i;
				else
					itarget = i;
				for( j = 0; j < NY_COMPLEX; j++) {
					if(j!=0)
						jtarget = NY_COMPLEX - j;
					else
						jtarget = j;
					
					q0 = zplane[ j + i * NY_COMPLEX];					
					q1 = conj(zplane[ jtarget  + itarget * NY_COMPLEX]);
					q0 = 0.5*(q0+q1);
					zplane[ j + i * NY_COMPLEX] = q0;
					zplane[ jtarget  + itarget * NY_COMPLEX] = conj(q0);
				}
			}
		}
		
		// everyone wait here
		MPI_Barrier(MPI_COMM_WORLD);
		//Send it back
		MPI_Scatter( zplane, 2 * NX_COMPLEX * NY_COMPLEX / NPROC, MPI_DOUBLE,
					 w1,     2 * NX_COMPLEX * NY_COMPLEX / NPROC, MPI_DOUBLE,
					 0, MPI_COMM_WORLD);
		
		// put it back
		index=0;
		for(i=0 ; i< NX_COMPLEX/NPROC ; i++) {
			for(j=0 ; j < NY_COMPLEX ; j++ ) {
				wi[ IDX3D ] = w1[index];
				index++;
			}
		}
		// end of k-loop
	}
	
	if(rank==0) {
		free(zplane);
	}

#else
	
	for( i = 0; i <= NX_COMPLEX / 2; i++) {
		if(i!=0)
			itarget = NX_COMPLEX - i;
		else
			itarget = i;
			
		//MPI_Printf("i=%d, itarget=%d\n",i,itarget);
		for( j = 0; j < NY_COMPLEX; j++) {
			if(j!=0)
				jtarget = NY_COMPLEX - j;
			else
				jtarget = j;

			for( k = 0; k < NZ_COMPLEX; k++) {
				
				//MPI_Printf("kx1=%g, kx2=%g, ky1=%g, ky2=%g, kz1=%g, kz2=%g\n",kx[ k + j * NZ_COMPLEX + i * NZ_COMPLEX * NY_COMPLEX],kx[ k + j * NZ_COMPLEX + itarget * NZ_COMPLEX * NY_COMPLEX],ky[ k + j * NZ_COMPLEX + i * NZ_COMPLEX * NY_COMPLEX],ky[ k + j * NZ_COMPLEX + itarget * NZ_COMPLEX * NY_COMPLEX],kz[ k + j * NZ_COMPLEX + i * NZ_COMPLEX * NY_COMPLEX],kz[ k + j * NZ_COMPLEX + itarget * NZ_COMPLEX * NY_COMPLEX]);
				q0=wi[ k + j * NZ_COMPLEX + i * NZ_COMPLEX * NY_COMPLEX];
				q1=conj(wi[ k + jtarget * NZ_COMPLEX + itarget * NZ_COMPLEX * NY_COMPLEX]);
				q0 = 0.5*(q0+q1);
				wi[ k + j * NZ_COMPLEX + i * NZ_COMPLEX * NY_COMPLEX]= q0;
				wi[ k + jtarget * NZ_COMPLEX + itarget * NZ_COMPLEX * NY_COMPLEX] = conj(q0);
			}
		}
	}
#endif

	DEBUG_END_FUNC;
	return;
}

/*****************************************/
/** Symmetrize the complex space assuming
**  wi is a cosine in the x direction
/*******************************************/

void symm_cos_x(double complex wi[]) {
	
	int i,j,k,itarget;
	double complex q0, q1;
	
	DEBUG_START_FUNC;
	
#ifdef MPI_SUPPORT
	transpose_complex_XY(wi,wi);
	for( j = 0; j < NY_COMPLEX/NPROC; j++) {
		for( i = 0; i < NX_COMPLEX / 2; i++) {
			if(i!=0)
				itarget = NX_COMPLEX - i;
			else
				itarget = i;
			for( k = 0; k < NZ_COMPLEX; k++) {
				q0=wi[ k + i * NZ_COMPLEX + j * NZ_COMPLEX * NX_COMPLEX];
				q1=wi[ k + itarget * NZ_COMPLEX + j * NZ_COMPLEX * NX_COMPLEX];
				q0 = 0.5*(q0+q1);
				wi[ k + i * NZ_COMPLEX + j * NZ_COMPLEX * NX_COMPLEX]= q0;
				wi[ k + itarget * NZ_COMPLEX + j * NZ_COMPLEX * NX_COMPLEX] = q0;
			}
		}
	}
	transpose_complex_YX(wi,wi);
#else
	for( i = 0; i < NX_COMPLEX / 2; i++) {
		if(i!=0)
			itarget = NX_COMPLEX - i;
		else
			itarget = i;
			
		//MPI_Printf("i=%d, itarget=%d\n",i,itarget);
		for( j = 0; j < NY_COMPLEX; j++) {
			for( k = 0; k < NZ_COMPLEX; k++) {
				
				//MPI_Printf("kx1=%g, kx2=%g, ky1=%g, ky2=%g, kz1=%g, kz2=%g\n",kx[ k + j * NZ_COMPLEX + i * NZ_COMPLEX * NY_COMPLEX],kx[ k + j * NZ_COMPLEX + itarget * NZ_COMPLEX * NY_COMPLEX],ky[ k + j * NZ_COMPLEX + i * NZ_COMPLEX * NY_COMPLEX],ky[ k + j * NZ_COMPLEX + itarget * NZ_COMPLEX * NY_COMPLEX],kz[ k + j * NZ_COMPLEX + i * NZ_COMPLEX * NY_COMPLEX],kz[ k + j * NZ_COMPLEX + itarget * NZ_COMPLEX * NY_COMPLEX]);
				q0=wi[ k + j * NZ_COMPLEX + i * NZ_COMPLEX * NY_COMPLEX];
				q1=wi[ k + j * NZ_COMPLEX + itarget * NZ_COMPLEX * NY_COMPLEX];
				q0 = 0.5*(q0+q1);
				wi[ k + j * NZ_COMPLEX + i * NZ_COMPLEX * NY_COMPLEX]= q0;
				wi[ k + j * NZ_COMPLEX + itarget * NZ_COMPLEX * NY_COMPLEX] = q0;
			}
		}
	}
#endif
	
	DEBUG_END_FUNC;
	
	return;
}

/*****************************************/
/** Symmetrize the complex space assuming
**  wi is a cosine in the x direction
/*******************************************/

void symm_sin_x(double complex wi[]) {
	
	int i,j,k,itarget;
	double complex q0, q1;
	
	DEBUG_START_FUNC;
	
#ifdef MPI_SUPPORT
	transpose_complex_XY(wi,wi);
	for( j = 0; j < NY_COMPLEX/NPROC; j++) {
		for( i = 0; i < NX_COMPLEX / 2; i++) {
			if(i!=0)
				itarget = NX_COMPLEX - i;
			else
				itarget = i;
			for( k = 0; k < NZ_COMPLEX; k++) {
				q0=wi[ k + i * NZ_COMPLEX + j * NZ_COMPLEX * NX_COMPLEX];
				q1=wi[ k + itarget * NZ_COMPLEX + j * NZ_COMPLEX * NX_COMPLEX];
				q0 = 0.5*(q0+q1);
				wi[ k + i * NZ_COMPLEX + j * NZ_COMPLEX * NX_COMPLEX]= q0;
				wi[ k + itarget * NZ_COMPLEX + j * NZ_COMPLEX * NX_COMPLEX] = q0;
			}
		}
	}
	transpose_complex_YX(wi,wi);
#else
	for( i = 0; i < NX_COMPLEX / 2; i++) {
		if(i!=0)
			itarget = NX_COMPLEX - i;
		else
			itarget = i;
			
		//MPI_Printf("i=%d, itarget=%d\n",i,itarget);
		for( j = 0; j < NY_COMPLEX; j++) {
			for( k = 0; k < NZ_COMPLEX; k++) {
				
				//MPI_Printf("kx1=%g, kx2=%g, ky1=%g, ky2=%g, kz1=%g, kz2=%g\n",kx[ k + j * NZ_COMPLEX + i * NZ_COMPLEX * NY_COMPLEX],kx[ k + j * NZ_COMPLEX + itarget * NZ_COMPLEX * NY_COMPLEX],ky[ k + j * NZ_COMPLEX + i * NZ_COMPLEX * NY_COMPLEX],ky[ k + j * NZ_COMPLEX + itarget * NZ_COMPLEX * NY_COMPLEX],kz[ k + j * NZ_COMPLEX + i * NZ_COMPLEX * NY_COMPLEX],kz[ k + j * NZ_COMPLEX + itarget * NZ_COMPLEX * NY_COMPLEX]);
				q0=wi[ k + j * NZ_COMPLEX + i * NZ_COMPLEX * NY_COMPLEX];
				q1=wi[ k + j * NZ_COMPLEX + itarget * NZ_COMPLEX * NY_COMPLEX];
				q0 = 0.5*(q0-q1);
				wi[ k + j * NZ_COMPLEX + i * NZ_COMPLEX * NY_COMPLEX]= q0;
				wi[ k + j * NZ_COMPLEX + itarget * NZ_COMPLEX * NY_COMPLEX] = -q0;
			}
		}
	}
#endif
	
	DEBUG_END_FUNC;
	
	return;
}

/*****************************************/
/** Symmetrize the complex space assuming
**  we have walls in the radial direction
**	is equivalent to a plane Couette flow (but
**	no spectral accuracy...)
*/
/*******************************************/

void symmetrize_walls_x(struct Field fldi) {
	DEBUG_START_FUNC;
	
	symm_sin_x(fldi.vx);
	symm_cos_x(fldi.vy);
	symm_cos_x(fldi.vz);
#ifdef MHD
	symm_cos_x(fldi.bx);
	symm_sin_x(fldi.by);
	symm_sin_x(fldi.bz);
#endif
#ifdef BOUSSINESQ
	symm_sin_x(fldi.th);
#endif
	
	enforce_complex_symm(fldi);
	
	DEBUG_END_FUNC;
	return;
}				
	
/*****************************************/
/** Symmetrize the complex space assuming
**  we have walls in the vertical direction.
*/
/*******************************************/

void symmetrize_walls_z(struct Field fldi) {
	DEBUG_START_FUNC;
	
	symm_cos_z(fldi.vx);
	symm_cos_z(fldi.vy);
	symm_sin_z(fldi.vz);
#ifdef MHD
	symm_sin_z(fldi.bx);
	symm_sin_z(fldi.by);
	symm_cos_z(fldi.bz);
#endif
#ifdef BOUSSINESQ
	symm_sin_z(fldi.th);
#endif

	enforce_complex_symm(fldi);
	
	DEBUG_END_FUNC;
	return;
}				
	


