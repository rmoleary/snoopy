#include "common.h"
#ifdef MPI_SUPPORT
#include "transpose.h"
#endif

#ifdef WITH_SHEAR


PRECISION time_shift(PRECISION t) {
	PRECISION tremap;
#ifdef TIME_DEPENDANT_SHEAR
	tremap = sin(OMEGA_SHEAR * t);
#else
	tremap = fmod(t + LY / (2.0 * SHEAR * LX) , LY / (SHEAR * LX)) - LY / (2.0 * SHEAR * LX);
#endif
	return(tremap);
}

void remap(PRECISION complex qi[]) {
	int i, j, k;
	int nx, ny, nxtarget;
	
#ifdef DEBUG
	MPI_Printf("Remap called\n");
#endif
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

	
	return;
}

void kvolve(const PRECISION tremap) {
	int i, j, k;
#pragma omp parallel private(i,j,k) num_threads ( NTHREADS )
{
		/* Compute the convolution */
	#pragma omp for schedule(static) nowait	
	for( i = 0; i < NX_COMPLEX/NPROC; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			for( k = 0; k < NZ_COMPLEX; k++) {
				kxt[ IDX3D ] = kx[ IDX3D ] + tremap * SHEAR * ky[ IDX3D ];
			
				k2t[ IDX3D ] = kxt[IDX3D] * kxt[IDX3D] +
						       ky[IDX3D] * ky[IDX3D]+
							   kz[IDX3D] * kz[IDX3D];
							  
				if ( k2t[IDX3D] == 0.0 ) ik2t[IDX3D] = 1.0;
				else	ik2t[IDX3D] = 1.0 / k2t[IDX3D];
			}
		}
	}
}
	return;
}
#endif
