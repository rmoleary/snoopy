#include "common.h"
#include "gfft.h"

#ifdef DEBUG
#include "debug.h"
#endif


/** Allow one to init a structure in real space using ordinary defined x,y,z coordinates */

void init_spatial_structure() {
	PRECISION *x,*y,*z;
	int i,j,k;
	
	/*******************************************************************
	** This part does not need to be modified **************************
	********************************************************************/
	// Allocate coordinate arrays
	x = (PRECISION *) fftw_malloc( sizeof(PRECISION complex) * NTOTAL_COMPLEX);
	if (x == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for x allocation");
	
	y = (PRECISION *) fftw_malloc( sizeof(PRECISION complex) * NTOTAL_COMPLEX);
	if (y == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for y allocation");
	
	z = (PRECISION *) fftw_malloc( sizeof(PRECISION complex) * NTOTAL_COMPLEX);
	if (z == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for z allocation");

	// Initialize the arrays
	
	for(i = 0 ; i < NX/NPROC ; i++) {
		for(j = 0 ; j < NY ; j++) {
			for(k = 0 ; k < NZ ; k++) {
				x[k + (NZ + 2) * j + (NZ + 2) * NY * i] = - LX / 2 + (LX * (i + rank * NX / NPROC)) / NX;
				y[k + (NZ + 2) * j + (NZ + 2) * NY * i] = - LY / 2 + (LY * j ) / NY;
				z[k + (NZ + 2) * j + (NZ + 2) * NY * i] = - LZ / 2 + (LZ * k ) / NZ;
			}
		}
	}
	
	// Initialize the extra points (k=NZ and k=NZ+1) to zero to prevent stupid things from happening...
	for(i = 0 ; i < NX/NPROC ; i++) {
		for(j = 0 ; j < NY ; j++) {
			for(k = NZ ; k < NZ + 2 ; k++) {
				x[k + (NZ + 2) * j + (NZ + 2) * NY * i] = 0.0;
				y[k + (NZ + 2) * j + (NZ + 2) * NY * i] = 0.0;
				z[k + (NZ + 2) * j + (NZ + 2) * NY * i] = 0.0;
			}
		}
	}
	
	// Init work array to zero
	for(i = 0 ; i < NX/NPROC ; i++) {
		for(j = 0 ; j < NY ; j++) {
			for(k = 0 ; k < NZ + 2 ; k++) {
				wr1[k + (NZ + 2) * j + (NZ + 2) * NY * i] = 0.0;
				wr2[k + (NZ + 2) * j + (NZ + 2) * NY * i] = 0.0;
				wr3[k + (NZ + 2) * j + (NZ + 2) * NY * i] = 0.0;
				wr4[k + (NZ + 2) * j + (NZ + 2) * NY * i] = 0.0;
				wr5[k + (NZ + 2) * j + (NZ + 2) * NY * i] = 0.0;
				wr6[k + (NZ + 2) * j + (NZ + 2) * NY * i] = 0.0;
			}
		}
	}
	
	/*******************************************************************
	** This part can be modified              **************************
	********************************************************************/
	
	// The velocity field vx,vy,vz is stored in wr1,wr2,wr3
	// The magnetic field bx,by,bz is stored in wr4,wr5,wr6 (ignored if MHD is not set)
	
	for(i = 0 ; i < 2*NTOTAL_COMPLEX ; i++) {
		// Example: init a flux tube in the x direction+a vertical displacement
		wr4[i] = exp(-(y[i]*y[i]+z[i]*z[i])*20.0);
		wr3[i] = 0.5*cos(x[i] * 2.0 * M_PI);
	}
	
	/*******************************************************************
	** This part does not need to be modified **************************
	********************************************************************/
	// Fourier transform everything
	gfft_r2c(wr1);
	gfft_r2c(wr2);
	gfft_r2c(wr3);
	gfft_r2c(wr4);
	gfft_r2c(wr5);
	gfft_r2c(wr6);
	
	// Transfer data in the relevant array (including dealiasing mask)
	for(i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		fld.vx[i] = w1[i] * mask[i];
		fld.vy[i] = w2[i] * mask[i];
		fld.vz[i] = w3[i] * mask[i];
#ifdef MHD
		fld.bx[i] = w4[i] * mask[i];
		fld.by[i] = w5[i] * mask[i];
		fld.bz[i] = w6[i] * mask[i];
#endif
	}
	
	// Remove divergence (if any)
	projector(fld.vx,fld.vy,fld.vz);
#ifdef MHD
	projector(fld.bx,fld.by,fld.bz);
#endif
	// free memory
	fftw_free(x);
	fftw_free(y);
	fftw_free(z);
	
	//done
	return;
}


// Caution: Not coded for MPI!
void init_vortex(PRECISION complex wzf[]) {
	const PRECISION a = 0.04;
	const PRECISION b = 0.16;
	
	int i,j,k;
	
#ifdef SUPPORT_MPI
	ERROR_HANDLER( ERROR_CRITICAL, "No MPI Support for init_vortex");
#endif
	PRECISION w0, x, y;
	PRECISION chi;
	
	chi = b / a;
	w0 = 1.0/chi*(chi + 1)/(chi-1.0);			// According to Kida!
	
	for(i = 0 ; i < NX ; i++) {
		x = - LX / 2 + (LX * i) / NX;
		for(j = 0 ; j < NY ; j++) {
			y = - LY / 2 + (LY * j) / NY;
			for(k = 0 ; k < NZ ; k++) {
				if(x * x / (a * a) + y * y / (b * b) < 1) {
					// we are in the vortex
					wr1[k + j*(NZ+2) + (NZ+2) * NY * i] = -w0;
				}
				else {
					wr1[k + j*(NZ+2) + (NZ+2) * NY * i] = 0.0;
				}
			}
		}
	}
	
	// transform
	gfft_r2c(wr1);
	
	for(i = 0 ; i < NX_COMPLEX ; i++) {
		for(j = 0 ; j < NY_COMPLEX ; j++) {
			for(k = 0 ; k < NZ_COMPLEX ; k++) {
				fld.vx[ IDX3D ] +=  I * ky[i] * w1[i] * ik2t[i];
				fld.vy[ IDX3D ] += -I * kxt[i] * w1[i] * ik2t[i];
			}
		}
	}
	
	// done
	return;
}

/** Init the flow arrays... */	
void init_flow() {
	int i,j,k,k0;
	
	// Initialise vectors to 0
	
	for( i = 0; i < NX_COMPLEX/NPROC; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			for( k = 0; k < NZ_COMPLEX; k++) {
				fld.vx[ IDX3D ] = 0.0;
				fld.vy[ IDX3D ] = 0.0;
				fld.vz[ IDX3D ] = 0.0;
				
#ifdef BOUSSINESQ
				fld.th[ IDX3D ] = 0.0;
#endif
#ifdef MHD
				fld.bx[ IDX3D ] = 0.0;
				fld.by[ IDX3D ] = 0.0;
				fld.bz[ IDX3D ] = 0.0;
#endif
			}
		}
	}
		
// Init the spatial structure
	init_spatial_structure();	
		
	if(rank==0) k0=1;
	else k0=0;
// Add some noise	
	if(rank==0) {
		for( i = 0; i < 4; i++) {
			for( j = 0; j < 4; j++) {
				for( k = k0; k < 4; k++) {
					fld.vx[ IDX3D ] += PER_AMPLITUDE * mask[IDX3D] * randm() * cexp( I * 2.0*M_PI*randm() );
					fld.vy[ IDX3D ] += PER_AMPLITUDE * mask[IDX3D] * randm() * cexp( I * 2.0*M_PI*randm() );
					fld.vz[ IDX3D ] += PER_AMPLITUDE * mask[IDX3D] * randm() * cexp( I * 2.0*M_PI*randm() );
				}
			}
	
		}
#ifdef MHD
	
		for( i = 0; i < 4; i++) {
			for( j = 0; j < 4; j++) {
				for( k = k0; k < 4; k++) {
					fld.bx[ IDX3D ] += PER_AMPLITUDE * mask[IDX3D] * randm() * cexp( I * 2.0*M_PI*randm() );
					fld.by[ IDX3D ] += PER_AMPLITUDE * mask[IDX3D] * randm() * cexp( I * 2.0*M_PI*randm() );
					fld.bz[ IDX3D ] += PER_AMPLITUDE * mask[IDX3D] * randm() * cexp( I * 2.0*M_PI*randm() );
				}
			}
		}
		// Init the mean fields.
		
		/*
		fld.bx[0] = BX0 * ((double) NTOTAL);
		fld.by[0] = BY0 * ((double) NTOTAL);
		fld.bz[0] = BZ0 * ((double) NTOTAL);
		*/
#endif
	}




	
#ifdef DEBUG
	MPI_Printf("Initflow:\n");
	D_show_all(fld);
	MPI_Printf("**************************************************************************************\n");
#endif	
	
	projector(fld.vx,fld.vy,fld.vz);
#ifdef MHD
	projector(fld.bx,fld.by,fld.bz);
#endif	
	return;
}
	
	