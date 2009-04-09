#include "common.h"
#include "gfft.h"

#include "debug.h"



/** Allow one to init a structure in real space using ordinary defined x,y,z coordinates */

void init_SpatialStructure() {
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
		fld.vx[i] += w1[i] * mask[i];
		fld.vy[i] += w2[i] * mask[i];
		fld.vz[i] += w3[i] * mask[i];
#ifdef MHD
		fld.bx[i] += w4[i] * mask[i];
		fld.by[i] += w5[i] * mask[i];
		fld.bz[i] += w6[i] * mask[i];
#endif
	}
	
	// free memory
	fftw_free(x);
	fftw_free(y);
	fftw_free(z);
	
	//done
	return;
}


void init_KidaVortex() {
	const PRECISION a = VORTEX_A;
	const PRECISION b = VORTEX_B;
	
	int i,j,k;
	
	PRECISION w0, x, y;
	PRECISION chi;
	
	chi = b / a;
	w0 = 1.0/chi*(chi + 1)/(chi-1.0);			// According to Kida!
	
	for(i = 0 ; i < NX/NPROC ; i++) {
		x = - LX / 2 + (LX * (i + rank * NX / NPROC)) / NX;
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
	
	for(i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		fld.vx[ i ] +=  I * ky[i] * w1[i] * ik2t[i];
		fld.vy[ i ] += -I * kxt[i] * w1[i] * ik2t[i];
	}
	
	// done
	return;
}

void init_LargeScaleNoise() {
	int i,j,k;
	int num_force=0;
	int total_num_force;
	PRECISION fact;
	
	for( i = 0; i < NX_COMPLEX/NPROC; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			for( k = 0; k < NZ_COMPLEX; k++) {
				if(pow(k2t[ IDX3D ], 0.5) / ( 2.0*M_PI ) < 1.0 / NOISE_CUT_LENGTH) {
					fld.vx[ IDX3D ] += PER_AMPLITUDE_LARGE * mask[IDX3D] * randm() * cexp( I * 2.0*M_PI*randm() ) * NTOTAL;
					fld.vy[ IDX3D ] += PER_AMPLITUDE_LARGE * mask[IDX3D] * randm() * cexp( I * 2.0*M_PI*randm() ) * NTOTAL;
					fld.vz[ IDX3D ] += PER_AMPLITUDE_LARGE * mask[IDX3D] * randm() * cexp( I * 2.0*M_PI*randm() ) * NTOTAL;
#ifdef MHD
					fld.bx[ IDX3D ] += PER_AMPLITUDE_LARGE * mask[IDX3D] * randm() * cexp( I * 2.0*M_PI*randm() ) * NTOTAL;
					fld.by[ IDX3D ] += PER_AMPLITUDE_LARGE * mask[IDX3D] * randm() * cexp( I * 2.0*M_PI*randm() ) * NTOTAL;
					fld.bz[ IDX3D ] += PER_AMPLITUDE_LARGE * mask[IDX3D] * randm() * cexp( I * 2.0*M_PI*randm() ) * NTOTAL;
#endif
					if(mask[IDX3D] > 0) num_force++;
				}
			}
		}
	}
	
	// Get the total number of forced scales.
#ifdef MPI_SUPPORT
	MPI_Allreduce( &num_force, &total_num_force, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#else
	total_num_force=num_force;
#endif
	
	fact=pow(total_num_force,0.5);
	
	// Divide by the total number of modes
	for( i = 0; i < NX_COMPLEX/NPROC; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			for( k = 0; k < NZ_COMPLEX; k++) {
				fld.vx[ IDX3D ] = fld.vx[ IDX3D ] / fact;
				fld.vy[ IDX3D ] = fld.vy[ IDX3D ] / fact;
				fld.vz[ IDX3D ] = fld.vz[ IDX3D ] / fact;
#ifdef MHD
				fld.bx[ IDX3D ] = fld.bx[ IDX3D ] / fact;
				fld.by[ IDX3D ] = fld.by[ IDX3D ] / fact;
				fld.bz[ IDX3D ] = fld.bz[ IDX3D ] / fact;
#endif
			}
		}
	}
	
  symmetrize(fld.vx);
  if(rank==0) fld.vx[0]=0.0;
  symmetrize(fld.vy);
  if(rank==0) fld.vy[0]=0.0;
  symmetrize(fld.vz);
  if(rank==0) fld.vz[0]=0.0;
  
#ifdef MHD
  symmetrize(fld.bx);
  symmetrize(fld.by);
  symmetrize(fld.bz);
#endif
  
}

void init_WhiteNoise() {
	int i,j,k;
	int num_force=0;
	int total_num_force;
	PRECISION fact;
	
	// Excite (2/3)^3*NTOTAL modes
	fact = pow(27.0/8.0*NTOTAL, 0.5);
	
	for( i = 0; i < NX_COMPLEX/NPROC; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			for( k = 0; k < NZ_COMPLEX; k++) {
				fld.vx[ IDX3D ] += PER_AMPLITUDE_NOISE * mask[IDX3D] * randm() * cexp( I * 2.0*M_PI*randm() ) * fact;
				fld.vy[ IDX3D ] += PER_AMPLITUDE_NOISE * mask[IDX3D] * randm() * cexp( I * 2.0*M_PI*randm() ) * fact;
				fld.vz[ IDX3D ] += PER_AMPLITUDE_NOISE * mask[IDX3D] * randm() * cexp( I * 2.0*M_PI*randm() ) * fact;
#ifdef MHD
				fld.bx[ IDX3D ] += PER_AMPLITUDE_NOISE * mask[IDX3D] * randm() * cexp( I * 2.0*M_PI*randm() ) * fact;
				fld.by[ IDX3D ] += PER_AMPLITUDE_NOISE * mask[IDX3D] * randm() * cexp( I * 2.0*M_PI*randm() ) * fact;
				fld.bz[ IDX3D ] += PER_AMPLITUDE_NOISE * mask[IDX3D] * randm() * cexp( I * 2.0*M_PI*randm() ) * fact;
#endif
			}
		}
	}
	
  symmetrize(fld.vx);
  if(rank==0) fld.vx[0]=0.0;
  symmetrize(fld.vy);
  if(rank==0) fld.vy[0]=0.0;
  symmetrize(fld.vz);
  if(rank==0) fld.vz[0]=0.0;
  
#ifdef MHD
  symmetrize(fld.bx);
  symmetrize(fld.by);
  symmetrize(fld.bz);
#endif
  
}

void init_MeanField() {
#ifdef MHD
	if(rank==0) {
		fld.bx[0] = BX0 * ((double) NTOTAL);
		fld.by[0] = BY0 * ((double) NTOTAL);
		fld.bz[0] = BZ0 * ((double) NTOTAL);
	}
#endif
}
/** Init the flow arrays... */	
void init_flow() {
	int i;
	
	DEBUG_START_FUNC;
	// Initialise vectors to 0
	
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		fld.vx[ i ] = 0.0;
		fld.vy[ i ] = 0.0;
		fld.vz[ i ] = 0.0;
#ifdef BOUSSINESQ
		fld.th[ i ] = 0.0;
#endif
#ifdef MHD
		fld.bx[ i ] = 0.0;
		fld.by[ i ] = 0.0;
		fld.bz[ i ] = 0.0;
#endif
	}
	
#ifdef INIT_LARGE_SCALE_NOISE	
	init_LargeScaleNoise();
#endif

#ifdef INIT_VORTEX
	init_KidaVortex();
#endif

#ifdef INIT_SPATIAL_STRUCTURE
	init_SpatialStructure();
#endif

#ifdef INIT_WHITE_NOISE
	init_WhiteNoise();
#endif

#ifdef INIT_MEAN_FIELD
	init_MeanField();
#endif

	projector(fld.vx,fld.vy,fld.vz);
#ifdef MHD
	projector(fld.bx,fld.by,fld.bz);
#endif	

#ifdef DEBUG
	MPI_Printf("Initflow:\n");
	D_show_all(fld);
	MPI_Printf("**************************************************************************************\n");
#endif	
	
	DEBUG_END_FUNC;
	
	return;
}
	
	