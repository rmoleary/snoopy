#include "common.h"
#include "gfft.h"

#ifdef DEBUG
#include "debug.h"
#endif

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
			}
		}
	}
	// Put some noise on the large scales
	/*
	i=NX_COMPLEX-2;
	j=1;
	k=0;
	
	fld.vx[ IDX3D ] = NTOTAL;
	fld.vy[ IDX3D ] = NTOTAL;
	*/
	
	if(rank==0) k0=1;
	else k0=0;
	
	/*
	for( i = 0; i < NX_COMPLEX/NPROC; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			for( k = k0; k < NZ_COMPLEX; k++) {
				fld.vx[ IDX3D ] = PER_AMPLITUDE * mask[IDX3D] * randm() * cexp( I * 2.0*M_PI*randm() );
				fld.vy[ IDX3D ] = PER_AMPLITUDE * mask[IDX3D] * randm() * cexp( I * 2.0*M_PI*randm() );
				fld.vz[ IDX3D ] = PER_AMPLITUDE * mask[IDX3D] * randm() * cexp( I * 2.0*M_PI*randm() );
			}
		}
	}
	*/

	if(rank==0) {
	for( i = 0; i < 4; i++) {
		for( j = 0; j < 4; j++) {
			for( k = k0; k < 4; k++) {
				fld.vx[ IDX3D ] = PER_AMPLITUDE * mask[IDX3D] * randm() * cexp( I * 2.0*M_PI*randm() );
				fld.vy[ IDX3D ] = PER_AMPLITUDE * mask[IDX3D] * randm() * cexp( I * 2.0*M_PI*randm() );
				fld.vz[ IDX3D ] = PER_AMPLITUDE * mask[IDX3D] * randm() * cexp( I * 2.0*M_PI*randm() );
			}
		}
	}
	}

#ifdef DEBUG
	MPI_Printf("Initflow:\n");
	D_show_all(fld);
	MPI_Printf("**************************************************************************************\n");
#endif	
	
	projector(fld.vx,fld.vy,fld.vz);
	
	return;
}
	
	