#include "common.h"

void init_vortex(PRECISION complex wzf[]) {
	const PRECISION a = 0.04;
	const PRECISION b = 0.16;
	
	int i,j,k;
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
	fftw_execute_dft_r2c( r2cfft, wr1, w1);
	
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
	int i,j,k;
	
	// Initialise vectors to 0
	
	for( i = 0; i < NX_COMPLEX; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			for( k = 0; k < NY_COMPLEX; k++) {
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
	i=3;
	j=1;
	k=0;
	
	fld.vx[ IDX3D ] = NTOTAL;
	fld.vy[ IDX3D ] = NTOTAL;
	
	projector(fld.vx,fld.vy,fld.vz);
	
	return;
}
	
	