#include "common.h"

void init_vortex(PRECISION complex wzf[]) {
	const PRECISION a = 0.04;
	const PRECISION b = 0.16;
	
	int i,j;
	PRECISION w0, x, y;
	PRECISION chi;
	
	chi = b / a;
	w0 = 1.0/chi*(chi + 1)/(chi-1.0);			// According to Kida!
	
	for(i = 0 ; i < NX ; i++) {
		x = - LX / 2 + (LX * i) / NX;
		for(j = 0 ; j < NY ; j++) {
			y = - LY / 2 + (LY * j) / NY;
			if(x * x / (a * a) + y * y / (b * b) < 1) {
				// we are in the vortex
				wr1[j + NY * i] = w0;
			}
			else {
				wr1[j + NY * i] = 0.0;
			}
		}
	}
	
	// transform
	fftw_execute_dft_r2c( r2cfft, wr1, wzf);
	
	// Remove mean vorticity
	wzf[0] = 0.0;
	
	// done
	return;
}

	
void init_flow() {
	int i,j;
	
	// Initialise vectors to 0
	
	for( i = 0; i < NX_COMPLEX; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			fld.wz[ IDX2D ] = 0.0;
#ifdef BOUSSINESQ
			fld.th[ IDX2D ] = 0.0;
#endif
		}
	}
	// Put some noise on the large scales
	
	
	for( i = 0; i < 6; i++) {
		for( j = 0; j < 6; j++) {
			fld.wz[ IDX2D ] = PER_AMPLITUDE * randm() * cexp( I * 2.0 * M_PI * randm() ) * NTOTAL / 5;
		}
	}
	
	init_vortex(fld.wz);
	
	return;
}
	
	