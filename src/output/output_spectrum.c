#include <stdlib.h>

#include "../common.h"
#include "../gfft.h"
#include "../shear.h"
#include "../debug.h"
#include "../particles.h"

#define MAX_N_BIN					10000
#define OUTPUT_SPECTRUM_K_BIN		(2.0 * M_PI)
#define	OUTPUT_SPECTRUM_FILENAME	"spectrum.dat"


/***********************************************************/
/**
	compute a shell-integrated spectrum of the tensor real(wi*wj+)
	and write it on OUTPUT_SPECTRUM_FILENAME
	
	@param wi 1st double complex array from which the spectrum has to be deduced
	@param wj 2nd double complex array from which the spectrum has to be deduced
	@param ti Current time
*/
/***********************************************************/

void write_spectrum(const double complex wi[], const double complex wj[], const double ti) {
	DEBUG_START_FUNC;
	int i,j,k,m;
	int nbin;
	double spectrum[ MAX_N_BIN ];
	FILE *ht;
	
	nbin = (int) ceil(kmax / OUTPUT_SPECTRUM_K_BIN);
	
	for( i = 0; i < MAX_N_BIN; i++ )
		spectrum[ i ] = 0.0;
		
	for( i = 0; i < NX_COMPLEX/NPROC; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			for( k = 0; k < NZ_COMPLEX; k++) {
				m = (int) floor( pow( k2t[ IDX3D ], 0.5 ) / OUTPUT_SPECTRUM_K_BIN + 0.5 );
				if ( m < nbin) {
#ifdef WITH_2D
					if( j == 0)
#else
					if( k == 0) 
#endif
						// k=0, we have all the modes.
						spectrum[ m ] = spectrum[ m ] + creal( wi[ IDX3D ] * conj( wj[ IDX3D ] ) ) / ((double) NTOTAL*NTOTAL);
					else
						// k>0, only half of the complex plane is represented.
						spectrum[ m ] = spectrum[ m ] + creal( 2.0 * wi[ IDX3D ] * conj( wj[ IDX3D ] ) ) / ((double) NTOTAL*NTOTAL);
				}
			}
		}
	}
	
#ifdef MPI_SUPPORT
	// Reduce each component
	for( m=0; m < nbin; m++)
		reduce(&spectrum[m], 1);
#endif

	if(rank==0) {
		ht = fopen(OUTPUT_SPECTRUM_FILENAME,"a");
		fprintf(ht,"%08e\t", ti);
		for( i = 0; i < nbin; i++) 
			fprintf(ht,"%08e\t", spectrum[i]);
	
		fprintf(ht,"\n");
		
		if(ferror(ht)) ERROR_HANDLER( ERROR_CRITICAL, "Error writing spectrum file");
		
		fclose(ht);
	}

		
	DEBUG_END_FUNC;
	return;
}

/**********************************************************/
/**
	Output the transport spectrum in a file (OUTPUT_SPECTRUM_FILENAME)
	This routine is called only when shear is present.
	
	@param fldi: field from which the transport is computed
	@param ti: current time
*/
/*********************************************************/

void output1Dspectrum(const struct Field fldi, const double ti) {
	int i;
	
	DEBUG_START_FUNC;
	
	// The zero array is to be used for dummy spectrums
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i] = 0.0;
	}
	
	// V,B and theta spectrums
	write_spectrum(fldi.vx, fldi.vx, ti);
	write_spectrum(fldi.vy, fldi.vy, ti);
	write_spectrum(fldi.vz, fldi.vz, ti);
	
#ifdef MHD
	write_spectrum(fldi.bx, fldi.bx, ti);
	write_spectrum(fldi.by, fldi.by, ti);
	write_spectrum(fldi.bz, fldi.bz, ti);
#else
#ifdef WITH_LINEAR_TIDE
	write_spectrum(fldi.tvx, fldi.tvx,ti);		// When Linear tide is on, tide statistics replaces B statistics
	write_spectrum(fldi.tvy, fldi.tvy,ti);
	write_spectrum(fldi.tvz, fldi.tvz,ti);
#else
	write_spectrum(w1, w1, ti);
	write_spectrum(w1, w1, ti);
	write_spectrum(w1, w1, ti);
#endif
#endif

#ifdef BOUSSINESQ
	write_spectrum(fldi.th, fldi.th, ti);
#else
	write_spectrum(w1, w1, ti);
#endif

	// Transport spectrums
	write_spectrum(fldi.vx,fldi.vy, ti);
#ifdef MHD
	write_spectrum(fldi.bx,fldi.by, ti);
#else
#ifdef WITH_LINEAR_TIDE
	write_spectrum(fldi.vx,fldi.tvy, ti);
	write_spectrum(fldi.vy,fldi.tvx, ti);
#else
	write_spectrum(w1, w1, ti);
#endif
#endif

	
	// Transfer spectrums
	// Kinetic energy transfer
	
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i] =  fldi.vx[i];
		w2[i] =  fldi.vy[i];
		w3[i] =  fldi.vz[i];
	}

	gfft_c2r_t(w1);
	gfft_c2r_t(w2);
	gfft_c2r_t(w3);
	
		/* Compute the convolution for the advection process */
	for( i = 0 ; i < 2*NTOTAL_COMPLEX ; i++) {
		wr4[i] = wr1[i] * wr1[i] / ((double) NTOTAL*NTOTAL);
		wr5[i] = wr2[i] * wr2[i] / ((double) NTOTAL*NTOTAL);
		wr6[i] = wr3[i] * wr3[i] / ((double) NTOTAL*NTOTAL);
		wr7[i] = wr1[i] * wr2[i] / ((double) NTOTAL*NTOTAL);
		wr8[i] = wr1[i] * wr3[i] / ((double) NTOTAL*NTOTAL);
		wr9[i] = wr2[i] * wr3[i] / ((double) NTOTAL*NTOTAL);
	}
	
	gfft_r2c_t(wr4);
	gfft_r2c_t(wr5);
	gfft_r2c_t(wr6);
	gfft_r2c_t(wr7);
	gfft_r2c_t(wr8);
	gfft_r2c_t(wr9);

	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i] = - I * mask[i] * (
					kxt[i] * w4[i] + ky[i] * w7[i] + kz[i] * w8[i] );
		w2[i] = - I * mask[i] * (
					kxt[i] * w7[i] + ky[i] * w5[i] + kz[i] * w9[i] );
		w3[i] = - I * mask[i] * (
					kxt[i] * w8[i] + ky[i] * w9[i] + kz[i] * w6[i] );
	}
	
	write_spectrum(fldi.vx, w1, ti);
	write_spectrum(fldi.vy, w2, ti);
	write_spectrum(fldi.vz, w3, ti);
	
#ifdef MHD

	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i] =  fldi.vx[i];
		w2[i] =  fldi.vy[i];
		w3[i] =  fldi.vz[i];
		w4[i] =  fldi.bx[i];
		w5[i] =  fldi.by[i];
		w6[i] =  fldi.bz[i];
	}

	gfft_c2r_t(w1);
	gfft_c2r_t(w2);
	gfft_c2r_t(w3);
	gfft_c2r_t(w4);
	gfft_c2r_t(w5);
	gfft_c2r_t(w6);
	
	// (vx,vy,vz) is in w1-w3 and (bx,by,bz) is in (w4-w6). It is now time to compute the emfs in w7-w9...

	for( i = 0 ; i < 2*NTOTAL_COMPLEX ; i++) {
		wr7[i] = (wr2[i] * wr6[i] - wr3[i] * wr5[i]) / ((double) NTOTAL*NTOTAL);
		wr8[i] = (wr3[i] * wr4[i] - wr1[i] * wr6[i]) / ((double) NTOTAL*NTOTAL);
		wr9[i] = (wr1[i] * wr5[i] - wr2[i] * wr4[i]) / ((double) NTOTAL*NTOTAL);
	}

	// Compute the curl of the emf involved in the induction equation.
	
	gfft_r2c_t(wr7);
	gfft_r2c_t(wr8);
	gfft_r2c_t(wr9);
	

	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i] = I * mask[i] * (ky[i] * w9[i] - kz[i] * w8[i]);
		w2[i] = I * mask[i] * (kz[i] * w7[i] - kxt[i]* w9[i]);
		w3[i] = I * mask[i] * (kxt[i]* w8[i] - ky[i] * w7[i]);
	}

	write_spectrum(fldi.bx, w1, ti);
	write_spectrum(fldi.by, w2, ti);
	write_spectrum(fldi.bz, w3, ti);
	
	// Let's do the Lorentz Force
	// We already have (bx,by,bz) in w4-w6. No need to compute them again...

	for( i = 0 ; i < 2*NTOTAL_COMPLEX ; i++) {
		wr1[i] = wr4[i] * wr4[i] / ((double) NTOTAL*NTOTAL);
		wr2[i] = wr5[i] * wr5[i] / ((double) NTOTAL*NTOTAL);
		wr3[i] = wr6[i] * wr6[i] / ((double) NTOTAL*NTOTAL);
		wr7[i] = wr4[i] * wr5[i] / ((double) NTOTAL*NTOTAL);
		wr8[i] = wr4[i] * wr6[i] / ((double) NTOTAL*NTOTAL);
		wr9[i] = wr5[i] * wr6[i] / ((double) NTOTAL*NTOTAL);
	}

	gfft_r2c_t(wr1);
	gfft_r2c_t(wr2);
	gfft_r2c_t(wr3);
	gfft_r2c_t(wr7);
	gfft_r2c_t(wr8);
	gfft_r2c_t(wr9);


	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w4[i] = I * mask[i] * (kxt[i] * w1[i] + ky[i] * w7[i] + kz[i] * w8[i]);
		w5[i] = I * mask[i] * (kxt[i] * w7[i] + ky[i] * w2[i] + kz[i] * w9[i]);
		w6[i] = I * mask[i] * (kxt[i] * w8[i] + ky[i] * w9[i] + kz[i] * w3[i]);
	}
	
	write_spectrum(fldi.vx, w4, ti);
	write_spectrum(fldi.vy, w5, ti);
	write_spectrum(fldi.vz, w6, ti);
	
	// Helicity spectrums
	
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i] = I * ik2t[i] * (ky[i] * fldi.bz[i] - kz[i] * fldi.by[i] );
		w2[i] = I * ik2t[i] * (kz[i] * fldi.bx[i] - kxt[i]* fldi.bz[i] );
		w3[i] = I * ik2t[i] * (kxt[i]* fldi.by[i] - ky[i] * fldi.bx[i] );
	}
	
	write_spectrum(fldi.bx, w1, ti);
	write_spectrum(fldi.by, w2, ti);
	write_spectrum(fldi.bz, w3, ti);
	
#else
	// The zero array is to be used for dummy spectrums
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i] = 0.0;
	}
	
	write_spectrum(w1, w1, ti);
	write_spectrum(w1, w1, ti);
	write_spectrum(w1, w1, ti);
	write_spectrum(w1, w1, ti);
	write_spectrum(w1, w1, ti);
	write_spectrum(w1, w1, ti);
	write_spectrum(w1, w1, ti);
	write_spectrum(w1, w1, ti);
	write_spectrum(w1, w1, ti);
#endif
	
	DEBUG_END_FUNC;
	
	return;
}

/**********************************************************/
/**
	Initialise the 1D spectrum output routine, used
	to output the spectrum
	This routine print the mode ks in the first line
	It also counts the number of mode in each shell and 
	output it in the second line of OUTPUT_SPECTRUM_FILENAME
*/
/*********************************************************/
void init1Dspectrum() {
	int i,j,k,m;
	int nbin;
	FILE * ht;
	double spectrum[ MAX_N_BIN ];
	
	DEBUG_START_FUNC;
	
	nbin = (int) ceil(kmax / OUTPUT_SPECTRUM_K_BIN);
	
	if(rank==0) {
		ht = fopen(OUTPUT_SPECTRUM_FILENAME,"w");
		
		for( m=0; m < nbin; m++) 
			fprintf(ht,"%08e\t", m * OUTPUT_SPECTRUM_K_BIN);
	
		fprintf(ht,"\n");
	}
	
	for( i = 0; i < MAX_N_BIN ; i++ )
		spectrum[ i ] = 0.0;
		
	for( i = 0; i < NX_COMPLEX/NPROC ; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			for( k = 0; k < NZ_COMPLEX; k++) {
				m = (int) floor( pow( k2t[ IDX3D ], 0.5 ) / OUTPUT_SPECTRUM_K_BIN + 0.5 );
				if ( m < nbin)
					spectrum[ m ] = spectrum[ m ] + 1.0;
			}
		}
	}
	
#ifdef MPI_SUPPORT
	// Reduce each component
	for( m=0; m < nbin ; m++)
		reduce(&spectrum[m], 1);
#endif

	if(rank==0) {
		for( i = 0; i < nbin ; i++) 
			fprintf(ht,"%08e\t", spectrum[i]);
	
		fprintf(ht,"\n");
	
		if(ferror(ht)) ERROR_HANDLER( ERROR_CRITICAL, "Error writing spectrum file");
		fclose(ht);
	}
	
	DEBUG_END_FUNC;
	
	return;
}
