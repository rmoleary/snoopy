#include <stdlib.h>

// Timing support
#ifndef MPI_SUPPORT
#ifndef OPENMP_SUPPORT
#include <stdio.h>
#include <time.h>
#endif
#endif

#include <math.h>
#include <complex.h>
#include <fftw3.h>

#ifdef MPI_SUPPORT
#include <mpi.h>
#ifdef FFTW3_MPI_SUPPORT
#include <fftw3-mpi.h>
#endif
#endif

#ifdef OPENMP_SUPPORT
#include <omp.h>
#endif

#include "gvars.h"
#include "error.h"

// Structures

/**\brief General structure for fields in snoopy. Holds all the information for the field at a given time*/
struct Field {						
	PRECISION complex *vx;			/**< X velocity */
	PRECISION complex *vy;			/**< Y velocity */
	PRECISION complex *vz;			/**< Z velocity */
#ifdef BOUSSINESQ
	PRECISION complex *th;			/**< Entropy perturbations (when using boussinesq) */
#endif
#ifdef MHD
	PRECISION complex *bx;			/**< X magnetic field */
	PRECISION complex *by;			/**< Y magnetic field */
	PRECISION complex *bz;			/**< Z magnetic field */
#endif
};


// This are global variables used throughout the code
// Wave number pointers
PRECISION	*kx;	/**< x Wavevector */
PRECISION	*ky;	/**< y Wavevector */
PRECISION	*kz;	/**< z Wavevector */
PRECISION	*kxt;	/**< Time dependant x Wavevector. Different from kx only when SHEAR is present.*/
PRECISION	*k2t;	/**<  k squared Wavevector, function of time when SHEAR is present.*/
PRECISION	*ik2t;  /**< inverse of k2t Wavevector, function of time when SHEAR is present. set to 0 wheh k2t=0 to avoid singularity*/
PRECISION	kxmax,	kymax,  kzmax,	kmax;	/**< Maximum wavevectors */

fftw_plan	fft_1d_forward, fft_1d_backward;	/**< 1D FFT transforms. Used by remap routines.*/

// Mask for dealiasing
PRECISION   *mask;	/**< Deasliasing Mask*/

PRECISION	*wr1,	*wr2,	*wr3;		/**< Temporary real array (alias of complex w**) */
PRECISION	*wr4,	*wr5,	*wr6;		/**< Temporary real array (alias of complex w**) */
PRECISION	*wr7,	*wr8,	*wr9;		/**< Temporary real array (alias of complex w**) */
PRECISION   *wr10;						/**< Temporary real array (alias of complex w**) */

struct Field			fld;

PRECISION complex		*w1,	*w2,	*w3;	/**< Temporary complex array (alias of real wr**) */
PRECISION complex		*w4,	*w5,	*w6;	/**< Temporary complex array (alias of real wr**) */
PRECISION complex		*w7,	*w8,	*w9;	/**< Temporary complex array (alias of real wr**) */
PRECISION complex		*w10;					/**< Temporary complex array (alias of real wr**) */

PRECISION complex		*w1d, *w2d;

// Physics variables 
PRECISION	nu;
#ifdef BOUSSINESQ
PRECISION	nu_th;
#endif
#ifdef MHD
PRECISION	eta;
#endif

int		rank;

/** Init all global variables, aligning them in memory */
void init_common(void) {
	/* This routine will initialize everything */
	int i,j,k;
#ifdef MPI_SUPPORT
#ifdef FFTW3_MPI_SUPPORT	
	fftw_mpi_init();
#endif
#endif
#ifdef OPEN_SUPPORT	
	if( !(fftw_init_threads()) ) ERROR_HANDLER( ERROR_CRITICAL, "Threads initialisation failed");
#endif
	
	/* We start with the coordinate system */
	kx = (PRECISION *) fftw_malloc( sizeof(PRECISION) * NTOTAL_COMPLEX );
	if (kx == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for kx allocation");
	
	ky = (PRECISION *) fftw_malloc( sizeof(PRECISION) * NTOTAL_COMPLEX );
	if (ky == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for ky allocation");
	
	kz = (PRECISION *) fftw_malloc( sizeof(PRECISION) * NTOTAL_COMPLEX );
	if (kz == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for kz allocation");
	
	kxt = (PRECISION *) fftw_malloc( sizeof(PRECISION) * NTOTAL_COMPLEX );
	if (kxt == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for kxt allocation");
	
	k2t = (PRECISION *) fftw_malloc( sizeof(PRECISION) * NTOTAL_COMPLEX );
	if (k2t == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for k2t allocation");
	
	ik2t = (PRECISION *) fftw_malloc( sizeof(PRECISION) * NTOTAL_COMPLEX );
	if (ik2t == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for ik2t allocation");


	for( i = 0; i < NX_COMPLEX / NPROC; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			for( k = 0; k < NZ_COMPLEX; k++) {
				kx[ IDX3D ] = (2.0 * M_PI) / LX *
						(fmod( NX_COMPLEX * rank / NPROC  + i + (NX_COMPLEX / 2) ,  NX_COMPLEX ) - NX_COMPLEX / 2 );
					 
				ky[ IDX3D ] = (2.0 * M_PI) / LY *
						(fmod( j + (NY_COMPLEX / 2) ,  NY_COMPLEX ) - NY_COMPLEX / 2 );
					 
				kz[ IDX3D ] = (2.0 * M_PI) / LZ * k;

				kxt[ IDX3D ]= kx[IDX3D];
			
				k2t[ IDX3D ] = kxt[IDX3D] * kxt[IDX3D] +
								ky[IDX3D] * ky[IDX3D] +
								kz[IDX3D] * kz[IDX3D];
							  
				if ( k2t[IDX3D] == 0.0 ) ik2t[IDX3D] = 1.0;
				else	ik2t[IDX3D] = 1.0 / k2t[IDX3D];
			}
		}
	}
	
	kxmax = 2.0 * M_PI/ LX * ( (NX / 2) - 1);
	kymax = 2.0 * M_PI/ LY * ( (NY / 2) - 1);
	kzmax = 2.0 * M_PI/ LZ * ( (NZ / 2) - 1);
	
	if( (kxmax>kymax) && (kxmax>kzmax) ) kmax = kxmax;
	else if( (kymax>kxmax) && (kymax > kzmax)) kmax = kymax;
	else kmax=kzmax;
	
	/* Initialize the dealiazing mask Or the niquist frequency mask (in case dealiasing is not required) */
	
	mask = (PRECISION *) fftw_malloc( sizeof(PRECISION) * NTOTAL_COMPLEX );
	if (mask == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for mask allocation");
	
	for( i = 0; i < NX_COMPLEX/NPROC; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			for( k = 0; k < NZ_COMPLEX; k++) {

				mask[ IDX3D ] = 1.0;
#ifdef ANTIALIASING
				if( fabs( kx[ IDX3D ] ) > 2.0/3.0 * kxmax)
					mask[ IDX3D ] = 0.0;
				
				if( fabs( ky[ IDX3D ] ) > 2.0/3.0 * kymax)
					mask[ IDX3D ] = 0.0;
					
				if( fabs( kz[ IDX3D ] ) > 2.0/3.0 * kzmax)
					mask[ IDX3D ] = 0.0;
#else			
			if (  NX_COMPLEX / NPROC * rank + i == NX_COMPLEX / 2 ) 
				mask[ IDX3D ] = 0.0;
			if ( j == NY_COMPLEX / 2 )  
				mask[ IDX3D ] = 0.0;
			if ( k == NZ_COMPLEX ) 
				mask[ IDX3D ] = 0.0;
#endif
			}
		}
	}

#ifdef ANTIALIASING
	kxmax = kxmax * 2.0 / 3.0;
	kymax = kymax * 2.0 / 3.0;
	kzmax = kzmax * 2.0 / 3.0;
#endif

	

// Allocating the fields
// Complex fields
	
	fld.vx = (PRECISION complex *) fftw_malloc( sizeof(PRECISION complex) * NTOTAL_COMPLEX);
	if (fld.vx == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for fld.vx allocation");
	
	fld.vy = (PRECISION complex *) fftw_malloc( sizeof(PRECISION complex) * NTOTAL_COMPLEX);
	if (fld.vy == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for fld.vy allocation");
	
	fld.vz = (PRECISION complex *) fftw_malloc( sizeof(PRECISION complex) * NTOTAL_COMPLEX);
	if (fld.vz == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for fld.vz allocation");
	
#ifdef BOUSSINESQ
	fld.th = (PRECISION complex *) fftw_malloc( sizeof(PRECISION complex) * NTOTAL_COMPLEX);
	if (fld.th == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for fld.th allocation");
#endif
#ifdef MHD
	fld.bx = (PRECISION complex *) fftw_malloc( sizeof(PRECISION complex) * NTOTAL_COMPLEX);
	if (fld.bx == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for fld.vx allocation");
	
	fld.by = (PRECISION complex *) fftw_malloc( sizeof(PRECISION complex) * NTOTAL_COMPLEX);
	if (fld.by == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for fld.vy allocation");
	
	fld.bz = (PRECISION complex *) fftw_malloc( sizeof(PRECISION complex) * NTOTAL_COMPLEX);
	if (fld.bz == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for fld.vz allocation");
#endif

	w1 = (PRECISION complex *) fftw_malloc( sizeof(PRECISION complex) * NTOTAL_COMPLEX);
	if (w1 == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for w1 allocation");
	
	w2 = (PRECISION complex *) fftw_malloc( sizeof(PRECISION complex) * NTOTAL_COMPLEX);
	if (w2 == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for w2 allocation");
	
	w3 = (PRECISION complex *) fftw_malloc( sizeof(PRECISION complex) * NTOTAL_COMPLEX);
	if (w3 == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for w3 allocation");
	
	w4 = (PRECISION complex *) fftw_malloc( sizeof(PRECISION complex) * NTOTAL_COMPLEX);
	if (w4 == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for w4 allocation");
	
	w5 = (PRECISION complex *) fftw_malloc( sizeof(PRECISION complex) * NTOTAL_COMPLEX);
	if (w5 == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for w5 allocation");
	
	w6 = (PRECISION complex *) fftw_malloc( sizeof(PRECISION complex) * NTOTAL_COMPLEX);
	if (w6 == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for w6 allocation");
	
	w7 = (PRECISION complex *) fftw_malloc( sizeof(PRECISION complex) * NTOTAL_COMPLEX);
	if (w7 == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for w7 allocation");
	
	w8 = (PRECISION complex *) fftw_malloc( sizeof(PRECISION complex) * NTOTAL_COMPLEX);
	if (w8 == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for w8 allocation");
	
	w9 = (PRECISION complex *) fftw_malloc( sizeof(PRECISION complex) * NTOTAL_COMPLEX);
	if (w9 == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for w9 allocation");
	
	w10 = (PRECISION complex *) fftw_malloc( sizeof(PRECISION complex) * NTOTAL_COMPLEX);
	if (w10 == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for w10 allocation");
	
	/* Will use the same memory space for real and complex fields */
	
	wr1 = (PRECISION *) w1;
	wr2 = (PRECISION *) w2;
	wr3 = (PRECISION *) w3;
	wr4 = (PRECISION *) w4;
	wr5 = (PRECISION *) w5;
	wr6 = (PRECISION *) w6;
	wr7 = (PRECISION *) w7;
	wr8 = (PRECISION *) w8;
	wr9 = (PRECISION *) w9;
	wr10 = (PRECISION *) w10;

// 1D arrays
	w1d = (PRECISION complex *) fftw_malloc( sizeof(PRECISION complex) * NY);
	if (w1d == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for w1d allocation");
	
	w2d = (PRECISION complex *) fftw_malloc( sizeof(PRECISION complex) * NY);
	if (w2d == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for w2d allocation");

// FFT plans (we use dummy arrays since we use the "guru" interface of fft3 in the code)
// The in place/ out of place will be set automatically at this stage

#ifdef OPENMP_SUPPORT	
	fftw_plan_with_nthreads( 1 );
#endif

	fft_1d_forward = fftw_plan_dft_1d(NY, w1d, w2d, FFTW_FORWARD, FFT_PLANNING);
	if (fft_1d_forward == NULL) ERROR_HANDLER( ERROR_CRITICAL, "FFTW 1D forward plan creation failed");
	
	fft_1d_backward = fftw_plan_dft_1d(NY, w2d, w1d, FFTW_BACKWARD, FFT_PLANNING);
	if (fft_1d_backward == NULL) ERROR_HANDLER( ERROR_CRITICAL, "FFTW 1D backward plan creation failed");
	
// Physic initialisation

	nu = 1.0 / REYNOLDS;
#ifdef BOUSSINESQ	
	nu_th = 1.0 / REYNOLDS_TH;
#endif
#ifdef MHD
	eta = 1.0 / REYNOLDS_M;
#endif
	return;
}

void finish_common(void) {
	free(kx);
	free(ky);
	free(kz);
	free(mask);
	free(fld.vx);
	free(fld.vy);
	free(fld.vz);
#ifdef BOUSSINESQ
	free(fld.th);
#endif
#ifdef MHD
	free(fld.bx);
	free(fld.by);
	free(fld.bz);
#endif
	free(w1);
	free(w2);
	free(w3);
	free(w4);
	free(w5);
	free(w6);
	free(w7);
	free(w8);
	free(w9);
	free(w10);

	return;
}

PRECISION randm(void) {
	PRECISION result;
	result = ( (PRECISION) rand() )/( (double) RAND_MAX );

	return(result);
}

/**
	 * Normal distribution
	 * Algorithm by D.E. Knut, 1997, The Art of Computer Programmin, Addison-Wesley. 
	 */
	 
PRECISION randm_normal(void) {
	PRECISION v1, v2;
	PRECISION rsq=1.0;
	
	while(rsq>=1. || rsq==0.0) {
		v1=2.*randm()-1.0;
		v2=2.*randm()-1.0;
		rsq=v1*v1+v2*v2;
	}
	
	return( v1*sqrt(-2.0 * log(rsq) / rsq));
}

void projector( PRECISION complex qx[],
			    PRECISION complex qy[],
			    PRECISION complex qz[]) {
				
	int i;
	PRECISION complex q0;
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		q0 = kxt[i] * qx[i] + ky[i] * qy[i] + kz[i] * qz[i];
		qx[i] = qx[i] - kxt[i] * q0 * ik2t[i];
		qy[i] = qy[i] - ky[i] * q0 * ik2t[i];
		qz[i] = qz[i] - kz[i] * q0 * ik2t[i];
	}
	return;
}


// Compute the energy of a given field.
// gives the energy in the current process...

PRECISION energy(const PRECISION complex q[]) {
	
	int i,j,k;
	PRECISION energ_tot;
	
	energ_tot=0.0;
	
	for( i = 0; i < NX_COMPLEX/NPROC; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			for( k=0; k < NZ_COMPLEX; k++) {
				if( k == 0) 
					// k=0, we have all the modes.
					energ_tot = energ_tot + creal( 0.5 * q[ IDX3D ] * conj( q[ IDX3D ] ) ) / ((PRECISION) NTOTAL*NTOTAL);
				else
					// k>0, only half of the complex plane is represented.
					energ_tot = energ_tot + creal( q[ IDX3D ] * conj( q[ IDX3D ] ) ) / ((PRECISION) NTOTAL*NTOTAL);
			}
		}
	}
//	energ_tot = 0;
	return(energ_tot);
}
	
double get_c_time(void) {
#ifdef MPI_SUPPORT
	// We have MPI
	return(MPI_Wtime());
#else
#ifdef OPENMP_SUPPORT
	// We don't have MPI, but we have OpenMP
	return(omp_get_wtime());
#else
	// We really have nothing...
	clock_t now;
	now = clock();
	
	return( (double) now / ( (double) CLOCKS_PER_SEC));

#endif
#endif
}
	
	
void reduce(double *var, const int op) {
	// op=1 ADD
	// op=2 Max
	// op=3 Min
	
#ifdef MPI_SUPPORT
	double mpi_temp;
	
	mpi_temp=*var;
	
	if(op==1) MPI_Allreduce( &mpi_temp, var, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	if(op==2) MPI_Allreduce( &mpi_temp, var, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	if(op==3) MPI_Allreduce( &mpi_temp, var, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	
#endif	

	// If no MPI, then this routine does nothing...
	return;
}