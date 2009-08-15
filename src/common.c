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

// Timing support
#ifndef MPI_SUPPORT
#ifndef _OPENMP
#include <stdio.h>
#include <time.h>
#endif
#endif

#include "error.h"
#include "gfft.h"
#include "debug.h"

// This are global variables used throughout the code
// Wave number pointers
double	*kx;	/**< x Wavevector */
double	*ky;	/**< y Wavevector */
double	*kz;	/**< z Wavevector */
double	*kxt;	/**< Time dependant x Wavevector. Different from kx only when SHEAR is present.*/
double	*k2t;	/**<  k squared Wavevector, function of time when SHEAR is present.*/
double	*ik2t;  /**< inverse of k2t Wavevector, function of time when SHEAR is present. set to 0 wheh k2t=0 to avoid singularity*/
double	kxmax,	kymax,  kzmax,	kmax;	/**< Maximum wavevectors */


// Mask for dealiasing
double   *mask;	/**< Deasliasing Mask*/

double	*wr1,	*wr2,	*wr3;		/** Temporary real array (alias of complex w**) */
double	*wr4,	*wr5,	*wr6;		/** Temporary real array (alias of complex w**) */
double	*wr7,	*wr8,	*wr9;		/** Temporary real array (alias of complex w**) */

struct Field			fld;

double complex		*w1,	*w2,	*w3;	/**< Temporary complex array (alias of real wr**) */
double complex		*w4,	*w5,	*w6;	/**< Temporary complex array (alias of real wr**) */
double complex		*w7,	*w8,	*w9;	/**< Temporary complex array (alias of real wr**) */

double complex		*pressure;				/**< Pressure field computed during the code evolution.
												 Is only initialized if param.output_pressure is true. */
												

// Parameters
struct Parameters			param;

// Physics variables 
double	nu;
#ifdef BOUSSINESQ
double	nu_th;
#ifdef N2PROFILE
double *N2_profile;
#endif
#endif
#ifdef MHD
double	eta;
#endif

#ifdef MPI_SUPPORT
int		NPROC;									/**< NPROC is a variable when MPI is on. Otherwise, it is preprocessor macro in gvars.h */
#endif

int		rank;
int		nthreads;								/**< Number of OpenMP threads */


/* Function prototypes */
void init_N2_profile();

/** Init all global variables, aligning them in memory */
void init_common(void) {
	/* This routine will initialize everything */
	int i,j,k;
	
	DEBUG_START_FUNC;
	
#ifdef MPI_SUPPORT
#ifdef FFTW3_MPI_SUPPORT	
	fftw_mpi_init();
#endif
#endif
#ifdef _OPENMP
	if( !(fftw_init_threads()) ) ERROR_HANDLER( ERROR_CRITICAL, "Threads initialisation failed");
#endif
	
	/* We start with the coordinate system */
	kx = (double *) fftw_malloc( sizeof(double) * NTOTAL_COMPLEX );
	if (kx == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for kx allocation");
	
	ky = (double *) fftw_malloc( sizeof(double) * NTOTAL_COMPLEX );
	if (ky == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for ky allocation");
	
	kz = (double *) fftw_malloc( sizeof(double) * NTOTAL_COMPLEX );
	if (kz == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for kz allocation");
	
	kxt = (double *) fftw_malloc( sizeof(double) * NTOTAL_COMPLEX );
	if (kxt == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for kxt allocation");
	
	k2t = (double *) fftw_malloc( sizeof(double) * NTOTAL_COMPLEX );
	if (k2t == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for k2t allocation");
	
	ik2t = (double *) fftw_malloc( sizeof(double) * NTOTAL_COMPLEX );
	if (ik2t == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for ik2t allocation");


	for( i = 0; i < NX_COMPLEX / NPROC; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			for( k = 0; k < NZ_COMPLEX; k++) {
				kx[ IDX3D ] = (2.0 * M_PI) / param.lx *
						(fmod( NX_COMPLEX * rank / NPROC  + i + (NX_COMPLEX / 2) ,  NX_COMPLEX ) - NX_COMPLEX / 2 );
					 
				ky[ IDX3D ] = (2.0 * M_PI) / param.ly *
						(fmod( j + (NY_COMPLEX / 2) ,  NY_COMPLEX ) - NY_COMPLEX / 2 );
					 
				kz[ IDX3D ] = (2.0 * M_PI) / param.lz * k;

				kxt[ IDX3D ]= kx[IDX3D];
			
				k2t[ IDX3D ] = kxt[IDX3D] * kxt[IDX3D] +
								ky[IDX3D] * ky[IDX3D] +
								kz[IDX3D] * kz[IDX3D];
							  
				if ( k2t[IDX3D] == 0.0 ) ik2t[IDX3D] = 1.0;
				else	ik2t[IDX3D] = 1.0 / k2t[IDX3D];
			}
		}
	}
	
	kxmax = 2.0 * M_PI/ param.lx * ( (NX / 2) - 1);
	kymax = 2.0 * M_PI/ param.ly * ( (NY / 2) - 1);
	kzmax = 2.0 * M_PI/ param.lz * ( (NZ / 2) - 1);
	
	if( (kxmax>kymax) && (kxmax>kzmax) ) kmax = kxmax;
	else if( (kymax>kxmax) && (kymax > kzmax)) kmax = kymax;
	else kmax=kzmax;
	
	/* Initialize the dealiazing mask Or the nyquist frequency mask (in case dealiasing is not required) */
	
	mask = (double *) fftw_malloc( sizeof(double) * NTOTAL_COMPLEX );
	if (mask == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for mask allocation");
	
	for( i = 0; i < NX_COMPLEX/NPROC; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			for( k = 0; k < NZ_COMPLEX; k++) {

				mask[ IDX3D ] = 1.0;
				if(param.antialiasing) {
					if( fabs( kx[ IDX3D ] ) > 2.0/3.0 * kxmax)
						mask[ IDX3D ] = 0.0;
				
					if( fabs( ky[ IDX3D ] ) > 2.0/3.0 * kymax)
						mask[ IDX3D ] = 0.0;
					
					if( fabs( kz[ IDX3D ] ) > 2.0/3.0 * kzmax)
						mask[ IDX3D ] = 0.0;
				}
				else {			
					if (  NX_COMPLEX / NPROC * rank + i == NX_COMPLEX / 2 ) 
						mask[ IDX3D ] = 0.0;
					if ( j == NY_COMPLEX / 2 )  
						mask[ IDX3D ] = 0.0;
					if ( k == NZ_COMPLEX ) 
						mask[ IDX3D ] = 0.0;
				}
			}
		}
	}

	if(param.antialiasing) {
		kxmax = kxmax * 2.0 / 3.0;
		kymax = kymax * 2.0 / 3.0;
		kzmax = kzmax * 2.0 / 3.0;
		kmax = kmax * 2.0 / 3.0;
	}

	

// Allocate fields
// Complex fields
	
	fld.vx = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (fld.vx == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for fld.vx allocation");
	
	fld.vy = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (fld.vy == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for fld.vy allocation");
	
	fld.vz = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (fld.vz == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for fld.vz allocation");
	
#ifdef BOUSSINESQ
	fld.th = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (fld.th == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for fld.th allocation");
#endif
#ifdef MHD
	fld.bx = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (fld.bx == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for fld.vx allocation");
	
	fld.by = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (fld.by == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for fld.vy allocation");
	
	fld.bz = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (fld.bz == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for fld.vz allocation");
#endif

	w1 = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (w1 == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for w1 allocation");
	
	w2 = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (w2 == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for w2 allocation");
	
	w3 = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (w3 == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for w3 allocation");
	
	w4 = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (w4 == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for w4 allocation");
	
	w5 = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (w5 == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for w5 allocation");
	
	w6 = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (w6 == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for w6 allocation");
	
	w7 = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (w7 == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for w7 allocation");
	
	w8 = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (w8 == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for w8 allocation");
	
	w9 = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (w9 == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for w9 allocation");
	
	if(param.output_pressure) {
		pressure = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
		if (pressure == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for pressure allocation");
	}
	else
		pressure = NULL;
			
	/* Will use the same memory space for real and complex fields */
	
	wr1 = (double *) w1;
	wr2 = (double *) w2;
	wr3 = (double *) w3;
	wr4 = (double *) w4;
	wr5 = (double *) w5;
	wr6 = (double *) w6;
	wr7 = (double *) w7;
	wr8 = (double *) w8;
	wr9 = (double *) w9;

#ifdef BOUSSINESQ
#ifdef N2PROFILE
	N2_profile = (double *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (N2_profile == NULL) ERROR_HANDLER( ERROR_CRITICAL, "no memory for N2 profile allocation");
#endif
#endif
	
// Physic initialisation

#ifdef N2PROFILE
	init_N2_profile();
#endif
	nu = 1.0 / param.reynolds;
#ifdef BOUSSINESQ	
	nu_th = 1.0 / param.reynolds_th;
#endif
#ifdef MHD
	eta = 1.0 / param.reynolds_m;
#endif
	DEBUG_END_FUNC;
	return;
}

void finish_common(void) {
	free(kx);
	free(ky);
	free(kz);
	free(kxt);
	free(k2t);
	free(ik2t);
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

	return;
}
/*********************************************/
/**
Customized random number generator
Allow one to have consistant random numbers
generators on different architectures.
**/
/*********************************************/
double randm(void) {
	const int a	=	16807;
	const int m =	2147483647;
	static int in0 = 13763;
	int q;
	
	// When using mpi, this allows us to have different number series in each process...
	if(in0 == 13763) in0 += 2543 * rank;
	
	/* find random number  */
	q= (int) fmod((double) a * in0, m);
	in0=q;
	
	return((double)q/(double)m);
}
/*********************************************/
/**
	 * Normal distribution
	 * Algorithm by D.E. Knut, 1997, The Art of Computer Programmin, Addison-Wesley. 
	 */
/*********************************************/
	 
double randm_normal(void) {
	double v1, v2;
	double rsq=1.0;
	
	while(rsq>=1. || rsq==0.0) {
		v1=2.*randm()-1.0;
		v2=2.*randm()-1.0;
		rsq=v1*v1+v2*v2;
	}
	
	return( v1*sqrt(-2.0 * log(rsq) / rsq));
}

/****************************************************/
/**
	Remove the divergence from a 3D field using
	the projector operator:
	q=q-k.q/k^2
	
	@param qx: x component of the field
	@param qy: y component of the field
	@param qz: z component of the field
*/
/****************************************************/
	
void projector( double complex qx[],
			    double complex qy[],
			    double complex qz[]) {
				
	int i;
	double complex q0;
	
	DEBUG_START_FUNC;
	
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		q0 = kxt[i] * qx[i] + ky[i] * qy[i] + kz[i] * qz[i];
		qx[i] = qx[i] - kxt[i] * q0 * ik2t[i];
		qy[i] = qy[i] - ky[i] * q0 * ik2t[i];
		qz[i] = qz[i] - kz[i] * q0 * ik2t[i];
	}
	
	DEBUG_END_FUNC;
	
	return;
}

/*********************************************/
/** Compute the energy of a given field.
	@param q complex array containing the field for which
				  we want the total energy
*/
/*********************************************/

double energy(const double complex q[]) {
	
	int i,j,k;
	double energ_tot;
	
	energ_tot=0.0;
	
	for( i = 0; i < NX_COMPLEX/NPROC; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			for( k=0; k < NZ_COMPLEX; k++) {
				if( k == 0) 
					// k=0, we have all the modes.
					energ_tot = energ_tot + creal( 0.5 * q[ IDX3D ] * conj( q[ IDX3D ] ) ) / ((double) NTOTAL*NTOTAL);
				else
					// k>0, only half of the complex plane is represented.
					energ_tot = energ_tot + creal( q[ IDX3D ] * conj( q[ IDX3D ] ) ) / ((double) NTOTAL*NTOTAL);
			}
		}
	}
//	energ_tot = 0;
	return(energ_tot);
}

/********************************************/
/**
Return the localtime in seconds. Use different
implementation depending on the avaiable
libraries  **/
/********************************************/
double get_c_time(void) {
#ifdef MPI_SUPPORT
	// We have MPI
	return(MPI_Wtime());
#else
#ifdef _OPENMP
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

/******************************************/
/**
	Reduce a variable over all the avaiable processes
	Can add a value on all the process, find a maximum
	or a minimum.
	
	NOTE: This routine makes sense only when MPI_SUPPORT
	is set. If not, this routine does nothing.
	
	@param *var: variable to be reduced
	@param op: operation needed to be done. Can be set to:
		1= Sum over all the processes
		2= Find the maximum over all the processes
		3= Find the minimum over all the processes
*/
/*******************************************/
	
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

/*****************************************/
/** Symmetrize the complex space according
**  to the symetries of the real tranform
** @param wi  Field to be symmetrized
*/
/*******************************************/

void symmetrize(double complex wi[]) {
	int i;
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i] = wi[i];
	}
	gfft_c2r(w1);
	gfft_r2c(wr1);
	
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		wi[i] = w1[i] * mask[ i ] / ((double) NTOTAL);
	}
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

	DEBUG_END_FUNC;
	return;
}				

	
	
/***********************************
** Bounadry conditions call
************************************/

void boundary_c(struct Field fldi) {
	DEBUG_START_FUNC;
	
	symmetrize_walls_z(fldi);
	
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

void enforce_symm(struct Field fldi) {
	DEBUG_START_FUNC;
	
	// Enforce symmetries of complex plane
	symmetrize(fldi.vx);
	symmetrize(fldi.vy);
	symmetrize(fldi.vz);
#ifdef BOUSSINESQ
	symmetrize(fldi.th);
#endif
	
	// Remove mean field (noise is generated by the symmetrize funtion)
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

/******************************************
** init the brunt vaissala frequency profile
******************************************/
#ifdef BOUSSINESQ
#ifdef N2PROFILE
void init_N2_profile() {
	double *x,*y,*z;
	int i,j,k;
	
	/*******************************************************************
	** This part does not need to be modified **************************
	********************************************************************/
	// Allocate coordinate arrays
	x = (double *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (x == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for x allocation");
	
	y = (double *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (y == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for y allocation");
	
	z = (double *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (z == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for z allocation");

	// Initialize the (transposed!) arrays
	for(i = 0 ; i < NX ; i++) {
		for(j = 0 ; j < NY/NPROC ; j++) {
			for(k = 0 ; k < NZ ; k++) {
				x[k + (NZ + 2) * i + (NZ + 2) * NX * j] = - param.lx / 2 + (param.lx * i ) / NX;
				y[k + (NZ + 2) * i + (NZ + 2) * NX * j] = - param.ly / 2 + (param.ly * (j + rank * NY / NPROC)) / NY;
				z[k + (NZ + 2) * i + (NZ + 2) * NX * j] = - param.lz / 2 + (param.lz * k ) / NZ;
			}
		}
	}
	
	// Initialize the extra points (k=NZ and k=NZ+1) to zero to prevent stupid things from happening...
	for(i = 0 ; i < NX ; i++) {
		for(j = 0 ; j < NY/NPROC ; j++) {
			for(k = NZ ; k < NZ + 2 ; k++) {
				x[k + (NZ + 2) * i + (NZ + 2) * NX * j] = 0.0;
				y[k + (NZ + 2) * i + (NZ + 2) * NX * j] = 0.0;
				z[k + (NZ + 2) * i + (NZ + 2) * NX * j] = 0.0;
			}
		}
	}
	
	// Init array to zero
	for(i = 0 ; i < NX ; i++) {
		for(j = 0 ; j < NY/NPROC ; j++) {
			for(k = 0 ; k < NZ + 2 ; k++) {
				N2_profile[ k + (NZ + 2) * i + (NZ + 2) * NX * j ] = 0.0;
			}
		}
	}
	
	// Init array
	for(i = 0 ; i < 2*NTOTAL_COMPLEX ; i++) {
		N2_profile[i] = - cos( 2.0 * M_PI * z[i] / param.lz);
	}
	return;
}
#endif
#endif
	