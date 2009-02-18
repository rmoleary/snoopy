#include <stdlib.h>

#include "common.h"
#include "timestep.h"
#include "output.h"
#include "interface.h"
#include "gfft.h"

struct Field			dfld, fld1;

PRECISION complex		gammaRK[3];
PRECISION complex		xiRK[2];

PRECISION forcing_last_time;

#ifdef WITH_SHEAR
void remap(PRECISION complex qi[]) {
	int i, j, k;
	int nx, ny, nxtarget;
	
#ifdef DEBUG
	printf("remap called\n");
#endif
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i]=0.0;
	}

	for( i = 0; i < NX_COMPLEX; i++) {
		nx = fmod( i + (NX_COMPLEX / 2) ,  NX_COMPLEX ) - NX_COMPLEX / 2 ;
		for( j = 0; j < NY_COMPLEX; j++) {
			ny = fmod( j + (NY_COMPLEX / 2) ,  NY_COMPLEX ) - NY_COMPLEX / 2 ;
			
			nxtarget = nx + ny;		// We have a negative shear, hence nx plus ny
			
			if ( nxtarget <= - NX_COMPLEX / 2 ) break;
			if ( nxtarget >=   NX_COMPLEX / 2 ) break;
			
			if ( nxtarget <0 ) nxtarget = nxtarget + NX_COMPLEX;
			
			for( k = 0; k < NZ_COMPLEX; k++) {
				w1[k + NZ_COMPLEX * ny + NZ_COMPLEX * NY_COMPLEX * nxtarget] = qi[ IDX3D ];
			}
		}
	}
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		qi[i] = w1[i] * mask[i];
	}
	return;
}

void kvolve(const PRECISION tremap) {
	int i, j, k;
#pragma omp parallel private(i,j) num_threads ( NTHREADS )
{
		/* Compute the convolution */
	#pragma omp for schedule(static) nowait	
	for( i = 0; i < NX_COMPLEX; i++) {
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

PRECISION newdt(PRECISION tremap) {

	int i;
	PRECISION gamma_v;
	PRECISION maxfx   , maxfy, maxfz;
	PRECISION dt;
	
#pragma omp parallel private(i) num_threads ( NTHREADS )
{
		/* Compute the convolution */
	#pragma omp for schedule(static) nowait	
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i] =  fld.vx[i];
		w2[i] =  fld.vy[i];
		w3[i] =  fld.vz[i];
	}
}

	gfft_c2r_t(w1);
	gfft_c2r_t(w2);
	gfft_c2r_t(w3);
	
	maxfx=0.0;
	maxfy=0.0;
	maxfz=0.0;

//#pragma omp parallel private(i) num_threads ( NTHREADS )
//{
		/* Compute the convolution */
//	#pragma omp for schedule(static) nowait		
	for( i = 0 ; i < NTOTAL ; i++) {
		if( fabs( wr1[i] ) > maxfx ) maxfx = fabs( wr1[i] );
		if( fabs( wr2[i] ) > maxfy ) maxfy = fabs( wr2[i] );
		if( fabs( wr3[i] ) > maxfz ) maxfz = fabs( wr3[i] );
		
	}
//}

	maxfx = maxfx / ((double) NTOTAL);
	maxfy = maxfy / ((double) NTOTAL);
	maxfz = maxfz / ((double) NTOTAL);
	
	gamma_v = (kxmax + fabs(tremap)*kymax) * maxfx + kymax * maxfy + kzmax * maxfz + fabs(OMEGA);
#ifdef WITH_SHEAR
	gamma_v += fabs(SHEAR);
#endif
#ifdef BOUSSINESQ
	gamma_v += pow(fabs(N2), 0.5);
#endif
	
	dt = CFL / gamma_v;

#ifdef DEBUG
	printf("newdt: maxfx=%e, maxfy=%e, dt=%e\n",maxfx,maxfy,dt);
#endif

	return(dt);
}			   			   
		
void init_mainloop() {
	dfld.vx = (PRECISION complex *) fftw_malloc( sizeof(PRECISION complex) * NTOTAL_COMPLEX);
	if (dfld.vx == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for dfld.vx allocation");
	
	dfld.vy = (PRECISION complex *) fftw_malloc( sizeof(PRECISION complex) * NTOTAL_COMPLEX);
	if (dfld.vy == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for dfld.vy allocation");
	
	dfld.vz = (PRECISION complex *) fftw_malloc( sizeof(PRECISION complex) * NTOTAL_COMPLEX);
	if (dfld.vz == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for dfld.vz allocation");
	
	fld1.vx = (PRECISION complex *) fftw_malloc( sizeof(PRECISION complex) * NTOTAL_COMPLEX);
	if (fld1.vx == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for fld1.vx allocation");
	
	fld1.vy = (PRECISION complex *) fftw_malloc( sizeof(PRECISION complex) * NTOTAL_COMPLEX);
	if (fld1.vy == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for fld1.vy allocation");
	
	fld1.vz = (PRECISION complex *) fftw_malloc( sizeof(PRECISION complex) * NTOTAL_COMPLEX);
	if (fld1.vz == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for fld1.vz allocation");

#ifdef BOUSSINESQ
	dfld.th = (PRECISION complex *) fftw_malloc( sizeof(PRECISION complex) * NTOTAL_COMPLEX);
	if (dfld.th == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for dfld.th allocation");
	
	fld1.th = (PRECISION complex *) fftw_malloc( sizeof(PRECISION complex) * NTOTAL_COMPLEX);
	if (fld1.th == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for fld1.th allocation");
#endif
	
// Init the Runge-Kutta timestepping

	gammaRK[0] = 8.0 / 15.0;
	gammaRK[1] = 5.0 / 12.0;
	gammaRK[2] = 3.0 / 4.0;
	
	xiRK[0] = -17.0 / 60.0;
	xiRK[1] = -5.0 / 12.0;

/*
	gammaRK[0] = 1.0;
	gammaRK[1] = 0.0;
	gammaRK[2] = 0.0;
	
	xiRK[0] = 0.0;
	xiRK[1] = 0.0;
*/
	
	return;
}

void finish_mainloop() {
	free(fld1.vx);
	free(fld1.vy);
	free(fld1.vz);
	
	free(dfld.vx);
	free(dfld.vy);
	free(dfld.vz);
	
#ifdef BOUSSINESQ
	free(fld1.th);
	free(dfld.th);
#endif
	return;
}


void mainloop() {
	PRECISION		dt = 0.0;
	PRECISION	    t = 0.0;
	PRECISION		tremap = 0.0;
	int i,nloop;
	
	init_mainloop();
	nloop=0;
	
#ifdef RESTART
	read_dump(fld,&t);
#else
	t = T_INITIAL;
	clear_timevar();
	output(t);
#endif

#ifdef WITH_SHEAR	
	tremap = fmod(t + LY / (2.0 * SHEAR * LX) , LY / (SHEAR * LX)) - LY / (2.0 * SHEAR * LX);
	kvolve(tremap);
#else
	tremap = 0.0;
#endif
	
	while (t < T_FINAL) {
		nloop++;
		if(!(nloop % INTERFACE_CHECK)) check_interface(fld,t,dt,nloop);
		
		dt = newdt(tremap);
		
		// This is an order 3 runge Kutta scheme with low storage
		
		// 1st RK3 step
		
		timestep(dfld, fld, t, dt);
#pragma omp parallel private(i) num_threads ( NTHREADS )
{
	#pragma omp for schedule(static) nowait	
		for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
			fld.vx[i]  = fld.vx[i] + gammaRK[0] * dfld.vx[i] * dt;
			fld.vy[i]  = fld.vy[i] + gammaRK[0] * dfld.vy[i] * dt;
			fld.vz[i]  = fld.vz[i] + gammaRK[0] * dfld.vz[i] * dt;

			fld1.vx[i] = fld.vx[i] + xiRK[0] * dfld.vx[i] * dt;
			fld1.vy[i] = fld.vy[i] + xiRK[0] * dfld.vy[i] * dt;
			fld1.vz[i] = fld.vz[i] + xiRK[0] * dfld.vz[i] * dt;
			
#ifdef BOUSSINESQ
			fld.th[i]  = fld.th[i] + gammaRK[0] * dfld.th[i] * dt;

			fld1.th[i] = fld.th[i] + xiRK[0] * dfld.th[i] * dt;
#endif
			
		}
}		
		// 2nd RK3 step
#ifdef WITH_SHEAR
		kvolve(tremap+gammaRK[0]*dt);
#endif
		
		timestep(dfld, fld, t, dt);

#pragma omp parallel private(i) num_threads ( NTHREADS )
{
	#pragma omp for schedule(static) nowait		
		for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
			fld.vx[i]  = fld1.vx[i] + gammaRK[1] * dfld.vx[i] * dt;
			fld.vy[i]  = fld1.vy[i] + gammaRK[1] * dfld.vy[i] * dt;
			fld.vz[i]  = fld1.vz[i] + gammaRK[1] * dfld.vz[i] * dt;
			
			fld1.vx[i] = fld.vx[i] + xiRK[1] * dfld.vx[i] * dt;
			fld1.vy[i] = fld.vy[i] + xiRK[1] * dfld.vy[i] * dt;
			fld1.vz[i] = fld.vz[i] + xiRK[1] * dfld.vz[i] * dt;
			
#ifdef BOUSSINESQ
			fld.th[i]  = fld1.th[i] + gammaRK[1] * dfld.th[i] * dt;
			fld1.th[i] = fld.th[i] + xiRK[1] * dfld.th[i] * dt;
#endif
		}
}
				
		// 3rd RK3 Step
#ifdef WITH_SHEAR
		kvolve(tremap + (gammaRK[0] + xiRK[0] + gammaRK[1]) * dt );
#endif
		
		timestep(dfld, fld, t, dt);

#pragma omp parallel private(i) num_threads ( NTHREADS )
{
	#pragma omp for schedule(static) nowait			
		for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
			fld.vx[i]  = fld1.vx[i] + gammaRK[2] * dfld.vx[i] * dt;
			fld.vy[i]  = fld1.vy[i] + gammaRK[2] * dfld.vy[i] * dt;
			fld.vz[i]  = fld1.vz[i] + gammaRK[2] * dfld.vz[i] * dt;
#ifdef BOUSSINESQ
			fld.th[i]  = fld1.th[i] + gammaRK[2] * dfld.th[i] * dt;
#endif
		}
}
		// Runge Kutta finished
		
		// Implicit step
		implicitstep(fld, t, dt);
		
		// evolving the frame
		t = t + dt;

#ifdef WITH_SHEAR		
		tremap = tremap + dt;
		
		if(tremap > LY / (2.0 * SHEAR * LX)) {
			tremap = fmod(t + LY / (2.0 * SHEAR * LX) , LY /  (LX * SHEAR)) - LY / (2.0 * SHEAR * LX);
			remap(fld.vx);
			remap(fld.vy);
			remap(fld.vz);
#ifdef BOUSSINESQ
			remap(fld.th);
#endif
		}
		
		kvolve(tremap);
#endif
		
		// Divergence cleaning
		projector(fld.vx,fld.vy,fld.vz);
				
		output(t);
	}
	printf("mainloop finished\n");
	finish_mainloop();
	return;

}
