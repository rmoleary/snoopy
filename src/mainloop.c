#include <stdlib.h>

#include "common.h"
#include "timestep.h"
#include "output.h"
#include "interface.h"
#include "gfft.h"

#ifdef MPI_SUPPORT
#include "transpose.h"
#endif

#ifdef DEBUG
#include "debug.h"
#endif

struct Field			dfld, fld1;

PRECISION complex		gammaRK[3];
PRECISION complex		xiRK[2];

PRECISION forcing_last_time;

#ifdef WITH_SHEAR
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

PRECISION newdt(PRECISION tremap) {

	int i;
	PRECISION gamma_v;
	PRECISION maxfx   , maxfy, maxfz;
#ifdef MHD
	PRECISION gamma_b;
	PRECISION maxbx   , maxby, maxbz;
#endif
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
	for( i = 0 ; i < NTOTAL_COMPLEX * 2 ; i++) {
		if( fabs( wr1[i] ) > maxfx ) maxfx = fabs( wr1[i] );
		if( fabs( wr2[i] ) > maxfy ) maxfy = fabs( wr2[i] );
		if( fabs( wr3[i] ) > maxfz ) maxfz = fabs( wr3[i] );
		
	}
//}

	maxfx = maxfx / ((double) NTOTAL);
	maxfy = maxfy / ((double) NTOTAL);
	maxfz = maxfz / ((double) NTOTAL);
	
#ifdef MPI_SUPPORT
	reduce(&maxfx,2);
	reduce(&maxfy,2);
	reduce(&maxfz,2);
#endif
	
	gamma_v = (kxmax + fabs(tremap)*kymax) * maxfx + kymax * maxfy + kzmax * maxfz + fabs(OMEGA);
#ifdef WITH_SHEAR
	gamma_v += fabs(SHEAR);
#endif
#ifdef BOUSSINESQ
	gamma_v += pow(fabs(N2), 0.5);
#endif
	
#ifdef MHD
#pragma omp parallel private(i) num_threads ( NTHREADS )
{
		/* Compute the magnetic CFL condition */
	#pragma omp for schedule(static) nowait	
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i] =  fld.bx[i];
		w2[i] =  fld.by[i];
		w3[i] =  fld.bz[i];
	}
}

	gfft_c2r_t(w1);
	gfft_c2r_t(w2);
	gfft_c2r_t(w3);
	
	maxbx=0.0;
	maxby=0.0;
	maxbz=0.0;

	for( i = 0 ; i < NTOTAL_COMPLEX * 2 ; i++) {
		if( fabs( wr1[i] ) > maxbx ) maxbx = fabs( wr1[i] );
		if( fabs( wr2[i] ) > maxby ) maxby = fabs( wr2[i] );
		if( fabs( wr3[i] ) > maxbz ) maxbz = fabs( wr3[i] );
	}

	maxbx = maxbx / ((double) NTOTAL);
	maxby = maxby / ((double) NTOTAL);
	maxbz = maxbz / ((double) NTOTAL);
	
#ifdef MPI_SUPPORT
	reduce(&maxbx,2);
	reduce(&maxby,2);
	reduce(&maxbz,2);
#endif
	
	gamma_b = (kxmax + fabs(tremap)*kymax) * maxbx + kymax * maxby + kzmax * maxbz;
	
	dt = CFL / (gamma_v + gamma_b);
#else
	dt = CFL / gamma_v;
#endif

#ifdef DEBUG
#ifdef MHD
	MPI_Printf("newdt: maxbx=%e, maxby=%e, maxbz=%e\n",maxbx,maxby, maxbz);
#endif
	MPI_Printf("newdt: maxfx=%e, maxfy=%e, maxfz=%e, dt=%e\n",maxfx,maxfy, maxfz, dt);
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
#ifdef MHD
	dfld.bx = (PRECISION complex *) fftw_malloc( sizeof(PRECISION complex) * NTOTAL_COMPLEX);
	if (dfld.bx == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for dfld.bx allocation");
	
	dfld.by = (PRECISION complex *) fftw_malloc( sizeof(PRECISION complex) * NTOTAL_COMPLEX);
	if (dfld.by == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for dfld.by allocation");
	
	dfld.bz = (PRECISION complex *) fftw_malloc( sizeof(PRECISION complex) * NTOTAL_COMPLEX);
	if (dfld.bz == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for dfld.bz allocation");
	
	fld1.bx = (PRECISION complex *) fftw_malloc( sizeof(PRECISION complex) * NTOTAL_COMPLEX);
	if (fld1.bx == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for fld1.bx allocation");
	
	fld1.by = (PRECISION complex *) fftw_malloc( sizeof(PRECISION complex) * NTOTAL_COMPLEX);
	if (fld1.by == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for fld1.by allocation");
	
	fld1.bz = (PRECISION complex *) fftw_malloc( sizeof(PRECISION complex) * NTOTAL_COMPLEX);
	if (fld1.bz == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for fld1.bz allocation");
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
#ifdef MHD
	free(fld1.bx);
	free(fld1.by);
	free(fld1.bz);
	
	free(dfld.bx);
	free(dfld.by);
	free(dfld.bz);

	return;
}


void mainloop() {
	PRECISION		dt = 0.0;
	PRECISION	    t = 0.0;
	PRECISION		tremap = 0.0;
	
	double tstart, tend;
	int i,nloop;
	
	init_mainloop();
	nloop=0;
	
#ifdef RESTART
#ifdef DEBUG
	MPI_Printf("Reading dump file\n");
#endif
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
	
	tstart = get_c_time();
	
	while (t < T_FINAL) {
#ifdef DEBUG
		MPI_Printf("Begining of loop:\n");
		MPI_Printf("fld:\n");
		D_show_all(fld);
		MPI_Printf("**************************************************************************************\n");
#endif

		nloop++;
		if(!(nloop % INTERFACE_CHECK)) check_interface(fld,t,dt,nloop,tstart);
		
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
#ifdef MHD
			fld.bx[i]  = fld.bx[i] + gammaRK[0] * dfld.bx[i] * dt;
			fld.by[i]  = fld.by[i] + gammaRK[0] * dfld.by[i] * dt;
			fld.bz[i]  = fld.bz[i] + gammaRK[0] * dfld.bz[i] * dt;

			fld1.bx[i] = fld.bx[i] + xiRK[0] * dfld.bx[i] * dt;
			fld1.by[i] = fld.by[i] + xiRK[0] * dfld.by[i] * dt;
			fld1.bz[i] = fld.bz[i] + xiRK[0] * dfld.bz[i] * dt;
#endif

			
		}
}		
#ifdef DEBUG
		MPI_Printf("RK, 1st Step:\n");
		MPI_Printf("fld:\n");
		D_show_all(fld);
		MPI_Printf("fld1:\n");
		D_show_all(fld1);
		MPI_Printf("dfld:\n");
		D_show_all(dfld);
		MPI_Printf("**************************************************************************************\n");
#endif
			
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
#ifdef MHD
			fld.bx[i]  = fld1.bx[i] + gammaRK[1] * dfld.bx[i] * dt;
			fld.by[i]  = fld1.by[i] + gammaRK[1] * dfld.by[i] * dt;
			fld.bz[i]  = fld1.bz[i] + gammaRK[1] * dfld.bz[i] * dt;
			
			fld1.bx[i] = fld.bx[i] + xiRK[1] * dfld.bx[i] * dt;
			fld1.by[i] = fld.by[i] + xiRK[1] * dfld.by[i] * dt;
			fld1.bz[i] = fld.bz[i] + xiRK[1] * dfld.bz[i] * dt;
#endif
		}
}
#ifdef DEBUG
		MPI_Printf("RK, 2nd Step:\n");
		MPI_Printf("fld:\n");
		D_show_all(fld);
		MPI_Printf("fld1:\n");
		D_show_all(fld1);
		MPI_Printf("dfld:\n");
		D_show_all(dfld);
		MPI_Printf("**************************************************************************************\n");
#endif

				
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
#ifdef MHD
			fld.bx[i]  = fld1.bx[i] + gammaRK[2] * dfld.bx[i] * dt;
			fld.by[i]  = fld1.by[i] + gammaRK[2] * dfld.by[i] * dt;
			fld.bz[i]  = fld1.bz[i] + gammaRK[2] * dfld.bz[i] * dt;
#endif
		}
}
#ifdef DEBUG
		MPI_Printf("RK, 3rd Step:\n");
		MPI_Printf("fld:\n");
		D_show_all(fld);
		MPI_Printf("fld1:\n");
		D_show_all(fld1);
		MPI_Printf("dfld:\n");
		D_show_all(dfld);
		MPI_Printf("**************************************************************************************\n");
#endif

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
#ifdef MHD
			remap(fld.bx);
			remap(fld.by);
			remap(fld.bz);
		}
		
		kvolve(tremap);
#endif
		
		// Divergence cleaning
		projector(fld.vx,fld.vy,fld.vz);
		projector(fld.bx,fld.by,fld.bz);
				
		output(t);
	}
	tend=get_c_time();
	MPI_Printf("mainloop finished in %d loops and %f seconds (%f sec/loop)\n",nloop,tend-tstart,(tend-tstart)/nloop);
	finish_mainloop();
	return;

}
