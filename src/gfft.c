#include "gvars.h"
#include "common.h"

#ifdef MPI_SUPPORT
// That's a long MPI if def...

#include "transpose.h"

#include <pthread.h>

const int n_size2D[2] = {NX, NZ};
const int n_size1D[1] = {NY_COMPLEX};

int use_thread;

fftw_plan	r2c_2d, c2r_2d, r2c_1d, c2r_1d;

PRECISION complex *wi1, *wi2, *wi3;
PRECISION *wir1, *wir2, *wir3;

double time_1D, time_2D, time_transpose;

// Used by the transposition thread
pthread_t thread[1];
pthread_attr_t thread_attr;

/* GFFT (Like Geo's FFT as you might have guessed...) is an FFT wrapper for FFTW>=3.2
It takes care of the MPI part of the FFT while FFTW deals with the FFT themselves
GFFT can handle threaded FFTS (just ask...) and communication overlaps while doing FFT (requires POSIX thread support).
One concludes It's FAAAARRRR better than FFTW 2.1.5
*/

	
// This is an inplace real 2 complex transform
// Assumes wrin has the logical dimensions [NY/PROC, NX, NZ] of real positions
// physical dimensions [NY/NPROC, NX, NZ+2];
void gfft_r2c_t(PRECISION *wrin) {
	int i;
	PRECISION complex *win = (PRECISION complex *) wrin;
	
	//start transforming in 2D wrin
	time_2D=time_2D-get_walltime();
	fftw_execute_dft_r2c(r2c_2d, wrin, win);
	time_2D=time_2D+get_walltime();
	
	// The logical dimensions of win are [NY_COMPLEX/NPROC, NX_COMPLEX, NZ_COMPLEX]
	// transpose it
	time_transpose=time_transpose-get_walltime();
	transpose_complex_YX(win, win);
	time_transpose=time_transpose+get_walltime();
	
	// We now have an array with logical dimensions[NX_COMPLEX/NPROC, NY_COMPLEX, NZ_COMPLEX]
	time_1D=time_1D-get_walltime();
#pragma omp parallel private(i) num_threads ( NTHREADS )
{
	#pragma omp for schedule(static) nowait	
	for(i=0 ; i < NX_COMPLEX/NPROC ; i++) 
		fftw_execute_dft(r2c_1d, &win[i*NY_COMPLEX*NZ_COMPLEX],&win[i*NY_COMPLEX*NZ_COMPLEX]);
}
	time_1D=time_1D+get_walltime();
	// done...
	return;
}


void gfft_c2r_t(PRECISION complex *win) {
	int i;
	PRECISION *wrin = (PRECISION *) win;
	// We now have an array with logical dimensions[NX_COMPLEX/NPROC, NY_COMPLEX, NZ_COMPLEX]
	time_1D=time_1D-get_walltime();
#pragma omp parallel private(i) num_threads ( NTHREADS )
{
	#pragma omp for schedule(static) nowait	
	for(i=0 ; i < NX_COMPLEX/NPROC ; i++) 
		fftw_execute_dft(c2r_1d, &win[i*NY_COMPLEX*NZ_COMPLEX],&win[i*NY_COMPLEX*NZ_COMPLEX]);
}		
	time_1D=time_1D+get_walltime();
	// The logical dimensions of win are [NX_COMPLEX/NPROC, NY_COMPLEX, NZ_COMPLEX]
	// transpose it
	time_transpose=time_transpose-get_walltime();
	transpose_complex_XY(win, win);
	time_transpose=time_transpose+get_walltime();
	
	// The final 2D transform
	time_2D=time_2D-get_walltime();
	fftw_execute_dft_c2r(c2r_2d, win, wrin);
	time_2D=time_2D+get_walltime();
	 // and we're done !
	 return;
}
														  
	
// In place double transpose transforms
// Not Fast, but convenient...!

void gfft_r2c(PRECISION *wrin) {
	transpose_real(NX, NY, NZ+2, NPROC, wrin, wrin);
	gfft_r2c_t(wrin);
	return;
}

void gfft_c2r(PRECISION complex *win) {
	PRECISION *wrin = (PRECISION *) win;
	gfft_c2r_t(win);
	transpose_real(NY,NX,NZ+2,NPROC,wrin,wrin);
	return;
}


// This are doing 3FFTs at a time. Uses communication overlapping...
void gfft3_r2c_t(PRECISION *wrin1, PRECISION *wrin2, PRECISION *wrin3) {
	int i;	
	void *status;
	int rc;

	PRECISION complex *win1 = (PRECISION complex *) wrin1;
	PRECISION complex *win2 = (PRECISION complex *) wrin2;
	PRECISION complex *win3 = (PRECISION complex *) wrin3;
	
	
	//Begining of first transform
	time_2D=time_2D-get_walltime();
	fftw_execute_dft_r2c(r2c_2d, wrin1, win1);
	time_2D=time_2D+get_walltime();

	if(use_thread) 
		rc = pthread_create(thread, &thread_attr, transpose_complex_YX_thread, (void*) win1);
	else
		transpose_complex_YX(win1, win1);
	
	// Begining of second transform
	time_2D=time_2D-get_walltime();
	fftw_execute_dft_r2c(r2c_2d, wrin2, win2);
	time_2D=time_2D+get_walltime();

	time_transpose=time_transpose-get_walltime();
    if(use_thread) {	
		// Wait for the previous transform
		rc = pthread_join(*thread, &status);
		// Here we should wait for the first transpose to finish...
		rc = pthread_create(thread, &thread_attr, transpose_complex_YX_thread, (void*) win2);
	}
	else
		transpose_complex_YX(win2, win2);

	time_transpose=time_transpose+get_walltime();
	
	// En of first, begining of third...
	time_1D=time_1D-get_walltime();
#pragma omp parallel private(i) num_threads ( NTHREADS )
{
	#pragma omp for schedule(static) nowait	
	for(i=0 ; i < NX_COMPLEX/NPROC ; i++) 
		fftw_execute_dft(r2c_1d, &win1[i*NY_COMPLEX*NZ_COMPLEX],&win1[i*NY_COMPLEX*NZ_COMPLEX]);
}
	time_1D=time_1D+get_walltime();
	
	time_2D=time_2D-get_walltime();
	fftw_execute_dft_r2c(r2c_2d, wrin3, win3);
	time_2D=time_2D+get_walltime();

	time_transpose=time_transpose-get_walltime();
    if(use_thread) {
		//TRANSPOSE_WAIT
		// Wait for the previous transform
		rc = pthread_join(*thread, &status);
		// Here we should wait for the first transpose to finish...
		rc = pthread_create(thread, &thread_attr, transpose_complex_YX_thread, (void*) win3);
	}
	else
		transpose_complex_YX(win3, win3);

	time_transpose=time_transpose+get_walltime();
	
	// End of second
	time_1D=time_1D-get_walltime();
#pragma omp parallel private(i) num_threads ( NTHREADS )
{
	#pragma omp for schedule(static) nowait	
	for(i=0 ; i < NX_COMPLEX/NPROC ; i++) 
		fftw_execute_dft(r2c_1d, &win2[i*NY_COMPLEX*NZ_COMPLEX],&win2[i*NY_COMPLEX*NZ_COMPLEX]);
}
	time_1D=time_1D+get_walltime();
	
	if(use_thread) {
		time_transpose=time_transpose-get_walltime();
		// TRANSPOSE WAIT
		// Wait for the previous transform
		rc = pthread_join(*thread, &status);
		time_transpose=time_transpose+get_walltime();
	}

	time_1D=time_1D-get_walltime();
#pragma omp parallel private(i) num_threads ( NTHREADS )
{
	#pragma omp for schedule(static) nowait	
	for(i=0 ; i < NX_COMPLEX/NPROC ; i++) 
		fftw_execute_dft(r2c_1d, &win3[i*NY_COMPLEX*NZ_COMPLEX],&win3[i*NY_COMPLEX*NZ_COMPLEX]);
}
	time_1D=time_1D+get_walltime();
	// done...
	return;
}

void gfft3_c2r_t(PRECISION complex *win1, PRECISION complex *win2, PRECISION complex *win3) {
	int i;	
	void *status;
	int rc;
	
	PRECISION *wrin1 = (PRECISION *) win1;
	PRECISION *wrin2 = (PRECISION *) win2;
	PRECISION *wrin3 = (PRECISION *) win3;
	
	// Begining of first transform
	time_1D=time_1D-get_walltime();
#pragma omp parallel private(i) num_threads ( NTHREADS )
{
	#pragma omp for schedule(static) nowait
	for(i=0 ; i < NX_COMPLEX/NPROC ; i++) 
		fftw_execute_dft(c2r_1d, &win1[i*NY_COMPLEX*NZ_COMPLEX],&win1[i*NY_COMPLEX*NZ_COMPLEX]);
}
	time_1D=time_1D+get_walltime();

	time_transpose=time_transpose-get_walltime();
	if(use_thread)
		rc = pthread_create(thread, &thread_attr, transpose_complex_XY_thread, (void*) win1);
	else
		transpose_complex_XY(win1, win1);
	time_transpose=time_transpose+get_walltime();

	// Begining of second transform
	time_1D=time_1D-get_walltime();
#pragma omp parallel private(i) num_threads ( NTHREADS )
{
	#pragma omp for schedule(static) nowait
	for(i=0 ; i < NX_COMPLEX/NPROC ; i++) 
		fftw_execute_dft(c2r_1d, &win2[i*NY_COMPLEX*NZ_COMPLEX],&win2[i*NY_COMPLEX*NZ_COMPLEX]);
	
}
	time_1D=time_1D+get_walltime();

	time_transpose=time_transpose-get_walltime();
	if(use_thread) {
		// Wait for the previous transform
		rc = pthread_join(*thread, &status);
		// Here we should wait for the first transpose to finish...
		rc = pthread_create(thread, &thread_attr, transpose_complex_XY_thread, (void*) win2);
    }
	else 
		transpose_complex_XY(win2, win2);
	time_transpose=time_transpose + get_walltime();
	
	// End of first, begining of third
	time_2D=time_2D-get_walltime();
	fftw_execute_dft_c2r(c2r_2d, win1, wrin1);
	time_2D=time_2D+get_walltime();
	
	time_1D=time_1D-get_walltime();
#pragma omp parallel private(i) num_threads ( NTHREADS )
{
	#pragma omp for schedule(static) nowait
	for(i=0 ; i < NX_COMPLEX/NPROC ; i++) 
		fftw_execute_dft(c2r_1d, &win3[i*NY_COMPLEX*NZ_COMPLEX],&win3[i*NY_COMPLEX*NZ_COMPLEX]);
}
	time_1D=time_1D+get_walltime();

	// TRANSPOSE_WAIT
	time_transpose=time_transpose-get_walltime();
	if(use_thread) {
		// Wait for the previous transform
		rc = pthread_join(*thread, &status);
		// Here we should wait for the first transpose to finish...
		rc = pthread_create(thread, &thread_attr, transpose_complex_XY_thread, (void*) win3);
	}
	else
		transpose_complex_XY(win3, win3);
	time_transpose=time_transpose+get_walltime();
	
	// End of second
	time_2D=time_2D-get_walltime();
	fftw_execute_dft_c2r(c2r_2d, win2, wrin2);
	time_2D=time_2D+get_walltime();
	
    if(use_thread) {
		time_transpose=time_transpose-get_walltime();
		rc = pthread_join(*thread, &status);
		time_transpose=time_transpose+get_walltime();
    }

	// End of third
	time_2D=time_2D-get_walltime();
	fftw_execute_dft_c2r(c2r_2d, win3, wrin3);
	time_2D=time_2D+get_walltime();
	
	 // and we're done !
	 return;
}
	

void init_gfft() {
	// This will init the plans needed by gfft
	// Transform of NY/NPROC arrays of (logical) size [NX, NZ]
	// The physical size is [NX, NZ+2]
	// We use in-place transforms
	int i;
#ifndef FORCE_NO_THREADS
	PRECISION timeinit1, timeinit2;
#endif
	
	wi1 = (PRECISION complex *) fftw_malloc( sizeof(PRECISION complex) * NTOTAL_COMPLEX);
	if (wi1 == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for wi1 allocation");
	wi2 = (PRECISION complex *) fftw_malloc( sizeof(PRECISION complex) * NTOTAL_COMPLEX);
	if (wi2 == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for wi2 allocation");
	wi3 = (PRECISION complex *) fftw_malloc( sizeof(PRECISION complex) * NTOTAL_COMPLEX);
	if (wi3 == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for wi3 allocation");

	wir1 = (PRECISION *) wi1;
	wir2 = (PRECISION *) wi2;
	wir3 = (PRECISION *) wi3;
	
	for(i = 0 ; i < NTOTAL_COMPLEX; i++) {
		wi1[i]=1.0;
		wi2[i]=1.0;
		wi3[i]=1.0;
	}
	
#ifdef OPENMP_SUPPORT
	fftw_plan_with_nthreads( NTHREADS );
#endif
	r2c_2d = fftw_plan_many_dft_r2c(2, n_size2D, NY / NPROC, wir1, NULL, 1, (NZ+2)*NX,
															 wi1,  NULL, 1, (NZ+2)*NX/2, FFT_PLANNING);
	if (r2c_2d == NULL) ERROR_HANDLER( ERROR_CRITICAL, "FFTW R2C_2D plan creation failed");
														   
	c2r_2d = fftw_plan_many_dft_c2r(2, n_size2D, NY / NPROC, wi1,  NULL, 1, (NZ+2)*NX/2,
														    wir1, NULL, 1, (NZ+2)*NX  , FFT_PLANNING);
	if (c2r_2d == NULL) ERROR_HANDLER( ERROR_CRITICAL, "FFTW C2R_2D plan creation failed");
	
															 
	// 1D transforms: This are actually c2c transforms, but are used for global 3D transforms.
	// We will transform forward and backward an array of logical size [NX/NPROC, NY, (NZ+2)/2] along the 2nd dimension
	// We will do NZ_COMPLEX transforms along Y. Will need a loop on NX/NPROC
	// We use &w1[NZ_COMPLEX] so that alignement check is done properly (see SIMD in fftw Documentation)
#ifdef OPENMP_SUPPORT	
	fftw_plan_with_nthreads( 1 );
#endif	
	r2c_1d = fftw_plan_many_dft(1, n_size1D, NZ_COMPLEX, &wi1[NZ_COMPLEX], NULL, NZ_COMPLEX, 1,
														 &wi1[NZ_COMPLEX], NULL, NZ_COMPLEX, 1, FFTW_FORWARD, FFT_PLANNING);
	if (r2c_1d == NULL) ERROR_HANDLER( ERROR_CRITICAL, "FFTW R2C_1D plan creation failed");
																			  
	c2r_1d = fftw_plan_many_dft(1, n_size1D, NZ_COMPLEX, &wi1[NZ_COMPLEX], NULL, NZ_COMPLEX, 1,
														 &wi1[NZ_COMPLEX], NULL, NZ_COMPLEX, 1, FFTW_BACKWARD, FFT_PLANNING);
	if (c2r_1d == NULL) ERROR_HANDLER( ERROR_CRITICAL, "FFTW C2R_1D plan creation failed");

	// Initialize threads
	pthread_attr_init(&thread_attr);
	pthread_attr_setdetachstate(&thread_attr, PTHREAD_CREATE_JOINABLE);
	
	// init transpose routines
	
	init_transpose();
	// Let's see which method is faster (with our without threads)
	
#ifndef FORCE_NO_THREADS
	use_thread = 0;
	// INIT everything, to be fair...
	gfft3_r2c_t(wir1, wir2, wir3);
	gfft3_c2r_t(wi1, wi2, wi3);
	
	MPI_Barrier(MPI_COMM_WORLD);
	timeinit1=get_walltime();
	gfft3_r2c_t(wir1, wir2, wir3);
	gfft3_c2r_t(wi1, wi2, wi3);
	MPI_Barrier(MPI_COMM_WORLD);
	timeinit1=get_walltime()-timeinit1;
	
	use_thread = 1;
	gfft3_r2c_t(wir1, wir2, wir3);
	gfft3_c2r_t(wi1, wi2, wi3);
	
	MPI_Barrier(MPI_COMM_WORLD);
	timeinit2=get_walltime();
	gfft3_r2c_t(wir1, wir2, wir3);
	gfft3_c2r_t(wi1, wi2, wi3);
	MPI_Barrier(MPI_COMM_WORLD);
	timeinit2=get_walltime()-timeinit2;
	
	
	if (timeinit2>timeinit1) {
		MPI_Printf("Gfft will NOT use threads as it's %f times slower.\n", timeinit2/timeinit1);
		use_thread = 0;
	}
	else {
		MPI_Printf("Gfft will use threads as it's %f times faster.\n", timeinit1/timeinit2);
		use_thread = 1;
	}
#else
	MPI_Printf("Gfft forced NOT to use threads\n");
	use_thread = 0;
#endif
	
	fftw_free(wi1);
	fftw_free(wi2);
	fftw_free(wi3);
	
	time_1D=0.0;
	time_2D=0.0;
	time_transpose=0.0;
	// That's in tranpose.c
	MPI_Communication_time=0.0;

	return;
}

/**********************************************************************
***********************************************************************
*********     N O    M P I    R O U T I N E S        ******************
***********************************************************************
***********************************************************************/
// These routines are essentially wrappers for fftw3 routines

#else

// Here, we assume we don't have MPI

fftw_plan	r2cfft, c2rfft;

PRECISION complex *wi1;
PRECISION *wir1;

void gfft_r2c_t(PRECISION *wrin) {
	PRECISION complex *win = (PRECISION complex *) wrin;
	fftw_execute_dft_r2c(r2cfft, wrin, win);
	return;
}

void gfft_c2r_t(PRECISION complex *win){
	PRECISION *wrin = (PRECISION *) win;
	fftw_execute_dft_c2r(c2rfft, win, wrin);
	return;
}

void gfft_r2c(PRECISION *wrin) {
	PRECISION complex *win = (PRECISION complex *) wrin;
	fftw_execute_dft_r2c(r2cfft, wrin, win);
	return;
}

void gfft_c2r(PRECISION complex *win){
	PRECISION *wrin = (PRECISION *) win;
	fftw_execute_dft_c2r(c2rfft, win, wrin);
	return;
}


void gfft3_r2c_t(PRECISION *wrin1, PRECISION *wrin2, PRECISION *wrin3){
	gfft_r2c_t(wrin1);
	gfft_r2c_t(wrin2);
	gfft_r2c_t(wrin3);
	return;
}
	
void gfft3_c2r_t(PRECISION complex *win1, PRECISION complex *win2, PRECISION complex *win3){
	gfft_c2r_t(win1);
	gfft_c2r_t(win2);
	gfft_c2r_t(win3);
	return;
}

void init_gfft() {
	
	wi1 = (PRECISION complex *) fftw_malloc( sizeof(PRECISION complex) * NTOTAL_COMPLEX);
	if (wi1 == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for wi1 allocation");
	
	wir1 = (PRECISION *) wi1;
	
#ifdef OPENMP_SUPPORT
	fftw_plan_with_nthreads( NTHREADS );
#endif
	
	r2cfft = fftw_plan_dft_r2c_3d( NX, NY, NZ, wr1, w1,  FFT_PLANNING);
	if (r2cfft == NULL) ERROR_HANDLER( ERROR_CRITICAL, "FFTW R2C plan creation failed");

	c2rfft = fftw_plan_dft_c2r_3d( NX, NY, NZ, w1,  wr1, FFT_PLANNING);
	if (c2rfft == NULL) ERROR_HANDLER( ERROR_CRITICAL, "FFTW C2R plan creation failed");

	
	fftw_free(wi1);
	
	
	return;
}

#endif

														   
	