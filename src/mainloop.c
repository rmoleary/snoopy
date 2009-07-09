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



#include <stdlib.h>

#include "common.h"
#include "timestep.h"
#include "output.h"
#include "interface.h"
#include "gfft.h"
#include "shear.h"
#include "transpose.h"
#include "debug.h"

struct Field			dfld, fld1;

double complex		gammaRK[3];
double complex		xiRK[2];

double forcing_last_time;

double newdt(double tremap) {

	int i;
	double gamma_v;
	double maxfx   , maxfy, maxfz;
#ifdef MHD
	double gamma_b;
	double maxbx   , maxby, maxbz;
#endif
	double dt;
	
	DEBUG_START_FUNC;
	
#ifdef _OPENMP
	#pragma omp parallel for private(i) schedule(static)
#endif
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i] =  fld.vx[i];
		w2[i] =  fld.vy[i];
		w3[i] =  fld.vz[i];
	}

	gfft_c2r_t(w1);
	gfft_c2r_t(w2);
	gfft_c2r_t(w3);
	
	maxfx=0.0;
	maxfy=0.0;
	maxfz=0.0;

	for( i = 0 ; i < NTOTAL_COMPLEX * 2 ; i++) {
		if( fabs( wr1[i] ) > maxfx ) maxfx = fabs( wr1[i] );
		if( fabs( wr2[i] ) > maxfy ) maxfy = fabs( wr2[i] );
		if( fabs( wr3[i] ) > maxfz ) maxfz = fabs( wr3[i] );
		
	}

	maxfx = maxfx / ((double) NTOTAL);
	maxfy = maxfy / ((double) NTOTAL);
	maxfz = maxfz / ((double) NTOTAL);
	
#ifdef MPI_SUPPORT
	reduce(&maxfx,2);
	reduce(&maxfy,2);
	reduce(&maxfz,2);
#endif
	
	gamma_v = (kxmax + fabs(tremap)*kymax) * maxfx + kymax * maxfy + kzmax * maxfz;
#ifdef WITH_ROTATION
	gamma_v += fabs(param.omega) / param.safety_source;
#endif

#ifdef WITH_SHEAR
	gamma_v += fabs(param.shear) / param.safety_source;
#endif
#ifdef BOUSSINESQ
	gamma_v += pow(fabs(param.N2), 0.5) / param.safety_source;
#endif
#ifdef TIME_DEPENDANT_SHEAR
	gamma_v += fabs(param.omega_shear) / param.safety_source;
#endif
#ifdef MHD

	/* Compute the magnetic CFL condition */
#ifdef _OPENMP
	#pragma omp parallel for private(i) schedule(static)	
#endif
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i] =  fld.bx[i];
		w2[i] =  fld.by[i];
		w3[i] =  fld.bz[i];
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
	
	dt = param.cfl / (gamma_v + gamma_b);
#else
	dt = param.cfl / gamma_v;
#endif

#ifdef DEBUG
#ifdef MHD
	MPI_Printf("newdt: maxbx=%e, maxby=%e, maxbz=%e\n",maxbx,maxby, maxbz);
#endif
	MPI_Printf("newdt: maxfx=%e, maxfy=%e, maxfz=%e, dt=%e\n",maxfx,maxfy, maxfz, dt);
#endif
	DEBUG_END_FUNC;
	return(dt);
}			   			   
		
void init_mainloop() {

	DEBUG_START_FUNC;
	
	dfld.vx = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (dfld.vx == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for dfld.vx allocation");
	
	dfld.vy = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (dfld.vy == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for dfld.vy allocation");
	
	dfld.vz = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (dfld.vz == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for dfld.vz allocation");
	
	fld1.vx = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (fld1.vx == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for fld1.vx allocation");
	
	fld1.vy = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (fld1.vy == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for fld1.vy allocation");
	
	fld1.vz = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (fld1.vz == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for fld1.vz allocation");

#ifdef BOUSSINESQ
	dfld.th = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (dfld.th == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for dfld.th allocation");
	
	fld1.th = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (fld1.th == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for fld1.th allocation");
#endif
#ifdef MHD
	dfld.bx = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (dfld.bx == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for dfld.bx allocation");
	
	dfld.by = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (dfld.by == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for dfld.by allocation");
	
	dfld.bz = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (dfld.bz == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for dfld.bz allocation");
	
	fld1.bx = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (fld1.bx == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for fld1.bx allocation");
	
	fld1.by = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (fld1.by == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for fld1.by allocation");
	
	fld1.bz = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (fld1.bz == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for fld1.bz allocation");
#endif
	
// Init the Runge-Kutta timestepping

	gammaRK[0] = 8.0 / 15.0;
	gammaRK[1] = 5.0 / 12.0;
	gammaRK[2] = 3.0 / 4.0;
	
	xiRK[0] = -17.0 / 60.0;
	xiRK[1] = -5.0 / 12.0;

	DEBUG_END_FUNC;
	
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
#endif
	return;
}

/***************************************************************/
/**
	Integrate in time the field stored in fld (found and initialized in common.c/common.h)
	from t_start to t_end. Outputs are done according to gvars.h
	
	@param t_start: initial time of the simulation (usually 0...)
	@param t_end: final time of the simulation (will stop precisely at that time).
*/
/***************************************************************/
void mainloop(double t_start, double t_end) {
	double		dt = 0.0;
	double	    t = 0.0;
	double		tremap = 0.0;
	
	double timer_end, timer_start;
	int i,nloop;
	
	DEBUG_START_FUNC;
	
	init_mainloop();
	nloop=0;
	
	if(param.restart) {
#ifdef DEBUG
		MPI_Printf("Reading dump file\n");
#endif
		if(read_dump(fld,&t)) {
			MPI_Printf("Mainloop: No dump found, using normal initialization.\n");
			t = t_start;
			clear_timevar();
			output(t);
		}
	}
	else {
		t = t_start;
		clear_timevar();
		output(t);
	}

#ifdef WITH_SHEAR	
	tremap = time_shift(t);
	kvolve(tremap);
#else
	tremap = 0.0;
#endif
	
	timer_start = get_c_time();
	
	while (t < t_end) {
#ifdef DEBUG
		MPI_Printf("Begining of loop:\n");
		MPI_Printf("fld:\n");
		D_show_all(fld);
		MPI_Printf("**************************************************************************************\n");
#endif

		nloop++;
		if(!(nloop % param.interface_check)) check_interface(fld,t,dt,nloop,timer_start);
		
		dt = newdt(tremap);
		// Let's try to stop exactly at t_final
		if(dt > (t_end - t)) dt = t_end - t;
		
		// Stop if elpased time is larger than MAX_ELAPSED_TIME (in hours)
		if((get_c_time()-timer_start) > 3600 * param.max_t_elapsed) {
			MPI_Printf("Maximum elapsed time reached. Terminating.\n");
			dump_immediate(t);
			break;
		}
		
		// This is an order 3 runge Kutta scheme with low storage
		
		// 1st RK3 step
		
		timestep(dfld, fld, pressure, t, dt );
		
#ifdef _OPENMP
		#pragma omp parallel for private(i) schedule(static)	
#endif
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
#ifdef TIME_DEPENDANT_SHEAR
		kvolve(time_shift(t+gammaRK[0]*dt));
#else
		kvolve(tremap+gammaRK[0]*dt);
#endif
#endif
		
		timestep(dfld, fld, NULL, t+gammaRK[0]*dt, dt);

#ifdef _OPENMP
		#pragma omp parallel for private(i) schedule(static)	
#endif
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
#ifdef TIME_DEPENDANT_SHEAR
		kvolve(time_shift(t + (gammaRK[0] + xiRK[0] + gammaRK[1]) * dt));
#else
		kvolve(tremap + (gammaRK[0] + xiRK[0] + gammaRK[1]) * dt );
#endif
#endif
		
		timestep(dfld, fld, NULL, t + (gammaRK[0] + xiRK[0] + gammaRK[1]) * dt, dt);

#ifdef _OPENMP
		#pragma omp parallel for private(i) schedule(static)	
#endif			
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
#ifdef TIME_DEPENDANT_SHEAR	
		tremap = time_shift(t);
#else
		tremap = tremap + dt;
		
		// Check if a remap is needed
		if(tremap > param.ly / (2.0 * param.shear * param.lx)) {
			tremap = time_shift(t);    // Recompute tremap from current time, assuming all the remaps have been done
			remap(fld.vx);
			remap(fld.vy);
			remap(fld.vz);
			if(param.output_pressure)
				remap(pressure);
#ifdef BOUSSINESQ
			remap(fld.th);
#endif
#ifdef MHD
			remap(fld.bx);
			remap(fld.by);
			remap(fld.bz);
#endif
		}
#endif
		kvolve(tremap);
#endif
		// Symmetries cleaning
		if(param.force_symmetries) {
			if(!(nloop % param.symmetries_step)) enforce_symm(fld);
		}
		
		// Divergence cleaning
		projector(fld.vx,fld.vy,fld.vz);
#ifdef MHD
		projector(fld.bx,fld.by,fld.bz);
#endif
		// The boundary conditions arises naturally from the initial conditions (the relevant symmetries are conserved by the eq. of motion)
		// We keep this instruction here to enforce these boundary conditions at the end of each loop to remove numerical noise.
		// Nevertheless, it is not required to call it so often...
#ifdef BOUNDARY_C
		boundary_c(fld);
#endif

		output(t);
	}
	timer_end=get_c_time();
	MPI_Printf("mainloop finished in %d loops and %f seconds (%f sec/loop)\n",nloop,timer_end-timer_start,(timer_end-timer_start)/nloop);
#ifdef MPI_SUPPORT
	MPI_Printf("Time used for transpose: %f seconds, or %f pc of total computation time\n",read_transpose_timer(), read_transpose_timer()/(timer_end-timer_start)*100.0);
#endif

	finish_mainloop();
	return;

}
