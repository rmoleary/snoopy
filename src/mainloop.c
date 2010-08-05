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
#include "symmetries.h"
#ifdef BOUNDARY_C
#include "boundary.h"
#endif
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

#ifdef WITH_PARTICLES
	gamma_v += 1.0 / (fabs(param.particles_stime) * param.safety_source );
	
	maxfx=0.0;
	maxfy=0.0;
	maxfz=0.0;
	
	for(i = 0 ; i < param.particles_n/NPROC ; i++) {
		if( fabs( fld.part[i].vx ) > maxfx ) maxfx = fabs( fld.part[i].vx );
		if( fabs( fld.part[i].vy ) > maxfy ) maxfy = fabs( fld.part[i].vy );
		if( fabs( fld.part[i].vz ) > maxfz ) maxfz = fabs( fld.part[i].vz );
	}
#ifdef MPI_SUPPORT
	reduce(&maxfx,2);
	reduce(&maxfy,2);
	reduce(&maxfz,2);
#endif
	
	gamma_v += param.lx/(NX)*maxfx+param.ly/(NY)*maxfy+param.lz/(NZ)*maxfz;
	
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
	
	allocate_field(&dfld);
	allocate_field(&fld1);
	
// Init the Runge-Kutta timestepping
// Values coming from Brandenburg (2001) page 8

	gammaRK[0] = 8.0 / 15.0;
	gammaRK[1] = 5.0 / 12.0;
	gammaRK[2] = 3.0 / 4.0;
	
	xiRK[0] = -17.0 / 60.0;
	xiRK[1] = -5.0 / 12.0;

	DEBUG_END_FUNC;
	
	return;
}

void finish_mainloop() {
	deallocate_field(&fld1);
	deallocate_field(&dfld);
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
	int i,n,nloop;
	
	DEBUG_START_FUNC;
	
	init_mainloop();
	nloop=0;
	
	if(param.restart) {
#ifdef DEBUG
		MPI_Printf("Reading dump file\n");
#endif
		read_dump(fld,&t,OUTPUT_DUMP);
	}
	else {
		t = t_start;
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
		
		timestep(dfld, fld, pressure, t, tremap, dt );
		
#ifdef _OPENMP
		#pragma omp parallel private(i,n) 
		{
#endif
		for( n = 0 ; n < fld.nfield ; n++) {
#ifdef _OPENMP
		#pragma omp for schedule(static)	
#endif
			for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
				fld.farray[n][i] = fld.farray[n][i] + gammaRK[0] * dfld.farray[n][i] * dt;
				fld1.farray[n][i] = fld.farray[n][i] + xiRK[0] * dfld.farray[n][i] * dt;
			}
		}
#ifdef WITH_PARTICLES
#ifdef _OPENMP
		#pragma omp for schedule(static)	
#endif
		for( i = 0 ; i < param.particles_n ; i++) {
			fld.part[i].x = fld.part[i].x + gammaRK[0] * dfld.part[i].x * dt;
			fld.part[i].y = fld.part[i].y + gammaRK[0] * dfld.part[i].y * dt;
			fld.part[i].z = fld.part[i].z + gammaRK[0] * dfld.part[i].z * dt;
			
			fld.part[i].vx = fld.part[i].vx + gammaRK[0] * dfld.part[i].vx * dt;
			fld.part[i].vy = fld.part[i].vy + gammaRK[0] * dfld.part[i].vy * dt;
			fld.part[i].vz = fld.part[i].vz + gammaRK[0] * dfld.part[i].vz * dt;
			
			fld1.part[i].x = fld.part[i].x + xiRK[0] * dfld.part[i].x * dt;
			fld1.part[i].y = fld.part[i].y + xiRK[0] * dfld.part[i].y * dt;
			fld1.part[i].z = fld.part[i].z + xiRK[0] * dfld.part[i].z * dt;
			
			fld1.part[i].vx = fld.part[i].vx + xiRK[0] * dfld.part[i].vx * dt;
			fld1.part[i].vy = fld.part[i].vy + xiRK[0] * dfld.part[i].vy * dt;
			fld1.part[i].vz = fld.part[i].vz + xiRK[0] * dfld.part[i].vz * dt;

		}
#endif
#ifdef _OPENMP
		}
#endif
		
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
		timestep(dfld, fld, NULL, t+gammaRK[0]*dt, tremap+gammaRK[0]*dt, dt);

#ifdef _OPENMP
		#pragma omp parallel private(i,n) 
		{
#endif
		for( n = 0 ; n < fld.nfield ; n++) {
#ifdef _OPENMP
		#pragma omp for schedule(static)	
#endif
			for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
				fld.farray[n][i] = fld1.farray[n][i] + gammaRK[1] * dfld.farray[n][i] * dt;
				fld1.farray[n][i] = fld.farray[n][i] + xiRK[1] * dfld.farray[n][i] * dt;
			}
		}
#ifdef WITH_PARTICLES
#ifdef _OPENMP
		#pragma omp for schedule(static)	
#endif
		for( i = 0 ; i < param.particles_n ; i++) {
			fld.part[i].x = fld1.part[i].x + gammaRK[1] * dfld.part[i].x * dt;
			fld.part[i].y = fld1.part[i].y + gammaRK[1] * dfld.part[i].y * dt;
			fld.part[i].z = fld1.part[i].z + gammaRK[1] * dfld.part[i].z * dt;
			
			fld.part[i].vx = fld1.part[i].vx + gammaRK[1] * dfld.part[i].vx * dt;
			fld.part[i].vy = fld1.part[i].vy + gammaRK[1] * dfld.part[i].vy * dt;
			fld.part[i].vz = fld1.part[i].vz + gammaRK[1] * dfld.part[i].vz * dt;
			
			fld1.part[i].x = fld.part[i].x + xiRK[1] * dfld.part[i].x * dt;
			fld1.part[i].y = fld.part[i].y + xiRK[1] * dfld.part[i].y * dt;
			fld1.part[i].z = fld.part[i].z + xiRK[1] * dfld.part[i].z * dt;
			
			fld1.part[i].vx = fld.part[i].vx + xiRK[1] * dfld.part[i].vx * dt;
			fld1.part[i].vy = fld.part[i].vy + xiRK[1] * dfld.part[i].vy * dt;
			fld1.part[i].vz = fld.part[i].vz + xiRK[1] * dfld.part[i].vz * dt;
		}
#endif
#ifdef _OPENMP
		}
#endif

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
		timestep(dfld, fld, NULL, t + (gammaRK[0] + xiRK[0] + gammaRK[1]) * dt, tremap + (gammaRK[0] + xiRK[0] + gammaRK[1]) * dt, dt);

#ifdef _OPENMP
		#pragma omp parallel private(i,n) 
		{
#endif
		for( n = 0 ; n < fld.nfield ; n++) {
#ifdef _OPENMP
		#pragma omp for schedule(static)	
#endif
			for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
				fld.farray[n][i] = fld1.farray[n][i] + gammaRK[2] * dfld.farray[n][i] * dt;
			}
		}
#ifdef WITH_PARTICLES
#ifdef _OPENMP
		#pragma omp for schedule(static)	
#endif
		for( i = 0 ; i < param.particles_n ; i++) {
			fld.part[i].x = fld1.part[i].x + gammaRK[2] * dfld.part[i].x * dt;
			fld.part[i].y = fld1.part[i].y + gammaRK[2] * dfld.part[i].y * dt;
			fld.part[i].z = fld1.part[i].z + gammaRK[2] * dfld.part[i].z * dt;
			
			fld.part[i].vx = fld1.part[i].vx + gammaRK[2] * dfld.part[i].vx * dt;
			fld.part[i].vy = fld1.part[i].vy + gammaRK[2] * dfld.part[i].vy * dt;
			fld.part[i].vz = fld1.part[i].vz + gammaRK[2] * dfld.part[i].vz * dt;
		}
#endif

#ifdef _OPENMP
		}
#endif

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
			for( n = 0 ; n < fld.nfield ; n++) {
				remap(fld.farray[n]);
			}
			if(param.output_pressure)
				remap(pressure);
		}
#endif
		kvolve(tremap);
#endif
		// Symmetries cleaning
		if(param.force_symmetries) {
			if(!(nloop % param.symmetries_step)) enforce_complex_symm(fld);
		}
		
		// Divergence cleaning
#ifdef ANELASTIC
		projector_anelastic(fld.vx,fld.vy,fld.vz);
#else
		projector(fld.vx,fld.vy,fld.vz);
#endif

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
