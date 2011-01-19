#include <stdlib.h>
#include <string.h>

#include "../common.h"
#include "../gfft.h"
#include "../shear.h"
#include "../debug.h"
#include "../particles.h"

/***********************************************************/
/** 
	find the maximum of a real array of size 2*NTOTAL_COMPLEX
	@param wri array in which we want to know the maximum
*/
/***********************************************************/

double find_max(double *wri) {
	double q0;
	int i,j,k,idx;
	q0=wri[0];
	
	for(i=0 ; i < NX/NPROC ; i++) {
		for(j=0 ; j < NY ; j++) {
			for(k=0 ; k < NZ ; k++) {
#ifdef WITH_2D
				idx = j + i * (NY+2);
#else
				idx = k + j * (NZ + 2) + i * NY * (NZ + 2);
#endif
				if(q0<wri[idx]) q0=wri[idx];
			}
		}
	}
	
	return(q0);
}

/***********************************************************/
/** 
	find the minimum of a real array of size 2*NTOTAL_COMPLEX
	@param wri array in which we want to know the minimum
*/
/***********************************************************/

double find_min(double *wri) {
	double q0;
	int i,j,k,idx;
	q0=wri[0];
	
	for(i=0 ; i < NX/NPROC ; i++) {
		for(j=0 ; j < NY ; j++) {
			for(k=0 ; k < NZ ; k++) {
#ifdef WITH_2D
				idx = j + i * (NY+2);
#else
				idx = k + j * (NZ + 2) + i * NY * (NZ + 2);
#endif
				if(q0>wri[idx]) q0=wri[idx];
			}
		}
	}
	
	return(q0);
}

/***********************************************************/
/** 
	compute the correlation between 2 fields
*/
/***********************************************************/
double compute_2correlation(double *wri1, double *wri2) {
	double q0;
	int i,j,k;
	q0=0.0;
	
	for(i = 0 ; i < 2 * NTOTAL_COMPLEX ; i ++) {
		q0 += wri1[i] * wri2[i] / ((double) NTOTAL);
	}
	
	return(q0);
}

/***********************************************************/
/** 
	compute the correlation between 3 fields
*/
/***********************************************************/
double compute_3correlation(double *wri1, double *wri2, double *wri3) {
	double q0;
	int i,j,k;
	q0=0.0;
	
	for(i = 0 ; i < 2 * NTOTAL_COMPLEX ; i ++) {
		q0 += wri1[i] * wri2[i] * wri3[i] / ((double) NTOTAL);
	}
	
	return(q0);
}

/***********************************************************/
/** 
	Write statistical quantities using text format in the file
	timevar. 
	@param fldi Field structure from which the statistical quantities are derived.
	@param t Current time of the simulation
*/
/***********************************************************/

void output_timevar(const struct Field fldi,
					const double t) {
					
	FILE *ht;
	double output_var;
	static int warning_flag = 0;

	int i,j;
	
	DEBUG_START_FUNC;
	
	// Open the timevar file
	if(rank==0) {
		ht=fopen("timevar","a");
	}
	
// Transform a few basic variables.
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i] = fldi.vx[i];
		w2[i] = fldi.vy[i];
		w3[i] = fldi.vz[i];
#ifdef BOUSSINESQ
		w4[i] = fldi.th[i];
#endif
#ifdef COMPRESSIBLE
		w4[i] = fldi.d[i];
#endif
#ifdef MHD
		w5[i] = fldi.bx[i];
		w6[i] = fldi.by[i];
		w7[i] = fldi.bz[i];
#endif
	}

	gfft_c2r(w1);
	gfft_c2r(w2);
	gfft_c2r(w3);
#ifdef BOUSSINESQ	// Boussinesq and compressibility are in principle mutually exclusive.
	gfft_c2r(w4);
#endif
#ifdef COMPRESSIBLE
	gfft_c2r(w4);
	check_positivity(wr4);
#endif
#ifdef MHD
	gfft_c2r(w5);
	gfft_c2r(w6);
	gfft_c2r(w7);
#endif

	
	for( i = 0 ; i < 2*NTOTAL_COMPLEX ; i++) {
#ifndef COMPRESSIBLE
		wr1[i] = wr1[i] / ((double) NTOTAL );
		wr2[i] = wr2[i] / ((double) NTOTAL );
		wr3[i] = wr3[i] / ((double) NTOTAL );
#ifdef BOUSSINESQ
		wr4[i] = wr4[i] / ((double) NTOTAL );
#endif
#else
		wr1[i] = wr1[i] / wr4[i];		// Compute the real velocity field
		wr2[i] = wr2[i] / wr4[i];
		wr3[i] = wr3[i] / wr4[i];
		wr4[i] = wr4[i] / ((double) NTOTAL);
#endif
		
#ifdef MHD
		wr5[i] = wr5[i] / ((double) NTOTAL );
		wr6[i] = wr6[i] / ((double) NTOTAL );
		wr7[i] = wr7[i] / ((double) NTOTAL );
		
#endif
	}
	
	// loop on all the requested variables
	for( i = 0 ; i < param.timevar_vars.length ; i++ ) {
	
		if(!strcmp(param.timevar_vars.name[i],"t")) {
			// current time
			output_var=t;
		}
		else if(!strcmp(param.timevar_vars.name[i],"ev")) {
			// kinetic energy
#ifndef COMPRESSIBLE
			output_var = energy(fldi.vx) + energy(fldi.vy)+energy(fldi.vz);
#else
			output_var = 0.5 * ( compute_3correlation(wr1, wr1, wr4) +
								 compute_3correlation(wr2, wr2, wr4) +
								 compute_3correlation(wr3, wr3, wr4) );
#endif

			reduce(&output_var,1);
		}
		else if(!strcmp(param.timevar_vars.name[i],"vxmax")) {
			// maximum of vx component
			output_var=find_max(wr1);
			reduce(&output_var,2);
		}
		else if(!strcmp(param.timevar_vars.name[i],"vymax")) {
			// maximum of vy component
			output_var=find_max(wr2);
			reduce(&output_var,2);
		}
		else if(!strcmp(param.timevar_vars.name[i],"vzmax")) {
			// maximum of vz component
			output_var=find_max(wr3);
			reduce(&output_var,2);
		}
		else if(!strcmp(param.timevar_vars.name[i],"vxmin")) {
			// minimum of vx component
			output_var=find_min(wr1);
			reduce(&output_var,3);
		}
		else if(!strcmp(param.timevar_vars.name[i],"vymin")) {
			// minimum of vy component
			output_var=find_min(wr2);
			reduce(&output_var,3);
		}
		else if(!strcmp(param.timevar_vars.name[i],"vzmin")) {
			// minimum of vz component
			output_var=find_min(wr3);
			reduce(&output_var,3);
		}
		else if(!strcmp(param.timevar_vars.name[i],"vxvy")) {
			// incompressible reynolds stress
			output_var=compute_2correlation(wr1,wr2);
			reduce(&output_var,1);
		}
		else if(!strcmp(param.timevar_vars.name[i],"hv")) {
			// kinetic helicity
			// Compute vector potential
			for( j = 0 ; j < NTOTAL_COMPLEX ; j++) {
				w10[j] = I * ik2t[j] * (ky[j] * fldi.vz[j] - kz[j] * fldi.vy[j] );
				w11[j] = I * ik2t[j] * (kz[j] * fldi.vx[j] - kxt[j]* fldi.vz[j] );
				w12[j] = I * ik2t[j] * (kxt[j]* fldi.vy[j] - ky[j] * fldi.vx[j] );
			}
	
			gfft_c2r(w10);
			gfft_c2r(w11);
			gfft_c2r(w12);
	
			for( j = 0 ; j < 2*NTOTAL_COMPLEX ; j++) {
				wr10[j] = wr10[j] / ((double) NTOTAL );
				wr11[j] = wr11[j] / ((double) NTOTAL );
				wr12[j] = wr12[j] / ((double) NTOTAL );
			}

			output_var=compute_2correlation(wr10,wr1)+compute_2correlation(wr11,wr2)+compute_2correlation(wr12,wr3);
			reduce(&output_var,1);
		}
		
		else if(!strcmp(param.timevar_vars.name[i],"w2")) {
			// total enstrophy
			// Compute vorticity
			for( j = 0 ; j < NTOTAL_COMPLEX ; j++) {
				w10[j] = I * (ky[j]  * fldi.vz[j] - kz[j]  * fldi.vy[j]);
				w11[j] = I * (kz[j]  * fldi.vx[j] - kxt[j] * fldi.vz[j]);
				w12[j] = I * (kxt[j] * fldi.vy[j] - ky[j]  * fldi.vx[j]);
			}
			output_var=energy(w10)+energy(w11)+energy(w12);
			reduce(&output_var,1);
		}

		

#ifdef MHD
		else if(!strcmp(param.timevar_vars.name[i],"em")) {
			// magnetic energy
			output_var = energy(fldi.bx) + energy(fldi.by)+energy(fldi.bz);
			reduce(&output_var,1);
		}
		else if(!strcmp(param.timevar_vars.name[i],"bxmax")) {
			// maximum of bx component
			output_var=find_max(wr5);
			reduce(&output_var,2);
		}
		else if(!strcmp(param.timevar_vars.name[i],"bymax")) {
			// maximum of by component
			output_var=find_max(wr6);
			reduce(&output_var,2);
		}
		else if(!strcmp(param.timevar_vars.name[i],"bzmax")) {
			// maximum of bz component
			output_var=find_max(wr7);
			reduce(&output_var,2);
		}
		else if(!strcmp(param.timevar_vars.name[i],"bxmin")) {
			// minimum of bx component
			output_var=find_min(wr5);
			reduce(&output_var,3);
		}
		else if(!strcmp(param.timevar_vars.name[i],"bymin")) {
			// minimum of by component
			output_var=find_min(wr6);
			reduce(&output_var,3);
		}
		else if(!strcmp(param.timevar_vars.name[i],"bzmin")) {
			// minimum of bz component
			output_var=find_min(wr7);
			reduce(&output_var,3);
		}
		else if(!strcmp(param.timevar_vars.name[i],"bxby")) {
			// incompressible maxwell stress
			output_var=compute_2correlation(wr5,wr6);
			reduce(&output_var,1);
		}
		else if(!strcmp(param.timevar_vars.name[i],"hc")) {
			// cross helicity (u.b)
			output_var=compute_2correlation(wr1,wr5)+compute_2correlation(wr2,wr6)+compute_2correlation(wr3,wr7);
			reduce(&output_var,1);
		}
		else if(!strcmp(param.timevar_vars.name[i],"hm")) {
			// magnetic helicity
			// Compute vector potential
			for( j = 0 ; j < NTOTAL_COMPLEX ; j++) {
				w10[j] = I * ik2t[j] * (ky[j] * fldi.bz[j] - kz[j] * fldi.by[j] );
				w11[j] = I * ik2t[j] * (kz[j] * fldi.bx[j] - kxt[j]* fldi.bz[j] );
				w12[j] = I * ik2t[j] * (kxt[j]* fldi.by[j] - ky[j] * fldi.bx[j] );
			}
	
			gfft_c2r(w10);
			gfft_c2r(w11);
			gfft_c2r(w12);
	
			for( j = 0 ; j < 2*NTOTAL_COMPLEX ; j++) {
				wr10[j] = wr10[j] / ((double) NTOTAL );
				wr11[j] = wr11[j] / ((double) NTOTAL );
				wr12[j] = wr12[j] / ((double) NTOTAL );
			}

			output_var=compute_2correlation(wr10,wr5)+compute_2correlation(wr11,wr6)+compute_2correlation(wr12,wr7);
			reduce(&output_var,1);
		}
		else if(!strcmp(param.timevar_vars.name[i],"j2")) {
			// total current rms
			// Compute current
			for( j = 0 ; j < NTOTAL_COMPLEX ; j++) {
				w10[j] = I * (ky[j]  * fldi.bz[j] - kz[j]  * fldi.by[j]);
				w11[j] = I * (kz[j]  * fldi.bx[j] - kxt[j] * fldi.bz[j]);
				w12[j] = I * (kxt[j] * fldi.by[j] - ky[j]  * fldi.bx[j]);
			}
			output_var=energy(w10)+energy(w11)+energy(w12);
			reduce(&output_var,1);
		}
		
		else if(!strcmp(param.timevar_vars.name[i],"az2")) {
			// Square of the vertical component of the vector potential
			for( j = 0 ; j < NTOTAL_COMPLEX ; j++) {
				w12[j] = I * ik2t[j] * (kxt[j]* fldi.by[j] - ky[j] * fldi.bx[j] );
			}
			
			output_var=energy(w12);
			
			reduce(&output_var,1);
		}
		
		else if(!strcmp(param.timevar_vars.name[i],"vxaz")) {
			// Source term in the streamfunction equation
			for( j = 0 ; j < NTOTAL_COMPLEX ; j++) {
				w12[j] = I * ik2t[j] * (kxt[j]* fldi.by[j] - ky[j] * fldi.bx[j] );
			}
			
			gfft_c2r(w12);
			
			for( j = 0 ; j < 2*NTOTAL_COMPLEX ; j++) {
				wr12[j] = wr12[j] / ((double) NTOTAL );
			}
			
			output_var=compute_2correlation(wr12,wr1);
			reduce(&output_var,1);
			
		}
#endif
#ifdef BOUSSINESQ
		else if(!strcmp(param.timevar_vars.name[i],"et")) {
			// thermal energy
			output_var = energy(fldi.th);
			reduce(&output_var,1);
		}
		else if(!strcmp(param.timevar_vars.name[i],"thmax")) {
			// maximum of th
			output_var=find_max(wr4);
			reduce(&output_var,2);
		}
		else if(!strcmp(param.timevar_vars.name[i],"thmin")) {
			// minimum of th
			output_var=find_min(wr4);
			reduce(&output_var,3);
		}
		else if(!strcmp(param.timevar_vars.name[i],"thvx")) {
			// turbulent heat flux in x
			output_var=compute_2correlation(wr4,wr1);
			reduce(&output_var,1);
		}
		else if(!strcmp(param.timevar_vars.name[i],"thvz")) {
			// turbulent heat flux in z
			output_var=compute_2correlation(wr4,wr3);
			reduce(&output_var,1);
		}
#endif
#ifdef COMPRESSIBLE
		else if(!strcmp(param.timevar_vars.name[i],"dmax")) {
			// maximum of density
			output_var=find_max(wr4);
			reduce(&output_var,2);
		}
		else if(!strcmp(param.timevar_vars.name[i],"dmin")) {
			// minimum of density
			output_var=find_min(wr4);
			reduce(&output_var,3);
		}
		else if(!strcmp(param.timevar_vars.name[i],"dvxvy")) {
			// compressible reynolds stress
			output_var=compute_3correlation(wr1,wr2,wr4);
			reduce(&output_var,1);
		}
#endif

		else {
			if(!warning_flag) {
				ERROR_HANDLER(ERROR_WARNING,"Unable to produce the requested data\n");
				MPI_Printf("Timevar output string ''%s'' is unknown\n",param.timevar_vars.name[i]);
				warning_flag=1;
			}
			output_var=0.0;
		}
		
		// Write the variable in the timevar file
		if(rank==0) {
			fprintf(ht,"%08e\t",output_var);
			if(ferror(ht)) ERROR_HANDLER( ERROR_CRITICAL, "Error writing timevar file");
		}
		
#ifdef MPI_SUPPORT
		MPI_Barrier(MPI_COMM_WORLD);
#endif
	}
		
	if(rank==0) {
		fprintf(ht,"\n");
		fclose(ht);
	}

	DEBUG_END_FUNC;
	
	return;
}

/**************************************************************************************/
/** 
	Remove the timevar file (if exists) to start from a fresh one.
*/
/**************************************************************************************/
void init_timevar() {
	FILE* ht;
	int i;
	DEBUG_START_FUNC;
	
	if(rank==0) {
		ht=fopen("timevar","w");
		// print a line with the fields we are going to write
		
		for( i = 0 ; i < param.timevar_vars.length ; i++) {
			fprintf(ht,"%s\t\t",param.timevar_vars.name[i]);
		}
		fprintf(ht,"\n");
		
		fclose(ht);
#ifdef WITH_PARTICLES
		ht=fopen("partvar","w");
		fclose(ht);
#endif
	}

	DEBUG_END_FUNC;
	return;
}


