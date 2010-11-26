#include <stdlib.h>

#include "../common.h"
#include "../gfft.h"
#include "../shear.h"
#include "../debug.h"
#include "../particles.h"

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
	double vxmax, vxmin, vymax, vymin, vzmax, vzmin, thmin, thmax;
	double bxmax, bxmin, bymax, bymin, bzmax, bzmin;
	double vort2, curr2;
	double energy_total;	
	double energy_mag;
	double reynolds_stress;
	double maxwell_stress;
	double helicity;
#ifdef WITH_LINEAR_TIDE
	double tide_stress;
#endif

	int i;
	
	DEBUG_START_FUNC;
	
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i] = fldi.vx[i];
		w2[i] = fldi.vy[i];
		w3[i] = fldi.vz[i];
#ifdef BOUSSINESQ
		w4[i] = fldi.th[i];
#endif
#ifdef MHD
		w5[i] = fldi.bx[i];
		w6[i] = fldi.by[i];
		w7[i] = fldi.bz[i];
#endif
#ifdef WITH_LINEAR_TIDE
		w5[i] = fldi.tvx[i];	// When WITH_LINEAR_TIDE is on, the statistics of the tide replaces the statistics of B
		w6[i] = fldi.tvy[i];
		w7[i] = fldi.tvz[i];
#endif
	}

	energy_total = energy(w1) + energy(w2) + energy(w3);
	reduce(&energy_total,1);
#if defined(MHD) || defined(WITH_LINEAR_TIDE)
	energy_mag = energy(w5) + energy(w6) + energy(w7);
	reduce(&energy_mag,1);
#else
	energy_mag=0.0;
#endif
	
	gfft_c2r(w1);
	gfft_c2r(w2);
	gfft_c2r(w3);
#ifdef BOUSSINESQ
	gfft_c2r(w4);
#endif
#if defined(MHD) || defined(WITH_LINEAR_TIDE)
	gfft_c2r(w5);
	gfft_c2r(w6);
	gfft_c2r(w7);
#endif

	for( i = 0 ; i < 2*NTOTAL_COMPLEX ; i++) {
		wr1[i] = wr1[i] / ((double) NTOTAL );
		wr2[i] = wr2[i] / ((double) NTOTAL );
		wr3[i] = wr3[i] / ((double) NTOTAL );
#ifdef BOUSSINESQ
		wr4[i] = wr4[i] / ((double) NTOTAL );
#endif
#if defined(MHD) || defined(WITH_LINEAR_TIDE)
		wr5[i] = wr5[i] / ((double) NTOTAL );
		wr6[i] = wr6[i] / ((double) NTOTAL );
		wr7[i] = wr7[i] / ((double) NTOTAL );
		
#endif
		
	}
	
	// w1, w2, w3 contains vx vy vz
	// w4 contains the entropy perturbations
		
	vxmax=0.0;
	vxmin=0.0;
	vymax=0.0;
	vymin=0.0;
	vzmax=0.0;
	vzmin=0.0;
	thmax=0.0;
	thmin=0.0;
	reynolds_stress=0.0;
	bxmax=0.0;
	bxmin=0.0;
	bymax=0.0;
	bymin=0.0;
	bzmax=0.0;
	bzmin=0.0;
	maxwell_stress=0.0;

	
		
	for(i = 0 ; i < 2 * NTOTAL_COMPLEX ; i ++) {
		if(vxmax < wr1[i]) vxmax = wr1[i];
		if(vxmin > wr1[i]) vxmin = wr1[i];
		if(vymax < wr2[i]) vymax = wr2[i];
		if(vymin > wr2[i]) vymin = wr2[i];
		if(vzmax < wr3[i]) vzmax = wr3[i];
		if(vzmin > wr3[i]) vzmin = wr3[i];
		
		reynolds_stress += wr1[i] * wr2[i] / ((double) NTOTAL);
#ifdef BOUSSINESQ
		if(thmax < wr4[i]) thmax = wr4[i];
		if(thmin > wr4[i]) thmin = wr4[i];
#endif
#if defined(MHD) || defined(WITH_LINEAR_TIDE)
		if(bxmax < wr5[i]) bxmax = wr5[i];
		if(bxmin > wr5[i]) bxmin = wr5[i];
		if(bymax < wr6[i]) bymax = wr6[i];
		if(bymin > wr6[i]) bymin = wr6[i];
		if(bzmax < wr7[i]) bzmax = wr7[i];
		if(bzmin > wr7[i]) bzmin = wr7[i];
#endif
#ifdef MHD
		maxwell_stress += wr5[i] * wr6[i] / ((double) NTOTAL);
#endif
#ifdef WITH_LINEAR_TIDE
		maxwell_stress += (wr5[i] * wr2[i] + wr6[i] * wr1[i]) / ((double) NTOTAL);
#endif

	}
	
	reduce(&vxmax,2);
	reduce(&vxmin,3);
	reduce(&vymax,2);
	reduce(&vymin,3);
	reduce(&vzmax,2);
	reduce(&vzmin,3);
	
	reduce(&reynolds_stress,1);
#ifdef BOUSSINESQ
	reduce(&thmax,2);
	reduce(&thmin,3);
#endif
#ifdef MHD
	reduce(&bxmax,2);
	reduce(&bxmin,3);
	reduce(&bymax,2);
	reduce(&bymin,3);
	reduce(&bzmax,2);
	reduce(&bzmin,3);
	
	reduce(&maxwell_stress,1);
#endif

// compute the magnetic helicity
	helicity = 0.0;
#ifdef MHD
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i] = I * ik2t[i] * (ky[i] * fldi.bz[i] - kz[i] * fldi.by[i] );
		w2[i] = I * ik2t[i] * (kz[i] * fldi.bx[i] - kxt[i]* fldi.bz[i] );
		w3[i] = I * ik2t[i] * (kxt[i]* fldi.by[i] - ky[i] * fldi.bx[i] );
	}
	
	gfft_c2r(w1);
	gfft_c2r(w2);
	gfft_c2r(w3);
	
	for( i = 0 ; i < 2*NTOTAL_COMPLEX ; i++) {
		wr1[i] = wr1[i] / ((double) NTOTAL );
		wr2[i] = wr2[i] / ((double) NTOTAL );
		wr3[i] = wr3[i] / ((double) NTOTAL );
	}
	
	for(i = 0 ; i < 2 * NTOTAL_COMPLEX ; i ++) {
		helicity += (wr1[i] * wr5[i] + wr2[i] * wr6[i] + wr3[i] * wr7[i]) / ((double) NTOTAL);
	}
	
	reduce(&helicity,1);
#endif
	vort2=0.0;
	curr2=0.0;
	
// Compute vorticity and currents
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i] = I * (ky[i] * fldi.vz[i] - kz[i] * fldi.vy[i]);
		w2[i] = I * (kz[i] * fldi.vx[i] - kxt[i] * fldi.vz[i]);
		w3[i] = I * (kxt[i] * fldi.vy[i] - ky[i] * fldi.vx[i]);
#ifdef MHD
		w5[i] = I * (ky[i] * fldi.bz[i] - kz[i] * fldi.by[i]);
		w6[i] = I * (kz[i] * fldi.bx[i] - kxt[i] * fldi.bz[i]);
		w7[i] = I * (kxt[i] * fldi.by[i] - ky[i] * fldi.bx[i]);
#endif
	}
	
	vort2 = energy(w1) + energy(w2) + energy(w3);
	reduce(&vort2, 1);
	
#ifdef MHD
	curr2 = energy(w5) + energy(w6) + energy(w7);
	reduce(&curr2, 1);
#endif
	
	if(rank==0) {
		ht=fopen("timevar","a");
		fprintf(ht,"%08e\t",t);
		fprintf(ht,"%08e\t%08e\t",energy_total,energy_mag);
		fprintf(ht,"%08e\t%08e\t%08e\t%08e\t%08e\t%08e\t",vxmax,vxmin,vymax,vymin,vzmax,vzmin);
		fprintf(ht,"%08e\t",reynolds_stress);
		fprintf(ht,"%08e\t%08e\t%08e\t%08e\t%08e\t%08e\t",bxmax,bxmin,bymax,bymin,bzmax,bzmin);
		fprintf(ht,"%08e\t",maxwell_stress);
		fprintf(ht,"%08e\t%08e\t",thmax,thmin);
		fprintf(ht,"%08e\t%08e\t%08e",vort2,curr2,helicity);
#ifdef TIME_DEPENDANT_SHEAR
		fprintf(ht,"\t%08e", cos(param.omega_shear * t) / param.shear );		// This parameter times the Reynolds stress leads to the instantenous "turbulent" viscosity
#endif
		fprintf(ht,"\n");
		
		if(ferror(ht)) ERROR_HANDLER( ERROR_CRITICAL, "Error writing timevar file");
		
		fclose(ht);
#ifdef MPI_SUPPORT
		MPI_Barrier(MPI_COMM_WORLD);
	}
	else	MPI_Barrier(MPI_COMM_WORLD);
#else
	}
#endif

	DEBUG_END_FUNC;
	
	return;
}

/**************************************************************************************/
/** 
	Remove the timevar file (if exists) to start from a fresh one.
*/
/**************************************************************************************/
void clear_timevar() {
	FILE* ht;

	DEBUG_START_FUNC;
	
	if(rank==0) {
		ht=fopen("timevar","w");
		fclose(ht);
#ifdef WITH_PARTICLES
		ht=fopen("partvar","w");
		fclose(ht);
#endif
	}

	DEBUG_END_FUNC;
	return;
}


