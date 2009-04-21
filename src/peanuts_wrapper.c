// Wrapper to be used with peanuts.
#include "common.h"
#include "mainloop.h"
#include "output.h"
#include "gfft.h"

int SnoopyInitFlag = 0;

/*******************************************************************/
/** 
	Solve the DNS using the Snoopy code from t_start to t_end
	Snoopy will generate a serie of snapshots and timevar files
	according to gvars.h
	
	@param fldi: input field (initial condition) using the C structure
				 described in common.h (Field). This structure contains the field in spectral
				 representation, with the disposition used for kx, ky and kz (see common.c).
	
	@param t_start: start time of the simulation (presumably 0)
	
	@param t_end: end time of the simulation (whatever you want).
*/
/******************************************************************/

void SnoopySolveDNS(struct Field fldi, double t_start, double t_end) {
	int i;
	if(SnoopyInitFlag == 1 ) {
	
		init_output();		// Restart the output routine for a new integration
		for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
			fld.vx[i] = fldi.vx[i];
			fld.vy[i] = fldi.vy[i];
			fld.vz[i] = fldi.vz[i];
#ifdef BOUSSINESQ		
			fld.th[i] = fldi.th[i];
#endif
#ifdef MHD
			fld.bx[i] = fldi.bx[i];
			fld.by[i] = fldi.by[i];
			fld.bz[i] = fldi.bz[i];
#endif
		}
		
		mainloop(t_start, t_end);
		finish_output();
	
		for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
			fldi.vx[i] = fld.vx[i];
			fldi.vy[i] = fld.vy[i];
			fldi.vz[i] = fld.vz[i];
#ifdef BOUSSINESQ		
			fldi.th[i] = fld.th[i];
#endif
#ifdef MHD
			fldi.bx[i] = fld.bx[i];
			fldi.by[i] = fld.by[i];
			fldi.bz[i] = fld.bz[i];
#endif
		}
	}
	else {
		ERROR_HANDLER(ERROR_CRITICAL, "You have not initialized Snoopy properly. Please call SnoopyInitDNS().");
	}
	return;
}

/*******************************************************************/
/** 
	Initialize Snoopy so that SnoopySolveDNS can be called
*/
/******************************************************************/

void SnoopyInitDNS() {
	if(SnoopyInitFlag == 0 ) {
#ifdef MPI_SUPPORT
		MPI_Comm_rank(MPI_COMM_WORLD,&rank);
		MPI_Comm_size(MPI_COMM_WORLD,&NPROC);
		// Some consistancy check
		if(NX/NPROC < 1) ERROR_HANDLER( ERROR_CRITICAL, "NX should be a multiple of the number of process.");
		if(NY/NPROC < 1) ERROR_HANDLER( ERROR_CRITICAL, "NY should be a multiple of the number of process.");
#else
		rank=0;
#endif

		init_common();
		init_gfft();
		
		SnoopyInitFlag = 1;
	}
	else {
		ERROR_HANDLER(ERROR_WARNING, "SnoopyInitDNS has already been called, no need to call it several times!");
	}
		
	return;
}

/*******************************************************************/
/** 
	Free the ressources used by Snoopy. 
*/
/******************************************************************/
void SnoopyFreeDNS() {
	if(SnoopyInitFlag == 1 ) {
		finish_gfft();
		finish_common();
		
		SnoopyInitFlag = 0;
	}
	else {
		ERROR_HANDLER(ERROR_WARNING, "Snoopy is not initialized, no need to call SnoopyFreeDNS()!");
	}
	
	return;
}
	

