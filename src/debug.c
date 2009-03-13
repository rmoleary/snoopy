#include <stdlib.h>

#include "common.h"
#include "gfft.h"

#ifdef DEBUG
	
void D_show_field(PRECISION complex * field) {
	// Print several informations about the field field
	PRECISION complex * df;
	PRECISION * dfr;
	PRECISION maxfield, minfield, avgfield, avg2field;
	int i;
	
	df = (PRECISION complex *) fftw_malloc( sizeof(PRECISION complex) * NTOTAL_COMPLEX );
	if (df == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for df allocation");
	dfr = (PRECISION *) df;
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		df[i] =  field[i];
	}
	
	gfft_c2r(df);
	
	maxfield=dfr[0];
	minfield=dfr[0];
	avgfield=0.0;
	avg2field=0.0;
	
	for( i = 0 ; i < NTOTAL_COMPLEX * 2 ; i++) {
		if( dfr[i] > maxfield ) maxfield = dfr[i];
		if( dfr[i] < minfield ) minfield = dfr[i];
		avgfield+=dfr[i];
		avg2field+=dfr[i]*dfr[i];
	}
	
	maxfield=maxfield/ ((PRECISION) NTOTAL);
	minfield=minfield/ ((PRECISION) NTOTAL);
	avgfield=avgfield/ ((PRECISION)  NTOTAL*NTOTAL);
	avg2field=avg2field/ ((PRECISION) NTOTAL*NTOTAL*NTOTAL);
	
#ifdef MPI_SUPPORT
	reduce(&maxfield,2);
	reduce(&minfield,3);
	reduce(&avgfield,1);
	reduce(&avgfield,1);
#endif

	MPI_Printf("maxfield= %12e, minfield= %12e, avgfield= %12e, avg2field= %12e\n",maxfield, minfield, avgfield, avg2field);
	
	fftw_free(df);
}

void D_show_all(struct Field fldi) {
	MPI_Printf("   vx:");
	D_show_field(fldi.vx);
	MPI_Printf("   vy:");
	D_show_field(fldi.vy);
	MPI_Printf("   vz:");
	D_show_field(fldi.vz);
#ifdef BOUSSINESQ
	MPI_Printf("   th:");
	D_show_field(fldi.th);
#endif
#ifdef MHD
	MPI_Printf("   bx:");
	D_show_field(fldi.vx);
	MPI_Printf("   by:");
	D_show_field(fldi.vy);
	MPI_Printf("   bz:");
	D_show_field(fldi.vz);
#endif
	return;
}

#endif
	
	