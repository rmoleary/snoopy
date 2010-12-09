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
#include "gfft.h"

#ifdef DEBUG
	
void D_show_field(double complex * field) {
	// Print several informations about the field field
	double complex * df;
	double * dfr;
	double maxfield, minfield, avgfield, avg2field;
	int i;
	
	df = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX );
	if (df == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for df allocation");
	dfr = (double *) df;
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
	
	maxfield=maxfield/ ((double) NTOTAL);
	minfield=minfield/ ((double) NTOTAL);
	avgfield=avgfield/ ((double)  NTOTAL*NTOTAL);
	avg2field=avg2field/ ((double) NTOTAL*NTOTAL*NTOTAL);
	
#ifdef MPI_SUPPORT
	reduce(&maxfield,2);
	reduce(&minfield,3);
	reduce(&avgfield,1);
	reduce(&avgfield,1);
#endif

	MPI_Printf("maxfield= %12e, minfield= %12e, avgfield= %12e, avg2field= %12e\n",maxfield, minfield, avgfield, avg2field);
	
	fftw_free(df);
}

#ifdef WITH_PARTICLES
void D_show_part(struct Field fldi) {
	int i;
	double xmax, xmin, ymax, ymin, zmax, zmin, vxmax, vxmin, vymax, vymin, vzmax, vzmin;
	double vxmean, vymean, vzmean;
	xmax = fldi.part[0].x;
	xmin = fldi.part[0].x;
	ymax = fldi.part[0].y;
	ymin = fldi.part[0].y;
	zmax = fldi.part[0].z;
	zmin = fldi.part[0].z;
	vxmax = fldi.part[0].vx;
	vxmin = fldi.part[0].vx;
	vymax = fldi.part[0].vy;
	vymin = fldi.part[0].vy;
	vzmax = fldi.part[0].vz;
	vzmin = fldi.part[0].vz;
	vxmean = 0.0;
	vymean = 0.0;
	vzmean = 0.0;
	
	for( i = 0 ; i < param.particles_n/NPROC ; i++) {
		if(fldi.part[i].x > xmax) xmax = fldi.part[i].x;
		if(fldi.part[i].x < xmin) xmin = fldi.part[i].x;
		if(fldi.part[i].y > ymax) ymax = fldi.part[i].y;
		if(fldi.part[i].y < ymin) ymin = fldi.part[i].y;
		if(fldi.part[i].z > zmax) zmax = fldi.part[i].z;
		if(fldi.part[i].z < zmin) zmin = fldi.part[i].z;
		
		if(fldi.part[i].vx > vxmax) vxmax = fldi.part[i].vx;
		if(fldi.part[i].vx < vxmin) vxmin = fldi.part[i].vx;
		if(fldi.part[i].vy > vymax) vymax = fldi.part[i].vy;
		if(fldi.part[i].vy < vymin) vymin = fldi.part[i].vy;
		if(fldi.part[i].vz > vzmax) vzmax = fldi.part[i].vz;
		if(fldi.part[i].vz < vzmin) vzmin = fldi.part[i].vz;
			
		vxmean = vxmean + fldi.part[i].vx;
		vymean = vymean + fldi.part[i].vy;
		vzmean = vzmean + fldi.part[i].vz;
	}
	
	vxmean = vxmean / param.particles_n;
	vymean = vymean / param.particles_n;
	vzmean = vzmean / param.particles_n;
	
	reduce(&xmax, 2);
	reduce(&xmin, 3);
	reduce(&ymax, 2);
	reduce(&ymin, 3);
	reduce(&zmax, 2);
	reduce(&zmin, 3);
	
	reduce(&vxmax, 2);
	reduce(&vxmin, 3);
	reduce(&vymax, 2);
	reduce(&vymin, 3);
	reduce(&vzmax, 2);
	reduce(&vzmin, 3);
	
	reduce(&vxmean, 1);
	reduce(&vymean, 1);
	reduce(&vzmean, 1);
	
	MPI_Printf("xmax = %12e, xmin = %12e, ymax = %12e, ymin = %12e, zmax = %12e, zmin = %12e\n",xmax, xmin, ymax, ymin, zmax, zmin);
	MPI_Printf("vxmax= %12e, vxmin= %12e, vymax= %12e, vymin= %12e, vzmax= %12e, vzmin= %12e\n", vxmax, vxmin, vymax, vymin, vzmax, vzmin);
	MPI_Printf("vxm  = %12e, vym  = %12e, vzm  = %12e\n",vxmean, vymean, vzmean);
	
	CHECK_NAN(xmax);
	CHECK_NAN(xmin);
	CHECK_NAN(ymax);
	CHECK_NAN(ymin);
	CHECK_NAN(zmax);
	CHECK_NAN(zmin);
	
	CHECK_NAN(vxmax);
	CHECK_NAN(vxmin);
	CHECK_NAN(vymax);
	CHECK_NAN(vymin);
	CHECK_NAN(vzmax);
	CHECK_NAN(vzmin);

	CHECK_NAN(vxmean);
	CHECK_NAN(vymean);
	CHECK_NAN(vzmean);
	
	return;
}
#endif

void D_show_all(struct Field fldi) {
	int i;
	for( i = 0 ; i < fldi.nfield ; i++) {
		MPI_Printf("   %s:",fldi.fname[i]);
		D_show_field(fldi.farray[i]);
	}
#ifdef WITH_PARTICLES
	MPI_Printf("   Particles:\n");
	D_show_part(fldi);
#endif
	return;
}

void debug_start_f(const char ErrorRoutine[], const int line, const char ErrorFile[]) {
	MPI_Printf("Start of %s (Line %d of file %s)\n", ErrorRoutine, line, ErrorFile);
	return;
}

void debug_end_f(const char ErrorRoutine[], const int line, const char ErrorFile[]) {
	MPI_Printf("End of %s (Line %d of file %s)\n", ErrorRoutine, line, ErrorFile);
	return;
}

#endif
	
	