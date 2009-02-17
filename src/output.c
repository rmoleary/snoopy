#include <stdlib.h>

#include "common.h"
#include "timestep.h"
#include "gfft.h"
#include "vtk_writer.h"

#define OUTPUT_SPECTRUM_N_BIN		64
#define	OUTPUT_SPECTRUM_FILENAME	"spectrum.dat"

#define	OUTPUT_DUMP					"dump.dmp"
#define	OUTPUT_DUMP_VERSION			02

#define	DUMP_MARKER					1981

int noutput_flow;
PRECISION lastoutput_time, lastoutput_flow, lastoutput_spectrum, lastoutput_dump;

void remap_output(	PRECISION wri[], 
					const PRECISION t) {
					
	int i,j,k;
	PRECISION tvelocity;
	PRECISION tremap;
	complex PRECISION wexp;
	complex PRECISION phase;
	
	
	tvelocity = fmod(t, 2.0 * LX / (SHEAR * LY));
	tremap = fmod(t + LY / (2.0 * SHEAR * LX) , LY / (SHEAR * LX)) - LY / (2.0 * SHEAR * LX);
	
	for( i = 0 ; i < NX ; i++) {
		for( k = 0 ; k < NZ ; k++) {
			for( j = 0 ; j < NY ; j++) {
				w1d[j] = wri[k + j * (NZ + 2) + NY * (NZ + 2) * i];
			}
		
			// Transform w1d, which will be stored in w2d
			fftw_execute(fft_1d_forward);
					
			for( j = 0 ; j < NY ; j++) {
			// advection
				phase = (PRECISION complex) ((2.0 * M_PI) / LY * (fmod( j + (NY / 2) ,  NY ) - NY / 2 ) * 
											( ((double) i / (double) NX ) * tremap - tvelocity / 2.0 ) * LX );
				wexp = cexp( I * phase);
									
				w2d[ j ] = w2d[ j ] * wexp;
			}
			
			fftw_execute(fft_1d_backward);
		
			for( j = 0 ; j < NY ; j++) {
				wri[k + j * (NZ + 2) + NY * (NZ + 2) * i] = w1d[j] / NY;
			}
		}
	}
	return;
}
		
void write_snap(const PRECISION t, const char filename[], const PRECISION complex wi[]) {
	// Write the complex field wi[] in file filename. We need t for the remap thing.
	FILE *ht;
	int i,j,k;
	
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i] = wi[i];
	}
	
	gfft_c2r(w1);
	
	for( i = 0 ; i < 2 * NTOTAL_COMPLEX ; i++) {
		wr1[i] = wr1[i] / ((double) NTOTAL );
	}
	
#ifdef WITH_SHEAR
	remap_output(wr1,t);
#endif
	
	ht=fopen(filename,"w");
#ifdef FORTRAN_OUTPUT_ORDER	
	for( j = 0; j < NY; j++) {
		for( i = 0; i < NX; i++) {
			for( k = 0 ; k < NZ; k++) {
#else
	for( i = 0; i < NX; i++) {
		for( j = 0; j < NY; j++) {
			for( k = 0 ; k < NZ; k++) {
#endif
				fwrite(&wr1[k + j * (NZ + 2) + i * NY * (NZ + 2) ], sizeof(PRECISION), 1, ht);
			}
		}
	}
	fclose(ht);
}

// VTK output using the visit_writer code.
void output_vtk(const int n, const PRECISION t) {
	int i,j,k;
	char  filename[50];
	char  varname1[10];
	char  varname2[10];
	char  varname3[10];
	char  varname4[10];
	
	char* varnames[4];
	
	float* vxf = (float *) wr1;
	float* vyf = (float *) wr2;
	float* vzf = (float *) wr3;
	float* thf = (float *) wr4;
	
	float xcoord[NX];
	float ycoord[NY];
	float zcoord[NZ];
	
	float* v[4];
	int dims[3];
	int vardims[4];
	int centering[4];
	
	// Init arrays required by vtk writer
	dims[0] = NX;
	dims[1] = NY;
	dims[2] = NZ;

	// 
	v[0] = vxf;
	v[1] = vyf;
	v[2] = vzf;
	v[3] = thf;

	vardims[0] = 1;
	vardims[1] = 1;
	vardims[2] = 1;
	vardims[3] = 1;
	
	centering[0] = 1;
	centering[1] = 1;
	centering[2] = 1;
	centering[3] = 1;

	// Init varnames
	sprintf(varname1,"vx");
	sprintf(varname2,"vy");
	sprintf(varname3,"vz");
	sprintf(varname4,"theta");
	
	varnames[0] = varname1;
	varnames[1] = varname2;
	varnames[2] = varname3;
	varnames[3] = varname4;
	
	// Init coordinate system
	for( i = 0 ; i < NX ; i++) {
		xcoord[i] = ((float) i) / ((float) NX) * LX - LX / 2.0;
	}
	for( i = 0 ; i < NY ; i++) {
		ycoord[i] = ((float) i) / ((float) NY) * LY - LY / 2.0;
	}
	for( i = 0 ; i < NZ ; i++) {
		zcoord[i] = ((float) i) / ((float) NZ) * LZ - LZ / 2.0;
	}
	// Put variables in the right format.
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w5[i] = fld.vx[i];
		w6[i] = fld.vy[i];
		w7[i] = fld.vz[i];
#ifdef BOUSSINESQ
		w8[i] = fld.th[i];
#endif
	}
	gfft_c2r(w5);
	gfft_c2r(w6);
	gfft_c2r(w7);
#ifdef BOUSSINESQ
	gfft_c2r(w8);
#endif

	for( i = 0 ; i < 2*NTOTAL_COMPLEX ; i++) {
		wr5[i] = wr5[i] / ((double) NTOTAL );
		wr6[i] = wr6[i] / ((double) NTOTAL );
		wr7[i] = wr7[i] / ((double) NTOTAL );
#ifdef BOUSSINESQ
		wr8[i] = wr8[i] / ((double) NTOTAL );
#endif
	}

#ifdef WITH_SHEAR
	remap_output(wr5,t);
	remap_output(wr6,t);
	remap_output(wr7,t);
#ifdef BOUSSINESQ
	remap_output(wr8,t);
#endif
#endif

	for( i = 0; i < NX; i++) {
		for( j = 0; j < NY; j++) {
			for( k = 0 ; k < NZ; k++) {
				vxf[i + j * NX + k * NX * NY] = (float) wr5[k + j * (NZ + 2) + i * (NZ + 2) * NY];
				vyf[i + j * NX + k * NX * NY] = (float) wr6[k + j * (NZ + 2) + i * (NZ + 2) * NY];
				vzf[i + j * NX + k * NX * NY] = (float) wr7[k + j * (NZ + 2) + i * (NZ + 2) * NY];
#ifdef BOUSSINESQ
				thf[i + j * NX + k * NX * NY] = (float) wr8[k + j * (NZ + 2) + i * (NZ + 2) * NY];
#endif
			}
		}
	}
	
	sprintf(filename,"data/v%04i.vtk",n);

	// Output everything
#ifdef BOUSSINESQ
	write_rectilinear_mesh(filename, 1, dims, xcoord, ycoord, zcoord, 4, vardims, centering, varnames, v);
#else
	write_rectilinear_mesh(filename, 1, dims, xcoord, ycoord, zcoord, 3, vardims, centering, varnames, v);
#endif
}
	
void output_flow(const int n, const PRECISION t) {
	
	char filename[50];
	
	sprintf(filename,"data/vx%04i.raw",n);
	write_snap(t, filename, fld.vx);
	
	sprintf(filename,"data/vy%04i.raw",n);
	write_snap(t, filename, fld.vy);
	
	sprintf(filename,"data/vz%04i.raw",n);
	write_snap(t, filename, fld.vz);
	
#ifdef BOUSSINESQ
	sprintf(filename,"data/th%04i.raw",n);
	write_snap(t, filename, fld.th);
#endif

	return;
}
	
void output_timevar(const struct Field fldi,
					const PRECISION t) {
					
	FILE *ht;
	PRECISION vxmax, vxmin, vymax, vymin, vzmax, vzmin, thmin, thmax;
	PRECISION energy_total;
	PRECISION transport;
	int i;
	
		
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i] = fldi.vx[i];
		w2[i] = fldi.vy[i];
		w3[i] = fldi.vz[i];
#ifdef BOUSSINESQ
		w4[i] = fldi.th[i];
#endif
	}

	energy_total = energy(w1) + energy(w2) + energy(w3);
	
	gfft_c2r(w1);
	gfft_c2r(w2);
	gfft_c2r(w3);
#ifdef BOUSSINESQ
	gfft_c2r(w4);
#endif

	for( i = 0 ; i < 2*NTOTAL_COMPLEX ; i++) {
		wr1[i] = wr1[i] / ((double) NTOTAL );
		wr2[i] = wr2[i] / ((double) NTOTAL );
		wr3[i] = wr3[i] / ((double) NTOTAL );
#ifdef BOUSSINESQ
		wr4[i] = wr4[i] / ((double) NTOTAL );
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
	transport=0.0;
		
	for(i = 0 ; i < 2 * NTOTAL_COMPLEX ; i ++) {
		if(vxmax < wr1[i]) vxmax = wr1[i];
		if(vxmin > wr1[i]) vxmin = wr1[i];
		if(vymax < wr2[i]) vymax = wr2[i];
		if(vymin > wr2[i]) vymin = wr2[i];
		if(vzmax < wr3[i]) vzmax = wr3[i];
		if(vzmin > wr3[i]) vzmin = wr3[i];
#ifdef BOUSSINESQ
		if(thmax < wr4[i]) thmax = wr4[i];
		if(thmin > wr4[i]) thmin = wr4[i];
#endif
		transport = transport + wr1[i] * wr2[i] / ((double) NTOTAL);

	}
		
	ht=fopen("timevar","a");
	fprintf(ht,"%08e\t",t);
	fprintf(ht,"%08e\t",energy_total);
	fprintf(ht,"%08e\t%08e\t%08e\t%08e\t%08e\t%08e\t",vxmax,vxmin,vymax,vymin,vzmax,vzmin);
#ifdef BOUSSINESQ
	fprintf(ht,"%08e\t%08e\t",thmax,thmin);
#else
	fprintf(ht,"%08e\t%08e\t",0.0,0.0);
#endif
	fprintf(ht,"%08e\n",transport);
	fclose(ht);
	return;
}
		
	
					
	
void output_dump( const struct Field fldi,
				  const PRECISION t) {	
				  
	FILE *ht;
	int dump_version;
	int size_x,	size_y, size_z;
	int marker;
	
	size_x = NX;
	size_y = NY;
	size_z = NZ;
	
	// This is a check when we try to read a dump file
	dump_version = OUTPUT_DUMP_VERSION;
	
	// This is a hard coded marker to check that we have read correctly the file
	marker = DUMP_MARKER;
	
	ht=fopen(OUTPUT_DUMP,"w");
	if(ht==NULL) ERROR_HANDLER( ERROR_CRITICAL, "Error opening dump file.");
	
	fwrite(&dump_version, sizeof(int), 1, ht);
	
	fwrite(&size_x		, sizeof(int), 1, ht);
	fwrite(&size_y		, sizeof(int), 1, ht);
	fwrite(&size_z		, sizeof(int), 1, ht);
	
	fwrite(fldi.vx			, sizeof(PRECISION complex), NTOTAL_COMPLEX, ht);
	fwrite(fldi.vy			, sizeof(PRECISION complex), NTOTAL_COMPLEX, ht);
	fwrite(fldi.vz			, sizeof(PRECISION complex), NTOTAL_COMPLEX, ht);
	
#ifdef BOUSSINESQ
	fwrite(fldi.th			, sizeof(PRECISION complex), NTOTAL_COMPLEX, ht);
#endif
	fwrite(&t			, sizeof(PRECISION)		   , 1			   , ht);
	
	fwrite(&noutput_flow		, sizeof(int)			   , 1             , ht);
	fwrite(&lastoutput_time 	, sizeof(PRECISION)		   , 1			   , ht);
	fwrite(&lastoutput_flow 	, sizeof(PRECISION)		   , 1			   , ht);
	fwrite(&lastoutput_spectrum	, sizeof(PRECISION)		   , 1			   , ht);
	fwrite(&lastoutput_dump 	, sizeof(PRECISION)		   , 1			   , ht);
	
// Any extra information should be put here.	
	
	fwrite(&marker		, sizeof(int)			   , 1			   , ht);
	
	fclose(ht);
	
	return;
}

void read_dump(   struct Field fldo,
				  PRECISION *t) {	
				  
	FILE *ht;
	int dump_version;
	int size_x,	size_y, size_z;
	int marker;
		
	ht=fopen(OUTPUT_DUMP,"r");
	if(ht==NULL) ERROR_HANDLER( ERROR_CRITICAL, "Error opening dump file.");
	
	fread(&dump_version, sizeof(int), 1, ht);
	if( dump_version != OUTPUT_DUMP_VERSION) ERROR_HANDLER( ERROR_CRITICAL, "Incorrect dump file version.");
	
	fread(&size_x		, sizeof(int), 1, ht);
	fread(&size_y		, sizeof(int), 1, ht);
	fread(&size_z		, sizeof(int), 1, ht);
	
	if(size_x != NX) ERROR_HANDLER( ERROR_CRITICAL, "Incorrect X grid size in dump file.");
	if(size_y != NY) ERROR_HANDLER( ERROR_CRITICAL, "Incorrect Y grid size in dump file.");
	if(size_z != NZ) ERROR_HANDLER( ERROR_CRITICAL, "Incorrect Y grid size in dump file.");
	
	fread(fldo.vx	, sizeof(PRECISION complex), NTOTAL_COMPLEX, ht);
	fread(fldo.vy	, sizeof(PRECISION complex), NTOTAL_COMPLEX, ht);
	fread(fldo.vz	, sizeof(PRECISION complex), NTOTAL_COMPLEX, ht);
#ifdef BOUSSINESQ
	fread(fldo.th	, sizeof(PRECISION complex), NTOTAL_COMPLEX, ht);
#endif
	fread(t			, sizeof(PRECISION)		   , 1			   , ht);
	
	fread(&noutput_flow			, sizeof(int)			   , 1             , ht);
	fread(&lastoutput_time		, sizeof(PRECISION)		   , 1			   , ht);
	fread(&lastoutput_flow		, sizeof(PRECISION)		   , 1			   , ht);
	fread(&lastoutput_spectrum	, sizeof(PRECISION)		   , 1			   , ht);
	fread(&lastoutput_dump		, sizeof(PRECISION)		   , 1			   , ht);
	
	fread(&marker , sizeof(int)			   , 1, ht);
	
	if(marker != DUMP_MARKER) ERROR_HANDLER( ERROR_CRITICAL, "Incorrect marker. Probably an incorrect dump file!");	
	fclose(ht);
	
	printf("Restarting at t=%e...\n",*t);
	return;
}
	
void init_output() {
#ifndef RESTART
	noutput_flow=0;
	lastoutput_time = T_INITIAL - TOUTPUT_TIME;
	lastoutput_flow = T_INITIAL - TOUTPUT_FLOW;
	lastoutput_dump = T_INITIAL - TOUTPUT_DUMP;
	
#endif
	return;
}

void output(const PRECISION t) {
	// Very rough output function
	if( (t-lastoutput_time)>=TOUTPUT_TIME) {
		output_timevar(fld,t);
		lastoutput_time = lastoutput_time + TOUTPUT_TIME;
	}
	
	if( (t-lastoutput_flow)>=TOUTPUT_FLOW) {
		output_flow(noutput_flow,t);
		output_vtk(noutput_flow,t);
		noutput_flow++;
		lastoutput_flow = lastoutput_flow + TOUTPUT_FLOW;
	}
		
	if( (t-lastoutput_dump)>=TOUTPUT_DUMP) {
		lastoutput_dump=lastoutput_dump+TOUTPUT_DUMP;
		output_dump(fld,t);
	}
	return;
}

void output_status(FILE * iostream) {
	fprintf(iostream,"Next output in file n %d, at t=%e\n",noutput_flow, lastoutput_flow+TOUTPUT_FLOW);
	return;
}

void output_immediate(const PRECISION t) {
	// Very rough output function
	// Immediate output
	output_timevar(fld,t);
	output_flow(noutput_flow,t);
	noutput_flow++;
	output_dump(fld,t);
	return;
}

void dump_immediate(const PRECISION t) {
	output_dump(fld,t);
	return;
}

void clear_timevar() {
	FILE* ht;
	ht=fopen("timevar","w");
	fclose(ht);
}
	
	