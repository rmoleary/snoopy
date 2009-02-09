#include <stdlib.h>

#include "common.h"
#include "timestep.h"

#define OUTPUT_SPECTRUM_N_BIN		64
#define	OUTPUT_SPECTRUM_FILENAME	"spectrum.dat"

#define	OUTPUT_DUMP					"dump.dmp"
#define	OUTPUT_DUMP_VERSION			01

int noutput_flow;
PRECISION lastoutput_time, lastoutput_flow, lastoutput_spectrum, lastoutput_dump;

void remap_output(	PRECISION wri[], 
					const PRECISION t) {
					
	int i,j;
	PRECISION tvelocity;
	PRECISION tremap;
	complex PRECISION wexp;
	complex PRECISION phase;
	
	
	tvelocity = fmod(t, 2.0 * LX / LY);
	tremap = fmod(t + LY / (2.0 * LX) , LY /  LX) - LY / (2.0 * LX);
	
	for( i = 0 ; i < NX ; i++) {
		for( j = 0 ; j < NY ; j++) {
			w1d[j] = wri[j + NY * i];
		}
		
		// Transform w1d, which will be stored in w2d
		fftw_execute(fft_1d_forward);
					
		for( j = 0 ; j < NY ; j++) {
		// advection
			phase = (PRECISION complex) ((2.0 * M_PI) / LY * (fmod( j + (NY / 2) ,  NY ) - NY / 2 ) * 
										( ((double) i / (double) NX ) * tremap - tvelocity / 2.0 ) * LX );
			wexp = cexp( - I * phase);
									
			w2d[ j ] = w2d[ j ] * wexp;
		}
		
		fftw_execute(fft_1d_backward);
		
		for( j = 0 ; j < NY ; j++) {
			wri[j + NY * i] = w1d[j] / NY;
		}
	}
	return;
}

void compute_shear( PRECISION dwro[], PRECISION wri[]) {
	int i,j;
	complex PRECISION phase;
	
	for( i = 0 ; i < NX ; i++) {
		for( j = 0 ; j < NY ; j++) {
			w1d[j] = wri[j + NY * i];
		}
		
		// Transform w1d, which will be stored in w2d
		fftw_execute(fft_1d_forward);
					
		for( j = 0 ; j < NY ; j++) {
		// advection
			phase = (PRECISION complex) ((2.0 * M_PI) / LY * (fmod( j + (NY / 2) ,  NY ) - NY / 2 ) * 
										( ((double) i / (double) NX ) - 0.5 ) * LX );
									
			w2d[ j ] = - I * phase * w2d [ j ];
		}
		
		fftw_execute(fft_1d_backward);
		
		for( j = 0 ; j < NY ; j++) {
			dwro[j + NY * i] = dwro[j + NY * i] + w1d[j] / NY;
		}
	}
	return;
}


void init1Dspectrum() {
	int i,j,m;
	FILE * ht;
	PRECISION spectrum[ OUTPUT_SPECTRUM_N_BIN ];
		
	ht = fopen(OUTPUT_SPECTRUM_FILENAME,"w");
	for( m=0; m < OUTPUT_SPECTRUM_N_BIN; m++) 
		fprintf(ht,"%08e\t",kmax/OUTPUT_SPECTRUM_N_BIN*m);
	
	fprintf(ht,"\n");
	
	for( i = 0; i < OUTPUT_SPECTRUM_N_BIN; i++ )
		spectrum[ i ] = 0.0;
		
	for( i = 0; i < NX_COMPLEX; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			m = (int) ceil( ( pow( k2t[ IDX2D ], 0.5 ) * OUTPUT_SPECTRUM_N_BIN ) / kmax );
			if ( m < OUTPUT_SPECTRUM_N_BIN)
				spectrum[ m ] = spectrum[ m ] + 1.0;
		}
	}
	
	for( i = 0; i < OUTPUT_SPECTRUM_N_BIN; i++) 
		fprintf(ht,"%08e\t", spectrum[i]);
	
	fprintf(ht,"\n");
	
	fclose(ht);
	
	return;
}


void output1Dspectrum(const struct Field fldi) {
					  
	int i,j,m;
	PRECISION spectrum[ OUTPUT_SPECTRUM_N_BIN ];
	FILE *ht;
	
	for( i = 0; i < OUTPUT_SPECTRUM_N_BIN; i++ )
		spectrum[ i ] = 0.0;
		
	for( i = 0; i < NX_COMPLEX; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			m = (int) ceil( ( pow( k2t[ IDX2D ], 0.5 ) * OUTPUT_SPECTRUM_N_BIN ) / kmax );
			if ( m < OUTPUT_SPECTRUM_N_BIN) {
				if( j == 0) 
					// j=0, we have all the modes.
					spectrum[ m ] = spectrum[ m ] + creal( 0.5 * fldi.wz[ IDX2D ] * conj( fldi.wz[ IDX2D ] ) ) / ((PRECISION) NTOTAL*NTOTAL);
				else
					// j>0, only half of the complex plane is represented.
					spectrum[ m ] = spectrum[ m ] + creal( fldi.wz[ IDX2D ] * conj( fldi.wz[ IDX2D ] ) ) / ((PRECISION) NTOTAL*NTOTAL);
			}
		}
	}
	
	ht = fopen(OUTPUT_SPECTRUM_FILENAME,"a");

	for( i = 0; i < OUTPUT_SPECTRUM_N_BIN; i++) 
		fprintf(ht,"%08e\t", spectrum[i]);
	
	
	fprintf(ht,"\n");
	
	fclose(ht);
	return;
}
					  
void output_flow(const int n, const PRECISION t) {
	int i,j;
	struct Field dw;
	FILE *ht;
	char filenamewz[50];
	
	// Init the working field with a working field variable (we should not use boussinesq, otherwise crash!)
	dw.wz = w9;
	dw.th = w10;

	timestep(dw, fld, t, 0.0);
	
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i] = fld.wz[i];
		w2[i] = dw.wz[i];
		w3[i] = - nu * k2t[i] * fld.wz[i];
		w4[i] = I * N2 * ky[i] * fld.th[i];;
		w5[i] = w2[i] + w3[i];
	}
	
	fftw_execute_dft_c2r( c2rfft, w1, wr1);
	fftw_execute_dft_c2r( c2rfft, w2, wr2);
	fftw_execute_dft_c2r( c2rfft, w3, wr3);
	fftw_execute_dft_c2r( c2rfft, w4, wr4);
	fftw_execute_dft_c2r( c2rfft, w5, wr5);
	
	for( i = 0 ; i < NTOTAL ; i++) {
		wr1[i] = wr1[i] / ((double) NTOTAL );
		wr2[i] = wr2[i] / ((double) NTOTAL );
		wr3[i] = wr3[i] / ((double) NTOTAL );
		wr4[i] = wr4[i] / ((double) NTOTAL );
		wr5[i] = wr5[i] / ((double) NTOTAL );
	}
	
	// Compute the shear term hidden in dw (actually due to the sheared frame)
	compute_shear(wr2, wr1);
	compute_shear(wr5, wr1);
	
	remap_output(wr1,t);
	
	sprintf(filenamewz,"data/wz%04i.raw",n);

	ht=fopen(filenamewz,"w");
#ifdef FORTRAN_OUTPUT_ORDER	
	for( j = 0; j < NY; j++) {
		for( i = 0; i < NX; i++) {
#else
	for( i = 0; i < NX; i++) {
		for( j = 0; j < NY; j++) {
#endif
			fwrite(&wr1[j + i * NY], sizeof(PRECISION), 1, ht);
		}
	}
	fclose(ht);

	remap_output(wr2,t);
	
	sprintf(filenamewz,"data/dw%04i.raw",n);

	ht=fopen(filenamewz,"w");
#ifdef FORTRAN_OUTPUT_ORDER	
	for( j = 0; j < NY; j++) {
		for( i = 0; i < NX; i++) {
#else
	for( i = 0; i < NX; i++) {
		for( j = 0; j < NY; j++) {
#endif
			fwrite(&wr2[j + i * NY], sizeof(PRECISION), 1, ht);
		}
	}
	fclose(ht);

	remap_output(wr3,t);
	
	sprintf(filenamewz,"data/nw%04i.raw",n);

	ht=fopen(filenamewz,"w");
#ifdef FORTRAN_OUTPUT_ORDER	
	for( j = 0; j < NY; j++) {
		for( i = 0; i < NX; i++) {
#else
	for( i = 0; i < NX; i++) {
		for( j = 0; j < NY; j++) {
#endif
			fwrite(&wr3[j + i * NY], sizeof(PRECISION), 1, ht);
		}
	}
	fclose(ht);
	
	remap_output(wr4,t);
	
	sprintf(filenamewz,"data/bc%04i.raw",n);

	ht=fopen(filenamewz,"w");
#ifdef FORTRAN_OUTPUT_ORDER	
	for( j = 0; j < NY; j++) {
		for( i = 0; i < NX; i++) {
#else
	for( i = 0; i < NX; i++) {
		for( j = 0; j < NY; j++) {
#endif
			fwrite(&wr4[j + i * NY], sizeof(PRECISION), 1, ht);
		}
	}
	fclose(ht);


    remap_output(wr5,t);
	
	sprintf(filenamewz,"data/tt%04i.raw",n);

	ht=fopen(filenamewz,"w");
#ifdef FORTRAN_OUTPUT_ORDER	
	for( j = 0; j < NY; j++) {
		for( i = 0; i < NX; i++) {
#else
	for( i = 0; i < NX; i++) {
		for( j = 0; j < NY; j++) {
#endif
			fwrite(&wr5[j + i * NY], sizeof(PRECISION), 1, ht);
		}
	}
	fclose(ht);


	return;
}
	
void output_timevar(const struct Field fldi,
					const PRECISION t) {
					
	FILE *ht;
	PRECISION vxmax, vxmin, vymax, vymin, wzmax, wzmin, thmin, thmax;
	PRECISION energy_total, enstrophy_total;
	PRECISION transport, barocline_source, vortensity_dissip;
	int i;
	
		
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i] = fldi.wz[i];
#ifdef BOUSSINESQ
		w2[i] = fldi.th[i];
#endif
		w3[i] =   I * ky[i] * ik2t[i] * fldi.wz[i];
		w4[i] = - I * kxt[i] * ik2t[i] * fldi.wz[i];
		
		w5[i] = - nu * k2t[i] * fldi.wz[i];
#ifdef BOUSSINESQ
		w6[i] = I * N2 * ky[i] * fldi.th[i];
#endif
	}
	
	enstrophy_total = energy(fldi.wz);
	energy_total = energy(w3) + energy(w4);
	
	fftw_execute_dft_c2r( c2rfft, w1, wr1);
#ifdef BOUSSINESQ
	fftw_execute_dft_c2r( c2rfft, w2, wr2);
#endif
	fftw_execute_dft_c2r( c2rfft, w3, wr3);
	fftw_execute_dft_c2r( c2rfft, w4, wr4);
	fftw_execute_dft_c2r( c2rfft, w5, wr5);
#ifdef BOUSSINESQ	
	fftw_execute_dft_c2r( c2rfft, w6, wr6);
#endif

	for( i = 0 ; i < NTOTAL ; i++) {
		wr1[i] = wr1[i] / ((double) NTOTAL );
#ifdef BOUSSINESQ
		wr2[i] = wr2[i] / ((double) NTOTAL );
#endif
		wr3[i] = wr3[i] / ((double) NTOTAL );
		wr4[i] = wr4[i] / ((double) NTOTAL );
		
		wr5[i] = wr5[i] / ((double) NTOTAL );
#ifdef BOUSSINESQ
		wr6[i] = wr6[i] / ((double) NTOTAL );
#endif		
		
	}
	
	// w1 contains the vorticity
	// w2 contains the entropy perturbation (if appliable)
	// w3 contains vx
	// w4 contains vy
	// w5 contains nu * Delta w
	// w6 contains the baroclinic term
	
	vxmax=0.0;
	vxmin=0.0;
	vymax=0.0;
	vymin=0.0;
	wzmax=0.0;
	wzmin=0.0;
	thmax=0.0;
	thmin=0.0;
	transport=0.0;
	barocline_source=0.0;
	vortensity_dissip=0.0;
	
	for(i = 0 ; i < NTOTAL ; i ++) {
		if(vxmax < wr3[i]) vxmax = wr3[i];
		if(vxmin > wr3[i]) vxmin = wr3[i];
		if(vymax < wr4[i]) vymax = wr4[i];
		if(vymin > wr4[i]) vymin = wr4[i];
		if(wzmax < wr1[i]) wzmax = wr1[i];
		if(wzmin > wr1[i]) wzmin = wr1[i];
#ifdef BOUSSINESQ
		if(thmax < wr2[i]) thmax = wr2[i];
		if(thmin > wr2[i]) thmin = wr2[i];
#endif
		
		transport = transport + wr3[i] * wr4[i] / ((double) NTOTAL);
		vortensity_dissip = vortensity_dissip + wr1[i] * wr5[i] / ((double) NTOTAL);

#ifdef BOUSSINESQ
		barocline_source = barocline_source + wr1[i] * wr6[i] / ((double) NTOTAL);
#endif
	}
		
	ht=fopen("timevar","a");
	fprintf(ht,"%08e\t",t);
	fprintf(ht,"%08e\t%08e\t",energy_total,enstrophy_total);
	fprintf(ht,"%08e\t%08e\t%08e\t%08e\t%08e\t%08e\t%08e\t%08e\t",vxmax,vxmin,vymax,vymin,wzmax,wzmin,thmax,thmin);
	fprintf(ht,"%08e\t%08e\t%08e\n",transport,barocline_source,vortensity_dissip);
	fclose(ht);
	return;
}
		
	
					
	
void output_dump( const struct Field fldi,
				  const PRECISION t) {	
				  
	FILE *ht;
	int dump_version;
	int size_x,	size_y;
	int marker;
	
	size_x = NX;
	size_y = NY;
	
	// This is a check when we try to read a dump file
	dump_version = OUTPUT_DUMP_VERSION;
	
	// This is a hard coded marker to check that we have read correctly the file
	marker = 90;
	
	ht=fopen(OUTPUT_DUMP,"w");
	if(ht==NULL) ERROR_HANDLER( ERROR_CRITICAL, "Error opening dump file.");
	
	fwrite(&dump_version, sizeof(int), 1, ht);
	
	fwrite(&size_x		, sizeof(int), 1, ht);
	fwrite(&size_y		, sizeof(int), 1, ht);
	
	fwrite(fldi.wz			, sizeof(PRECISION complex), NTOTAL_COMPLEX, ht);
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
	int size_x,	size_y;
	int marker;
		
	ht=fopen(OUTPUT_DUMP,"r");
	if(ht==NULL) ERROR_HANDLER( ERROR_CRITICAL, "Error opening dump file.");
	
	fread(&dump_version, sizeof(int), 1, ht);
	if( dump_version != OUTPUT_DUMP_VERSION) ERROR_HANDLER( ERROR_CRITICAL, "Incorrect dump file version.");
	
	fread(&size_x		, sizeof(int), 1, ht);
	fread(&size_y		, sizeof(int), 1, ht);
	
	if(size_x != NX) ERROR_HANDLER( ERROR_CRITICAL, "Incorrect X grid size in dump file.");
	if(size_y != NY) ERROR_HANDLER( ERROR_CRITICAL, "Incorrect Y grid size in dump file.");
	
	fread(fldo.wz	, sizeof(PRECISION complex), NTOTAL_COMPLEX, ht);
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
	
	if(marker != 90) ERROR_HANDLER( ERROR_CRITICAL, "Incorrect marker. Probably an incorrect dump file!");	
	fclose(ht);
	
	printf("Restarting at t=%e...\n",*t);
	return;
}
	
void init_output() {
#ifndef RESTART
	noutput_flow=0;
	lastoutput_time = T_INITIAL - TOUTPUT_TIME;
	lastoutput_flow = T_INITIAL - TOUTPUT_FLOW;
	lastoutput_spectrum = T_INITIAL - TOUTPUT_SPECTR;
	lastoutput_dump = T_INITIAL - TOUTPUT_DUMP;
	
	// Init the spectrum file
	init1Dspectrum();
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
		noutput_flow++;
		lastoutput_flow = lastoutput_flow + TOUTPUT_FLOW;
	}
	
	if( (t-lastoutput_spectrum)>=TOUTPUT_SPECTR) {	
		output1Dspectrum(fld);
		lastoutput_spectrum=lastoutput_spectrum + TOUTPUT_SPECTR;
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
	output1Dspectrum(fld);
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
	
	