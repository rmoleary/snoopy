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
#include "gfft.h"
#include "shear.h"
#include "debug.h"

#define MAX_N_BIN					10000
#define OUTPUT_SPECTRUM_K_BIN		(2.0 * M_PI)
#define	OUTPUT_SPECTRUM_FILENAME	"spectrum.dat"

#define	OUTPUT_DUMP					"dump.dmp"			/**< Dump files filename. */
#define OUTPUT_DUMP_SAV				"dump_sav.dmp"      /**< Previous (saved) output dump. */
#define OUTPUT_DUMP_WRITE			"dump_write.dmp"	/**< dump currently written. */

#define	OUTPUT_DUMP_VERSION			04					/**< Version of the dump files read and written by this code. */

#define	DUMP_MARKER					1981				/**< Marker used to signify the end of a dump file (it is also an excellent year...).*/

int noutput_flow;										/**< Next snapshot output number */
double lastoutput_time;								/**< Time when the last timevar output was done */
double lastoutput_flow;								/**< Time when the las snapshot output was done */
double lastoutput_dump;								/**< Time when the last dump output was done */

#ifdef WITH_SHEAR
double complex		*w1d, *w2d;						/** 1D arrays used by the remap methods */

fftw_plan	fft_1d_forward;								/**< 1D FFT transforms. Used by remap routines.*/
fftw_plan	fft_1d_backward;							/**< 1D FFT transforms. Used by remap routines.*/
#endif

#ifdef WITH_SHEAR
/***********************************************************/
/** 
	Remap a real field from the current sheared frame to the classical 
	cartesian frame. It could be more optimized, but since it is used
	in an output routine, I don't think such an optimization is worth doing.
	
	@param wri Real array to be remapped. The remapped array will be stored here.
	@param t Current time of the simulation (used to compute what the sheared frame is).
	
*/
/***********************************************************/

void remap_output(	double wri[], 
					const double t) {
					
	int i,j,k;
	double tvelocity;
	double tremap;
	complex double wexp;
	complex double phase;
	
	DEBUG_START_FUNC;
	
#ifdef TIME_DEPENDANT_SHEAR
	tremap = time_shift(t);
	tvelocity = 0.0;
#else	
	tremap = time_shift(t);
	tvelocity = fmod(t, 2.0 * param.ly / (param.shear * param.lx));
#endif	
	
	
	for( i = 0 ; i < NX / NPROC ; i++) {
		for( k = 0 ; k < NZ ; k++) {
			for( j = 0 ; j < NY ; j++) {
#ifdef WITH_2D
				w1d[j] = wri[j + (NY + 2) * i];
#else
				w1d[j] = wri[k + j * (NZ + 2) + NY * (NZ + 2) * i];
#endif
			}
		
			// Transform w1d, which will be stored in w2d
			fftw_execute(fft_1d_forward);
					
			for( j = 0 ; j < NY ; j++) {
			// advection phase = ky*
#ifdef TIME_DEPENDANT_SHEAR
				// There is no proper remap in this case
				phase = (double complex) ((2.0 * M_PI) / param.ly * (fmod( j + (NY / 2) ,  NY ) - NY / 2 ) * 
											( ((double) (i + rank * (NX/NPROC)) / (double) NX - 0.5 ) * tremap ) * param.lx * param.shear );
#else
				phase = (double complex) ((2.0 * M_PI) / param.ly * (fmod( j + (NY / 2) ,  NY ) - NY / 2 ) * 
											( ((double) (i + rank * (NX/NPROC)) / (double) NX ) * tremap - tvelocity / 2.0 ) * param.lx * param.shear);
#endif
				wexp = cexp( I * phase);
									
				w2d[ j ] = w2d[ j ] * wexp;
			}
			
			fftw_execute(fft_1d_backward);
		
			for( j = 0 ; j < NY ; j++) {
#ifdef WITH_2D
				wri[j + (NY + 2) * i] = w1d[j] / NY;
#else
				wri[k + j * (NZ + 2) + NY * (NZ + 2) * i] = w1d[j] / NY;
#endif
			}
		}
	}
	
	DEBUG_END_FUNC;
	
	return;
}


#endif

/***********************************************************/
/**
	compute a shell-integrated spectrum of the tensor real(wi*wj+)
	and write it on OUTPUT_SPECTRUM_FILENAME
	
	@param wi 1st double complex array from which the spectrum has to be deduced
	@param wj 2nd double complex array from which the spectrum has to be deduced
	@param ti Current time
*/
/***********************************************************/

void write_spectrum(const double complex wi[], const double complex wj[], const double ti) {
	DEBUG_START_FUNC;
	int i,j,k,m;
	int nbin;
	double spectrum[ MAX_N_BIN ];
	FILE *ht;
	
	nbin = (int) ceil(kmax / OUTPUT_SPECTRUM_K_BIN);
	
	for( i = 0; i < MAX_N_BIN; i++ )
		spectrum[ i ] = 0.0;
		
	for( i = 0; i < NX_COMPLEX/NPROC; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			for( k = 0; k < NZ_COMPLEX; k++) {
				m = (int) floor( pow( k2t[ IDX3D ], 0.5 ) / OUTPUT_SPECTRUM_K_BIN + 0.5 );
				if ( m < nbin) {
#ifdef WITH_2D
					if( j == 0)
#else
					if( k == 0) 
#endif
						// k=0, we have all the modes.
						spectrum[ m ] = spectrum[ m ] + creal( wi[ IDX3D ] * conj( wj[ IDX3D ] ) ) / ((double) NTOTAL*NTOTAL);
					else
						// k>0, only half of the complex plane is represented.
						spectrum[ m ] = spectrum[ m ] + creal( 2.0 * wi[ IDX3D ] * conj( wj[ IDX3D ] ) ) / ((double) NTOTAL*NTOTAL);
				}
			}
		}
	}
	
#ifdef MPI_SUPPORT
	// Reduce each component
	for( m=0; m < nbin; m++)
		reduce(&spectrum[m], 1);
#endif

	if(rank==0) {
		ht = fopen(OUTPUT_SPECTRUM_FILENAME,"a");
		fprintf(ht,"%08e\t", ti);
		for( i = 0; i < nbin; i++) 
			fprintf(ht,"%08e\t", spectrum[i]);
	
		fprintf(ht,"\n");
		fclose(ht);
	}

		
	DEBUG_END_FUNC;
	return;
}

/**********************************************************/
/**
	Output the transport spectrum in a file (OUTPUT_SPECTRUM_FILENAME)
	This routine is called only when shear is present.
	
	@param fldi: field from which the transport is computed
	@param ti: current time
*/
/*********************************************************/

void output1Dspectrum(const struct Field fldi, const double ti) {
	int i;
	
	DEBUG_START_FUNC;
	
	// The zero array is to be used for dummy spectrums
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i] = 0.0;
	}
	
	// V,B and theta spectrums
	write_spectrum(fldi.vx, fldi.vx, ti);
	write_spectrum(fldi.vy, fldi.vy, ti);
	write_spectrum(fldi.vz, fldi.vz, ti);
	
#ifdef MHD
	write_spectrum(fldi.bx, fldi.bx, ti);
	write_spectrum(fldi.by, fldi.by, ti);
	write_spectrum(fldi.bz, fldi.bz, ti);
#else
	write_spectrum(w1, w1, ti);
	write_spectrum(w1, w1, ti);
	write_spectrum(w1, w1, ti);
#endif

#ifdef BOUSSINESQ
	write_spectrum(fldi.th, fldi.th, ti);
#else
	write_spectrum(w1, w1, ti);
#endif

	// Transport spectrums
	write_spectrum(fldi.vx,fldi.vy, ti);
#ifdef MHD
	write_spectrum(fldi.bx,fldi.by, ti);
#else
	write_spectrum(w1, w1, ti);
#endif
	
	// Transfer spectrums
	// Kinetic energy transfer
	
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i] =  fldi.vx[i];
		w2[i] =  fldi.vy[i];
		w3[i] =  fldi.vz[i];
	}

	gfft_c2r_t(w1);
	gfft_c2r_t(w2);
	gfft_c2r_t(w3);
	
		/* Compute the convolution for the advection process */
	for( i = 0 ; i < 2*NTOTAL_COMPLEX ; i++) {
		wr4[i] = wr1[i] * wr1[i] / ((double) NTOTAL*NTOTAL);
		wr5[i] = wr2[i] * wr2[i] / ((double) NTOTAL*NTOTAL);
		wr6[i] = wr3[i] * wr3[i] / ((double) NTOTAL*NTOTAL);
		wr7[i] = wr1[i] * wr2[i] / ((double) NTOTAL*NTOTAL);
		wr8[i] = wr1[i] * wr3[i] / ((double) NTOTAL*NTOTAL);
		wr9[i] = wr2[i] * wr3[i] / ((double) NTOTAL*NTOTAL);
	}
	
	gfft_r2c_t(wr4);
	gfft_r2c_t(wr5);
	gfft_r2c_t(wr6);
	gfft_r2c_t(wr7);
	gfft_r2c_t(wr8);
	gfft_r2c_t(wr9);

	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i] = - I * mask[i] * (
					kxt[i] * w4[i] + ky[i] * w7[i] + kz[i] * w8[i] );
		w2[i] = - I * mask[i] * (
					kxt[i] * w7[i] + ky[i] * w5[i] + kz[i] * w9[i] );
		w3[i] = - I * mask[i] * (
					kxt[i] * w8[i] + ky[i] * w9[i] + kz[i] * w6[i] );
	}
	
	write_spectrum(fldi.vx, w1, ti);
	write_spectrum(fldi.vy, w2, ti);
	write_spectrum(fldi.vz, w3, ti);
	
#ifdef MHD

	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i] =  fldi.vx[i];
		w2[i] =  fldi.vy[i];
		w3[i] =  fldi.vz[i];
		w4[i] =  fldi.bx[i];
		w5[i] =  fldi.by[i];
		w6[i] =  fldi.bz[i];
	}

	gfft_c2r_t(w1);
	gfft_c2r_t(w2);
	gfft_c2r_t(w3);
	gfft_c2r_t(w4);
	gfft_c2r_t(w5);
	gfft_c2r_t(w6);
	
	// (vx,vy,vz) is in w1-w3 and (bx,by,bz) is in (w4-w6). It is now time to compute the emfs in w7-w9...

	for( i = 0 ; i < 2*NTOTAL_COMPLEX ; i++) {
		wr7[i] = (wr2[i] * wr6[i] - wr3[i] * wr5[i]) / ((double) NTOTAL*NTOTAL);
		wr8[i] = (wr3[i] * wr4[i] - wr1[i] * wr6[i]) / ((double) NTOTAL*NTOTAL);
		wr9[i] = (wr1[i] * wr5[i] - wr2[i] * wr4[i]) / ((double) NTOTAL*NTOTAL);
	}

	// Compute the curl of the emf involved in the induction equation.
	
	gfft_r2c_t(wr7);
	gfft_r2c_t(wr8);
	gfft_r2c_t(wr9);
	

	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i] = I * mask[i] * (ky[i] * w9[i] - kz[i] * w8[i]);
		w2[i] = I * mask[i] * (kz[i] * w7[i] - kxt[i]* w9[i]);
		w3[i] = I * mask[i] * (kxt[i]* w8[i] - ky[i] * w7[i]);
	}

	write_spectrum(fldi.bx, w1, ti);
	write_spectrum(fldi.by, w2, ti);
	write_spectrum(fldi.bz, w3, ti);
	
	// Let's do the Lorentz Force
	// We already have (bx,by,bz) in w4-w6. No need to compute them again...

	for( i = 0 ; i < 2*NTOTAL_COMPLEX ; i++) {
		wr1[i] = wr4[i] * wr4[i] / ((double) NTOTAL*NTOTAL);
		wr2[i] = wr5[i] * wr5[i] / ((double) NTOTAL*NTOTAL);
		wr3[i] = wr6[i] * wr6[i] / ((double) NTOTAL*NTOTAL);
		wr7[i] = wr4[i] * wr5[i] / ((double) NTOTAL*NTOTAL);
		wr8[i] = wr4[i] * wr6[i] / ((double) NTOTAL*NTOTAL);
		wr9[i] = wr5[i] * wr6[i] / ((double) NTOTAL*NTOTAL);
	}

	gfft_r2c_t(wr1);
	gfft_r2c_t(wr2);
	gfft_r2c_t(wr3);
	gfft_r2c_t(wr7);
	gfft_r2c_t(wr8);
	gfft_r2c_t(wr9);


	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w4[i] = I * mask[i] * (kxt[i] * w1[i] + ky[i] * w7[i] + kz[i] * w8[i]);
		w5[i] = I * mask[i] * (kxt[i] * w7[i] + ky[i] * w2[i] + kz[i] * w9[i]);
		w6[i] = I * mask[i] * (kxt[i] * w8[i] + ky[i] * w9[i] + kz[i] * w3[i]);
	}
	
	write_spectrum(fldi.vx, w4, ti);
	write_spectrum(fldi.vy, w5, ti);
	write_spectrum(fldi.vz, w6, ti);
	
	// Helicity spectrums
	
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i] = I * ik2t[i] * (ky[i] * fldi.bz[i] - kz[i] * fldi.by[i] );
		w2[i] = I * ik2t[i] * (kz[i] * fldi.bx[i] - kxt[i]* fldi.bz[i] );
		w3[i] = I * ik2t[i] * (kxt[i]* fldi.by[i] - ky[i] * fldi.bx[i] );
	}
	
	write_spectrum(fldi.bx, w1, ti);
	write_spectrum(fldi.by, w2, ti);
	write_spectrum(fldi.bz, w3, ti);
	
#else
	// The zero array is to be used for dummy spectrums
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i] = 0.0;
	}
	
	write_spectrum(w1, w1, ti);
	write_spectrum(w1, w1, ti);
	write_spectrum(w1, w1, ti);
	write_spectrum(w1, w1, ti);
	write_spectrum(w1, w1, ti);
	write_spectrum(w1, w1, ti);
	write_spectrum(w1, w1, ti);
	write_spectrum(w1, w1, ti);
	write_spectrum(w1, w1, ti);
#endif

	
	DEBUG_END_FUNC;
	
	return;
}

/**********************************************************/
/**
	Initialise the 1D spectrum output routine, used
	to output the spectrum
	This routine print the mode ks in the first line
	It also counts the number of mode in each shell and 
	output it in the second line of OUTPUT_SPECTRUM_FILENAME
*/
/*********************************************************/
void init1Dspectrum() {
	int i,j,k,m;
	int nbin;
	FILE * ht;
	double spectrum[ MAX_N_BIN ];
	
	DEBUG_START_FUNC;
	
	nbin = (int) ceil(kmax / OUTPUT_SPECTRUM_K_BIN);
	
	if(rank==0) {
		ht = fopen(OUTPUT_SPECTRUM_FILENAME,"w");
		
		for( m=0; m < nbin; m++) 
			fprintf(ht,"%08e\t", m * OUTPUT_SPECTRUM_K_BIN);
	
		fprintf(ht,"\n");
	}
	
	for( i = 0; i < MAX_N_BIN ; i++ )
		spectrum[ i ] = 0.0;
		
	for( i = 0; i < NX_COMPLEX/NPROC ; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			for( k = 0; k < NZ_COMPLEX; k++) {
				m = (int) floor( pow( k2t[ IDX3D ], 0.5 ) / OUTPUT_SPECTRUM_K_BIN + 0.5 );
				if ( m < nbin)
					spectrum[ m ] = spectrum[ m ] + 1.0;
			}
		}
	}
	
#ifdef MPI_SUPPORT
	// Reduce each component
	for( m=0; m < nbin ; m++)
		reduce(&spectrum[m], 1);
#endif

	if(rank==0) {
		for( i = 0; i < nbin ; i++) 
			fprintf(ht,"%08e\t", spectrum[i]);
	
		fprintf(ht,"\n");
	
		fclose(ht);
	}
	
	DEBUG_END_FUNC;
	
	return;
}



/*************************************************
** VTK For HD*************************************
**************************************************/

/* ****************************************************************************/
/** Determines if the machine is little-endian.  If so, 
    it will force the data to be big-endian. 
	@param in_number floating point number to be converted in big endian */
/* *************************************************************************** */

float big_endian(float in_number)
{
    static int doneTest = 0;
    static int shouldSwap = 0;
	
    if (!doneTest)
    {
        int tmp1 = 1;
        unsigned char *tmp2 = (unsigned char *) &tmp1;
        if (*tmp2 != 0)
            shouldSwap = 1;
        doneTest = 1;
    }

    if (shouldSwap)
    {
		unsigned char *bytes = (unsigned char*) &in_number;
        unsigned char tmp = bytes[0];
        bytes[0] = bytes[3];
        bytes[3] = tmp;
        tmp = bytes[1];
        bytes[1] = bytes[2];
        bytes[2] = tmp;
    }
	return(in_number);
}

/***********************************************************/
/** 
	Output a floating point big endian array from a complex array.
	Use Forran output format (required by VTK). This routine is essentially
	useful for VTK files (hence its name...).
 	
	@param ht File handler in which the data has to be written
	@param wi Complex array containing the field to be written. This array is transformed into real space
	and remapped (if SHEAR is present) before being written.
	@param t Current time of the simulation.
	
*/
/***********************************************************/
void write_vtk(FILE * ht, double complex wi[], const double t) {
	// Write the data in the file handler *ht
	int i,j,k;
	float q0;

	DEBUG_START_FUNC;

#ifdef MPI_SUPPORT
	double * chunk = NULL;
	if(rank==0) {
		chunk = (double *) malloc( NX * sizeof(double));
		if (chunk == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for chunk allocation");
	}
#endif

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

#ifdef BOUNDARY_C
	for( k = 0 ; k < NZ / 2 ; k++) {
#else
	for( k = 0 ; k < NZ; k++) {
#endif
		for( j = 0; j < NY; j++) {
#ifdef MPI_SUPPORT
			// We have to transpose manually to be Fortran compliant
			for(i = 0; i < NX/NPROC; i++) {
				wr2[i] = wr1[k + j * (NZ + 2) + i * NY * (NZ + 2) ];   // Transfer the chunk of data to wr2
			}
			MPI_Gather(wr2, NX/NPROC, MPI_DOUBLE,
					   chunk, NX/NPROC, MPI_DOUBLE, 0, MPI_COMM_WORLD); // Put the full chunk in chunk in the root process
#endif			
			for( i = 0; i < NX; i++) {
#ifdef MPI_SUPPORT
				if(rank==0) {
					q0 = big_endian( (float) chunk[ i ] );
					fwrite(&q0, sizeof(float), 1, ht);
				}
#else
#ifdef WITH_2D
				q0 = big_endian( (float) wr1[j + i * (NY + 2)] );
#else
				q0 = big_endian( (float) wr1[k + j * (NZ + 2) + i * NY * (NZ + 2) ] );
#endif
				fwrite(&q0, sizeof(float), 1, ht);
#endif
			}
#ifdef MPI_SUPPORT
			MPI_Barrier(MPI_COMM_WORLD);
#endif
		}
	}

#ifdef MPI_SUPPORT	
	if(rank==0) free(chunk);
#endif

	DEBUG_END_FUNC;
	
	return;
}

// Geo's personnal VTK writer, using structured data points	
/***********************************************************/
/** 
	Output a legacy VTK file readable by Paraview. This routine
	will output all the variables in files data/v****.vtk.
	
	@param n Number of the file in which the output will done.
	@param t Current time of the simulation.
*/
/***********************************************************/

void output_vtk(const int n, double t) {
	FILE *ht = NULL;
	char  filename[50];
	int num_remain_field;
	int array_size, i;
	
	DEBUG_START_FUNC;

	sprintf(filename,"data/v%04i.vtk",n);
#ifdef BOUNDARY_C
	array_size=NX*NY*NZ/2;	// Remove half of the vertical direction for symmetry reason when using walls in z
#else
	array_size=NX*NY*NZ;
#endif

	if(rank==0) {
		ht=fopen(filename,"w");
	
		fprintf(ht, "# vtk DataFile Version 2.0\n");
		fprintf(ht, "t= %015.15e Snoopy Code v5.0\n",t);
		fprintf(ht, "BINARY\n");
		fprintf(ht, "DATASET STRUCTURED_POINTS\n");
#ifdef BOUNDARY_C
		fprintf(ht, "DIMENSIONS %d %d %d\n", NX, NY, NZ / 2);
#else
		fprintf(ht, "DIMENSIONS %d %d %d\n", NX, NY, NZ);
#endif
		fprintf(ht, "ORIGIN %g %g %g\n", -param.lx/2.0, -param.ly/2.0, -param.lz/2.0);
		fprintf(ht, "SPACING %g %g %g\n", param.lx/NX, param.ly/NY, param.lz/NZ);
	
		// Write the primary scalar (f***ing VTK legacy format...)
		fprintf(ht, "POINT_DATA %d\n",array_size);
		fprintf(ht, "SCALARS %s float\n",fld.fname[0]);
		fprintf(ht, "LOOKUP_TABLE default\n");
	}
	write_vtk(ht,fld.farray[0],t);
	
	num_remain_field = fld.nfield - 1;		// we have already written the first one
	
	if(param.output_pressure)
		num_remain_field +=1;
		
	if(param.output_vorticity)
		num_remain_field +=3;
		
	if(rank==0) fprintf(ht, "FIELD FieldData %d\n",num_remain_field);
	
	// Write all the remaining fields
	
	for(i = 1 ; i < fld.nfield ; i++) {
		if(rank==0) fprintf(ht, "%s 1 %d float\n",fld.fname[i],array_size);
		write_vtk(ht,fld.farray[i],t);
	}
	
	if(param.output_pressure) {
		if(rank==0) fprintf(ht, "p 1 %d float\n",array_size);	// Output the pressure field when needed.
		write_vtk(ht,pressure,t);
	}
	if(param.output_vorticity) {
		// Compute the vorticity field
		for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
			w4[i] = I * (ky[i] * fld.vz[i] - kz[i] * fld.vy[i]);
			w5[i] = I * (kz[i] * fld.vx[i] - kxt[i] * fld.vz[i]);
			w6[i] = I * (kxt[i] * fld.vy[i] - ky[i] * fld.vx[i]);
		}
		if(rank==0) fprintf(ht, "wx 1 %d float\n",array_size);
		write_vtk(ht,w4,t);
		if(rank==0) fprintf(ht, "wy 1 %d float\n",array_size);
		write_vtk(ht,w5,t);
		if(rank==0) fprintf(ht, "wz 1 %d float\n",array_size);
		write_vtk(ht,w6,t);
	}
		
	if(rank==0) fclose(ht);
	
	DEBUG_END_FUNC;
	
	return;
	
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
	double vxmax, vxmin, vymax, vymin, vzmax, vzmin, thmin, thmax;
	double bxmax, bxmin, bymax, bymin, bzmax, bzmin;
	double vort2, curr2;
	double energy_total;	
	double energy_mag;
	double reynolds_stress;
	double maxwell_stress;
	double helicity;

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
	}

	energy_total = energy(w1) + energy(w2) + energy(w3);
	reduce(&energy_total,1);
#ifdef MHD
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
#ifdef MHD
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
#ifdef MHD
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
#ifdef MHD
		if(bxmax < wr5[i]) bxmax = wr5[i];
		if(bxmin > wr5[i]) bxmin = wr5[i];
		if(bymax < wr6[i]) bymax = wr6[i];
		if(bymin > wr6[i]) bymin = wr6[i];
		if(bzmax < wr7[i]) bzmax = wr7[i];
		if(bzmin > wr7[i]) bzmin = wr7[i];
		
		maxwell_stress += wr5[i] * wr6[i] / ((double) NTOTAL);
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
/**********************************************************
** Restart DUMP I/O routines ******************************
***********************************************************/

/***********************************************************/
/** 
	write a field information to a restart dump file, taking care of the MPI reduction
	@param handler handler to an already opened restart dump filed
	@param fldwrite Field structure pointing to the field needed to be saved
*/
/***********************************************************/
void write_field(FILE *handler, double complex *fldwrite) {
#ifdef MPI_SUPPORT
	MPI_Status status;
	// Write in rank order using the file opened if handler
	int current_rank;
#endif
	int i;

#ifdef MPI_SUPPORT	
	if(rank==0) {
		for(current_rank=0; current_rank < NPROC; current_rank++) {
			if(current_rank==0) {
				// Copy the dump in the rank=0 process
#endif
				for(i=0; i< NTOTAL_COMPLEX; i++) {
					w1[i]=fldwrite[i];
				}
#ifdef MPI_SUPPORT
			}
			else {
				MPI_Recv( w1, NTOTAL_COMPLEX * 2, MPI_DOUBLE, current_rank, 2, MPI_COMM_WORLD, &status);
			}
#endif
			fwrite(w1, sizeof(double complex), NTOTAL_COMPLEX, handler);
#ifdef MPI_SUPPORT
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
	else {
		MPI_Send(fldwrite, NTOTAL_COMPLEX * 2, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
	}
#endif
	return;
}

/**************************************************************************************/
/** 
	read a field information from a restart dump file, taking care of the MPI broadcast
	@param handler handler to an already opened restart dump filed
	@param fldread Field structure in which the restart information will be stored
*/
/**************************************************************************************/
	
void read_field(FILE *handler, double complex *fldread) {
#ifdef MPI_SUPPORT
	MPI_Status status;
	// Write in rank order using the file opened if handler
	int current_rank;
#endif
	int i;

#ifdef MPI_SUPPORT
	if(rank==0) {
		for(current_rank=0; current_rank < NPROC; current_rank++) {
#endif
			fread(w1, sizeof(double complex), NTOTAL_COMPLEX, handler);

#ifdef MPI_SUPPORT
			if(current_rank==0) {
#endif
				// Copy the dump in the rank=0 process
				for(i=0; i< NTOTAL_COMPLEX; i++) {
					fldread[i]=w1[i];
				}
#ifdef MPI_SUPPORT
			}
			else {
				MPI_Send( w1, NTOTAL_COMPLEX * 2, MPI_DOUBLE, current_rank, 3, MPI_COMM_WORLD);
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
	else {
		MPI_Recv(fldread, NTOTAL_COMPLEX * 2, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD,&status);
		MPI_Barrier(MPI_COMM_WORLD);
	}
#endif
	return;
}

/**************************************************************************************/
/** 
	Check if a file exists
	@param filename string containing the filename (including path) to be tested.
*/
/**************************************************************************************/
int file_exist(char filename[]) {
	FILE* ht;
	int file_status;
	
	ht=NULL;
	
	if(rank==0) {
		ht=fopen(filename,"r");
		if(ht==NULL) file_status = 0;
		else {
			file_status = 1;
			fclose(ht);
		}
	}
	
#ifdef MPI_SUPPORT
	MPI_Bcast(&file_status, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

	return(file_status);
}

/**************************************************************************************/
/** 
	write a full restart dump
	@param fldi Field structure pointing to the field needed to be savec
	@param t current time of the simulation (will be stored in the dump file)
*/
/**************************************************************************************/
void output_dump( const struct Field fldi,
				  const double t) {	
				  
	FILE *ht;
	int dump_version;
	int size_x,	size_y, size_z;
	int marker, included_field;
	
	ht=NULL;
	
	DEBUG_START_FUNC;
	
	size_x = NX;
	size_y = NY;
	size_z = NZ;
	
	// This is a check when we try to read a dump file
	dump_version = OUTPUT_DUMP_VERSION;
	
	// This is a hard coded marker to check that we have read correctly the file
	marker = DUMP_MARKER;
	
	if(rank==0) {
		ht=fopen(OUTPUT_DUMP_WRITE,"w");
		if(ht==NULL) ERROR_HANDLER( ERROR_CRITICAL, "Error opening dump file.");
	
		fwrite(&dump_version, sizeof(int), 1, ht);
	
		fwrite(&size_x		, sizeof(int), 1, ht);
		fwrite(&size_y		, sizeof(int), 1, ht);
		fwrite(&size_z		, sizeof(int), 1, ht);
		// Included fields
		// First bit is Boussinesq fields
		// Second bit is MHD fields
		// Other fields can be added from that stage...
		included_field=0;
#ifdef BOUSSINESQ
		included_field+=1;
#endif
#ifdef MHD
		included_field+=2;
#endif
		fwrite(&included_field, sizeof(int), 1, ht);
	}
	
	write_field(ht, fldi.vx);
	write_field(ht, fldi.vy);
	write_field(ht, fldi.vz);
	
#ifdef BOUSSINESQ
	write_field(ht, fldi.th);
#endif
#ifdef MHD
	write_field(ht, fldi.bx);
	write_field(ht, fldi.by);
	write_field(ht, fldi.bz);
#endif

	if(rank==0) {
		fwrite(&t			, sizeof(double)		   , 1			   , ht);
	
		fwrite(&noutput_flow		, sizeof(int)			   , 1             , ht);
		fwrite(&lastoutput_time 	, sizeof(double)		   , 1			   , ht);
		fwrite(&lastoutput_flow 	, sizeof(double)		   , 1			   , ht);
		fwrite(&lastoutput_dump 	, sizeof(double)		   , 1			   , ht);
	
// Any extra information should be put here.	
	
		fwrite(&marker		, sizeof(int)			   , 1			   , ht);
	
		fclose(ht);
	}
	
// This bit prevents the code from loosing all the dump files (this kind of thing happens sometimes...)
// With this routine, one will always have a valid restart dump, either in OUTPUT_DUMP_WRITE, OUTPUT_DUMP or OUTPUT_DUMP_SAV 
// (it should normally be in OUTPUT_DUMP)

	if(rank==0) {
		remove(OUTPUT_DUMP_SAV);				 // Delete the previously saved output dump
		rename(OUTPUT_DUMP, OUTPUT_DUMP_SAV);	 // Save the current dump file
		rename(OUTPUT_DUMP_WRITE, OUTPUT_DUMP);  // Move the new dump file to its final location
	}
	
#ifdef MPI_SUPPORT
	MPI_Barrier(MPI_COMM_WORLD);
#endif
	DEBUG_END_FUNC;
	
	return;
}

/**************************************************************************************/
/** 
	read a full restart dump
	@param fldo Field structure in which the restart dump will be stored.
	@param t time of the simulation, overwritten by this with the dump information.
*/
/**************************************************************************************/
void read_dump(   struct Field fldo,
				  double *t) {	
				  
	FILE *ht;
	int dump_version;
	int size_x,	size_y, size_z, included_field;
	int marker;

	DEBUG_START_FUNC;
	
	if( !file_exist(OUTPUT_DUMP) ) {
		// The file can't be openend
		ERROR_HANDLER(ERROR_CRITICAL, "Cannot open dump file.");
	}
	
	ht=fopen(OUTPUT_DUMP,"r");
	// The file has been opened by process 0
	
	if(rank==0) {
		fread(&dump_version, sizeof(int), 1, ht);
		if( dump_version != OUTPUT_DUMP_VERSION) ERROR_HANDLER( ERROR_CRITICAL, "Incorrect dump file version.");
		
		fread(&size_x		, sizeof(int), 1, ht);
		fread(&size_y		, sizeof(int), 1, ht);
		fread(&size_z		, sizeof(int), 1, ht);
		
		if(size_x != NX) ERROR_HANDLER( ERROR_CRITICAL, "Incorrect X grid size in dump file.");
		if(size_y != NY) ERROR_HANDLER( ERROR_CRITICAL, "Incorrect Y grid size in dump file.");
		if(size_z != NZ) ERROR_HANDLER( ERROR_CRITICAL, "Incorrect Y grid size in dump file.");
		
		fread(&included_field, sizeof(int), 1, ht);
	}
#ifdef MPI_SUPPORT
	MPI_Bcast( &included_field,		1, MPI_INT,		0, MPI_COMM_WORLD);
#endif
	
	MPI_Printf("Reading velocity field\n");
	read_field(ht, fldo.vx);
	read_field(ht, fldo.vy);
	read_field(ht, fldo.vz);
	
#ifdef BOUSSINESQ
	// Do we have Boussinesq Field in the dump?
	if(included_field & 1) {
		// Yes
		MPI_Printf("Reading Boussinesq field\n");
		read_field(ht, fldo.th);
	}
	else {
		// No
		ERROR_HANDLER( ERROR_WARNING, "No Boussinesq field in the dump, using initial conditions.");
	}
#endif
#ifdef MHD
	// Do we have MHD field in the dump?
	if(included_field & 2) {
		// Yes
		MPI_Printf("Reading MHD field\n");
		read_field(ht, fldo.bx);
		read_field(ht, fldo.by);
		read_field(ht, fldo.bz);
	}
	else {
		//No
		ERROR_HANDLER( ERROR_WARNING, "No MHD field in the dump, using initial conditions.");
	}
#endif
	
	if(param.restart) {		// If the dump is used to restart, we need these extra variables
		if(rank==0) {
			fread(t			, sizeof(double)		   , 1			   , ht);
	
			fread(&noutput_flow			, sizeof(int)			   , 1             , ht);
			fread(&lastoutput_time		, sizeof(double)		   , 1			   , ht);
			fread(&lastoutput_flow		, sizeof(double)		   , 1			   , ht);
			fread(&lastoutput_dump		, sizeof(double)		   , 1			   , ht);
	
			fread(&marker , sizeof(int)			   , 1, ht);
	
			if(marker != DUMP_MARKER) ERROR_HANDLER( ERROR_CRITICAL, "Incorrect marker. Probably an incorrect dump file!");	
			
		}
	
	// Transmit the values to all processes
#ifdef MPI_SUPPORT
		MPI_Bcast( t,					1, MPI_DOUBLE,	0, MPI_COMM_WORLD);
		MPI_Bcast( &noutput_flow,		1, MPI_INT,		0, MPI_COMM_WORLD);
		MPI_Bcast( &lastoutput_time,	1, MPI_DOUBLE,	0, MPI_COMM_WORLD);
		MPI_Bcast( &lastoutput_flow,	1, MPI_DOUBLE,	0, MPI_COMM_WORLD);
		MPI_Bcast( &lastoutput_dump,	1, MPI_DOUBLE,	0, MPI_COMM_WORLD);
#endif
	
		MPI_Printf("Restarting at t=%e...\n",*t);
	}

	fclose(ht);
	
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
	}

	DEBUG_END_FUNC;
	return;
}

/*********************************************************
*** General routine, callable from outside ***************
**********************************************************/
/**************************************************************************************/
/** 
	Initialize the output variables. Should be called only in the begining .
*/
/**************************************************************************************/

void init_output() {

	DEBUG_START_FUNC;
	
#ifdef WITH_SHEAR
// Initialize 1D arrays for remaps
	w1d = (double complex *) fftw_malloc( sizeof(double complex) * NY);
	if (w1d == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for w1d allocation");
	
	w2d = (double complex *) fftw_malloc( sizeof(double complex) * NY);
	if (w2d == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for w2d allocation");

// FFT plans (we use dummy arrays since we use the "guru" interface of fft3 in the code)
// The in place/ out of place will be set automatically at this stage

#ifdef _OPENMP
	fftw_plan_with_nthreads( 1 );
#endif

	fft_1d_forward = fftw_plan_dft_1d(NY, w1d, w2d, FFTW_FORWARD, FFT_PLANNING);
	if (fft_1d_forward == NULL) ERROR_HANDLER( ERROR_CRITICAL, "FFTW 1D forward plan creation failed");
	
	fft_1d_backward = fftw_plan_dft_1d(NY, w2d, w1d, FFTW_BACKWARD, FFT_PLANNING);
	if (fft_1d_backward == NULL) ERROR_HANDLER( ERROR_CRITICAL, "FFTW 1D backward plan creation failed");
	
#endif	
	// Check that the file restart exists
	
	if(param.restart) {
		if( !file_exist(OUTPUT_DUMP) ) {
			ERROR_HANDLER( ERROR_WARNING, "No restart dump found. I will set restart=false for this run.");
			param.restart = 0;
		}
	}
	
	if(!param.restart) {
		noutput_flow=0;
		lastoutput_time = param.t_initial - param.toutput_time;
		lastoutput_flow = param.t_initial - param.toutput_flow;
		lastoutput_dump = param.t_initial - param.toutput_dump;
	
		init1Dspectrum();
		clear_timevar();
	}
	
	DEBUG_END_FUNC;
	
	return;
}

/****************************************************************************/
/**
	Free the variables used by the output routines
*/
/****************************************************************************/
void finish_output() {
	DEBUG_START_FUNC;
	
#ifdef WITH_SHEAR
	fftw_free(w1d);
	fftw_free(w2d);
	
	fftw_destroy_plan(fft_1d_forward);
	fftw_destroy_plan(fft_1d_backward);	
#endif
	
	DEBUG_END_FUNC;
		
	return;
}
	
/**************************************************************************************/
/** 
	Check if an output (timevar, snapshot or dump) is required at t. If yes, call the 
	relevant routines.
	@param t Current time in the simulation
*/
/**************************************************************************************/

void output(const double t) {
	
	DEBUG_START_FUNC;
	
	// Very rough output function
	if( (t-lastoutput_time)>=param.toutput_time) {
		output_timevar(fld,t);
		output1Dspectrum(fld,t);
		lastoutput_time = lastoutput_time + param.toutput_time;
	}
	
	if( (t-lastoutput_flow)>=param.toutput_flow) {
		output_vtk(noutput_flow,t);
		noutput_flow++;
		lastoutput_flow = lastoutput_flow + param.toutput_flow;
	}
		
	if( (t-lastoutput_dump)>=param.toutput_dump) {
		lastoutput_dump=lastoutput_dump+param.toutput_dump;
		output_dump(fld,t);
	}

	DEBUG_END_FUNC;
	
	return;
}

/**************************************************************************************/
/** 
	Show the current status of the output routine.
	@param iostream Handler of the file in which the status is written.
*/
/**************************************************************************************/
void output_status(FILE * iostream) {
	if(rank==0)
		fprintf(iostream,"Next output in file n %d, at t=%e\n",noutput_flow, lastoutput_flow+param.toutput_flow);
	return;
}

/**************************************************************************************/
/** 
	Immediatly output a timevar, snapshot and dump, regardless of the output
	parameters.
	@param t Current time of the simulation.
*/
/**************************************************************************************/

void output_immediate(const double t) {
	// Very rough output function
	// Immediate output
	output_timevar(fld,t);
	output_vtk(noutput_flow,t);
	noutput_flow++;
	output_dump(fld,t);
	return;
}

/**************************************************************************************/
/** 
	Immediatly output a dump file, regardless of the output
	parameters.
	@param t Current time of the simulation.
*/
/**************************************************************************************/

void dump_immediate(const double t) {
	output_dump(fld,t);
	return;
}
