// routines for particle dynamic

#include <stdlib.h>

#include "common.h"
#include "gfft.h"
#include "debug.h"
#include "shear.h"
#include "transpose.h"

#ifdef WITH_PARTICLES
// init particle positions
double *x3D;	/**< x positions in space */
double *y3D;	/**< y positions in space */
double *z3D;	/**< z positions in space */

#ifdef WITH_SHEAR
fftw_plan	fft_particle_forward;								/**< 1D FFT transforms. Used by remap routines.*/
fftw_plan	fft_particle_backward;							/**< 1D FFT transforms. Used by remap routines.*/

/***********************************************************/
/** 
	Remap a real field from the current sheared frame to the classical 
	cartesian frame. This remap routine assumes the incoming field 
	has dimensions NX/NPROC+1, NY+1, NZ+1
	
	@param wri Real array to be remapped. The remapped array will be stored here.
	@param t Current time of the simulation (used to compute what the sheared frame is).
	
*/
/***********************************************************/

void remap_flow(	double wri[], 
					const double t,
					const double tremap) {
					
	int i,j,k;
	double tvelocity;
	complex double wexp;
	complex double phase;
	double complex		*wremap;						/** 1D arrays used by the remap methods */
	
	DEBUG_START_FUNC;
	
	tvelocity = fmod(t, 2.0 * param.ly / (param.shear * param.lx));
	
	wremap = (double complex *) fftw_malloc( sizeof(double complex) * (NY/2+1) * (NZ+1) );
	if (wremap == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for wremap allocation");
	
	for( i = 0 ; i < NX/NPROC+1 ; i++) {
		fftw_execute_dft_r2c(fft_particle_forward, wri + i*(NZ+1)*(NY+1), wremap);
		for( j = 0 ; j < NY/2+1 ; j++) {
			phase = (double complex) ((2.0 * M_PI) / param.ly *  ((double) j )  * 
									( ((double) (i + rank * (NX/NPROC)) / (double) NX ) * tremap - tvelocity / 2.0 ) * param.lx * param.shear);
			
			wexp = cexp( I * phase) / NY;
			
			for( k = 0 ; k < NZ+1; k++) {
				wremap[ k + j * (NZ+1) ] = wexp * wremap[ k + j * (NZ+1) ];
			}
		}
		fftw_execute_dft_c2r(fft_particle_backward, wremap, wri+i*(NZ+1)*(NY+1));
	}

	fftw_free(wremap);
	
	DEBUG_END_FUNC;
	
	return;
}

/***********************************************************/
/** 
	Remap a real field from the galilean frame to the classical 
	sheared frame. This remap routine assumes the incoming field 
	has dimensions NX/NPROC+1, NY+1, NZ+1
	
	@param wri Real array to be remapped. The remapped array will be stored here.
	@param t Current time of the simulation (used to compute what the sheared frame is).
	
*/
/***********************************************************/

void remap_flow_inv( double wri[], 
					 const double t,
					 const double tremap) {
					
	int i,j,k;
	double tvelocity;
	complex double wexp;
	complex double phase;
	double complex		*wremap;						/** 1D arrays used by the remap methods */
	
	DEBUG_START_FUNC;
	
	tvelocity = fmod(t, 2.0 * param.ly / (param.shear * param.lx));
	
	wremap = (double complex *) fftw_malloc( sizeof(double complex) * (NY/2+1) * (NZ+1) );
	if (wremap == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for wremap allocation");
	
	for( i = 0 ; i < NX/NPROC+1 ; i++) {
		fftw_execute_dft_r2c(fft_particle_forward, wri + i*(NZ+1)*(NY+1), wremap);
		for( j = 0 ; j < NY/2+1 ; j++) {
			phase = - (double complex) ((2.0 * M_PI) / param.ly *  ((double) j )  * 
									( ((double) (i + rank * (NX/NPROC)) / (double) NX ) * tremap - tvelocity / 2.0 ) * param.lx * param.shear);
			
			wexp = cexp( I * phase) / NY;
			
			for( k = 0 ; k < NZ+1; k++) {
				wremap[ k + j * (NZ+1) ] = wexp * wremap[ k + j * (NZ+1) ];
			}
		}
		fftw_execute_dft_c2r(fft_particle_backward, wremap, wri+i*(NZ+1)*(NY+1));
	}

	fftw_free(wremap);
	
	DEBUG_END_FUNC;
	
	return;
}



#endif

void output_partvar(struct Particle *part, double t) {
	int i;
	double vxmax,vxmin,vymax,vymin,vzmax,vzmin;
	double vxm,vym,vzm,vx2,vy2,vz2;
	FILE *ht;
	
	DEBUG_START_FUNC;
	
	vxm=0.0;
	vym=0.0;
	vzm=0.0;
	vx2=0.0;
	vy2=0.0;
	vz2=0.0;
	
	vxmax=part[0].vx;
	vymax=part[0].vy;
	vzmax=part[0].vz;
	
	vxmin=part[0].vx;
	vymin=part[0].vy;
	vzmin=part[0].vz;
	
	// Compute average velocity and max/min
	for(i = 0 ; i < param.particles_n/NPROC ; i++) {
		vxm += part[i].vx;
		vym += part[i].vy;
		vzm += part[i].vz;
		
		if(vxmax<part[i].vx) vxmax = part[i].vx;
		if(vymax<part[i].vy) vymax = part[i].vy;
		if(vzmax<part[i].vz) vzmax = part[i].vz;
		if(vxmin>part[i].vx) vxmin = part[i].vx;
		if(vymin>part[i].vy) vymin = part[i].vy;
		if(vzmin>part[i].vz) vzmin = part[i].vz;
	}
	
	reduce(&vxm, 1);
	reduce(&vym, 1);
	reduce(&vzm, 1);
	
	vxm=vxm/param.particles_n;
	vym=vym/param.particles_n;
	vzm=vzm/param.particles_n;
	
	reduce(&vxmax, 2);
	reduce(&vymax, 2);
	reduce(&vzmax, 2);
	
	reduce(&vxmin, 3);
	reduce(&vymin, 3);
	reduce(&vzmin, 3);
	
	// Compute mean deviation
	for(i = 0 ; i < param.particles_n/NPROC ; i++) {
		vx2 += (part[i].vx-vxm)*(part[i].vx-vxm);
		vy2 += (part[i].vy-vym)*(part[i].vy-vym);
		vz2 += (part[i].vz-vzm)*(part[i].vz-vzm);
	}
	
	reduce(&vx2,1);
	reduce(&vy2,1);
	reduce(&vz2,1);
	
	vx2=vx2/param.particles_n;
	vy2=vy2/param.particles_n;
	vz2=vz2/param.particles_n;
	
	// Root mean square
	vx2=pow(vx2,0.5);
	vy2=pow(vy2,0.5);
	vz2=pow(vz2,0.5);
	
	if(rank==0) {
		ht=fopen("partvar","a");
		fprintf(ht,"%08e\t",t);
		fprintf(ht,"%08e\t%08e\t%08e\t%08e\t%08e\t%08e\t",vxmax,vxmin,vymax,vymin,vzmax,vzmin);
		fprintf(ht,"%08e\t%08e\t%08e\t",vxm,vym,vzm);
		fprintf(ht,"%08e\t%08e\t%08e",vx2,vy2,vz2);
		fprintf(ht,"\n");
		
		if(ferror(ht)) ERROR_HANDLER( ERROR_CRITICAL, "Error writing partvar file");
		
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
	
void write_particles_mass(FILE *ht, struct Particle *part) {
	int i;
	float q0;
	
	for(i = 0 ; i < param.particles_n/NPROC ; i++) {
		q0 = big_endian( (float) part[i].mass);
		fwrite( &q0, sizeof(float), 1, ht);
	}
	return;
}


void write_particles_position(FILE *ht, struct Particle *part) {
	int i;
	float q0;
	
	for(i = 0 ; i < param.particles_n/NPROC ; i++) {
		q0 = big_endian( (float) part[i].x);
		fwrite( &q0, sizeof(float), 1, ht);
		
		q0 = big_endian( (float) part[i].y);
		fwrite( &q0, sizeof(float), 1, ht);
		
		q0 = big_endian( (float) part[i].z);
		fwrite( &q0, sizeof(float), 1, ht);
	}
	return;
}
	
/***********************************************************/
/** 
	Output the particle positions and mass to a VTK file
	using a POLYDATA format. The particle position can then be 
	easily represented in paraview using glyphs.
	
	@param n file number to be created
	@param t Current time of the simulation.
	
*/
/***********************************************************/

void output_particles(struct Field fldi, const int n, double t) {

	FILE *ht = NULL;
	char  filename[50];
	int i;
#ifdef MPI_SUPPORT
	struct Particle *part_chunk;
	MPI_Status status;
#endif
	
	DEBUG_START_FUNC;

#ifdef MPI_SUPPORT
	if(rank==0)
		part_chunk = (struct Particle *) malloc( param.particles_n/NPROC * sizeof(struct Particle) );
#endif
	
	sprintf(filename,"data/p%04i.vtk",n);

	if(rank==0) {
		ht=fopen(filename,"w");
	
		fprintf(ht, "# vtk DataFile Version 2.0\n");
		fprintf(ht, "t= %015.15e Snoopy Code v5.0\n",t);
		fprintf(ht, "BINARY\n");
		fprintf(ht, "DATASET POLYDATA\n");
		fprintf(ht, "POINTS %d float\n",param.particles_n);
	}
	// Write Particle position...
	
#ifdef MPI_SUPPORT
	if(rank == 0) {
		for(i = 0 ; i < NPROC ; i++) {
			if(i==0) {
				// Write down the local process
				write_particles_position(ht, fldi.part);
			}
			else {
				// Then gather the other processes
				MPI_Recv( part_chunk, sizeof(struct Particle) * param.particles_n/NPROC, MPI_BYTE, i, 2, MPI_COMM_WORLD, &status);
				write_particles_position(ht, part_chunk);
			}
		}
	}
	else {
		MPI_Send( fldi.part, sizeof(struct Particle) * param.particles_n/NPROC, MPI_BYTE, 0, 2, MPI_COMM_WORLD);
	}
	
#else
		
	write_particles_position(ht, fldi.part);

#endif
	
	if(rank==0) {
		fprintf(ht, "POINT_DATA %d\n",param.particles_n);
		fprintf(ht, "SCALARS %s float\n","mass");
		fprintf(ht, "LOOKUP_TABLE default\n");
	}
		
#ifdef MPI_SUPPORT
	if(rank == 0) {
		for(i = 0 ; i < NPROC ; i++) {
			if(i==0) {
				// Write down the local process
				write_particles_mass(ht, fldi.part);
			}
			else {
				// Then gather the other processes
				MPI_Recv( part_chunk, sizeof(struct Particle) * param.particles_n/NPROC, MPI_BYTE, i, 2, MPI_COMM_WORLD, &status);
				write_particles_mass(ht, part_chunk);
			}
		}
	}
	else {
		MPI_Send( fldi.part, sizeof(struct Particle) * param.particles_n/NPROC, MPI_BYTE, 0, 2, MPI_COMM_WORLD);
	}
	
#else

	write_particles_mass(ht,fldi.part);
	
#endif
	
	if(rank==0) {
		if(ferror(ht)) ERROR_HANDLER( ERROR_CRITICAL, "Error writing particles-vtk file");
		fclose(ht);
	}
	
#ifdef MPI_SUPPORT
	if(rank==0) free(part_chunk);
#endif

	DEBUG_END_FUNC; 
	
	return;
}

/***********************************************************/
/** 
	Output a density field for the particle positions
	in a VTK file. This routine assumes the VTK header has already
	be written in file *ht. 
	
	
	@param ht handler to the VTK file
	@param t Current time of the simulation.
	
*/
/***********************************************************/

void write_vtk_particles(struct Field fldi, FILE * ht, const double t) {
	int i,j,k;
	float q0;
	
	
	double *density;
	double *MassArray;
	int *indexArray;

	DEBUG_START_FUNC;
	
	density = (double *) fftw_malloc( (NX/NPROC+1) * (NY+1) * (NZ+1) * sizeof(double ));
	if (density == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for density allocation");

	indexArray = (int *) malloc( param.particles_n/NPROC * 3 * sizeof(int) );	//3 spatial coordinates per particle
	MassArray = (double *) malloc( param.particles_n/NPROC * sizeof(double) ); // interpolated vx values at particle location
	
	for(i=0 ; i < param.particles_n/NPROC ; i++) {
		MassArray[i] = fldi.part[i].mass;
	}
	
	build_indexArray(fldi, t, indexArray);
	
	for( i = 0 ; i < (NX/NPROC+1) * (NY+1) * (NZ+1) ; i++) {
		density[ i ] = 0.0;
	}
	
	particles_to_flow(indexArray, fldi, MassArray, density);
	
	
	for( k = 0 ; k < NZ; k++) {
		for( j = 0; j < NY; j++) {	
			for( i = 0; i < NX; i++) {
				q0 = big_endian( (float) (density[k + j * (NZ + 1) + i * (NY + 1) * (NZ + 1) ] * NX * NY * NZ / (param.lx * param.ly * param.lz)));
				fwrite(&q0, sizeof(float), 1, ht);
			}
		}
	}
	
	free(density);
	free(indexArray);
	free(MassArray);
	
	DEBUG_END_FUNC;
	
	return;
}

void read_particle_dump(FILE *ht, struct Particle *part) {

	int q0;
	
	#ifdef MPI_SUPPORT
	MPI_Status status;
	int current_rank;
	struct Particle * part_chunk;

	if(rank==0) {
		part_chunk = (struct Particle *) malloc( param.particles_n/NPROC * sizeof(struct Particle) );

		for(current_rank=0; current_rank < NPROC; current_rank++) {
			if(current_rank==0) {
				// read the rank=0 process dump
				fread(&q0, sizeof(int), 1, ht);
				if(q0 != param.particles_n) ERROR_HANDLER( ERROR_CRITICAL, "A different number of particles was saved in the dump.");
				fread(part, sizeof(struct Particle), param.particles_n/NPROC, ht);
			}
			
			else {
				fread(part_chunk, sizeof(struct Particle), param.particles_n/NPROC, ht);
				MPI_Send( part_chunk, param.particles_n/NPROC*sizeof(struct Particle), MPI_BYTE, current_rank, 2, MPI_COMM_WORLD);
			}

		}
		MPI_Barrier(MPI_COMM_WORLD);
		
		free(part_chunk);
	}
	else {
		MPI_Recv(part, param.particles_n/NPROC*sizeof(struct Particle), MPI_BYTE, 0, 2, MPI_COMM_WORLD, &status);
		MPI_Barrier(MPI_COMM_WORLD);
	}

#else
	fread(&q0, sizeof(int), 1, ht);
	if(q0 != param.particles_n) ERROR_HANDLER( ERROR_CRITICAL, "A different number of particles was saved in the dump.");
	fread(part, sizeof(struct Particle), param.particles_n/NPROC, ht);
#endif
}

void write_particle_dump(FILE *ht, struct Particle *part) {

	int q0;
	
#ifdef MPI_SUPPORT
	MPI_Status status;
	int current_rank;
	struct Particle * part_chunk;

	if(rank==0) {
		part_chunk = (struct Particle *) malloc( param.particles_n/NPROC * sizeof(struct Particle) );

		for(current_rank=0; current_rank < NPROC; current_rank++) {
			if(current_rank==0) {
				// write the rank=0 process dump
				q0=param.particles_n;
				fwrite(&q0, sizeof(int), 1, ht);
				
				fwrite(part, sizeof(struct Particle), param.particles_n/NPROC, ht);
			}
			
			else {
				MPI_Recv( part_chunk, param.particles_n/NPROC*sizeof(struct Particle), MPI_BYTE, current_rank, 2, MPI_COMM_WORLD, &status);
				fwrite(part_chunk, sizeof(struct Particle), param.particles_n/NPROC, ht);
			}

		}
		MPI_Barrier(MPI_COMM_WORLD);
		
		free(part_chunk);
	}
	else {
		MPI_Send(part, param.particles_n/NPROC*sizeof(struct Particle), MPI_BYTE, 0, 2, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
	}

#else
	q0=param.particles_n;
	fwrite(&q0, sizeof(int), 1, ht);
	fwrite(part, sizeof(struct Particle), param.particles_n/NPROC, ht);
#endif


}
/***********************************************************/
/** 
	Init the particle module. This routine has 3 main purposes:
	- It initializes the "grid" which is used later on for interpolations
	- It initializes the particle positions and velocities
	- It initializes the ffts used in the remap routines when shear
	   is present.
	
*/
/***********************************************************/

void init_particles(struct Field fldi) {
	int i,j,k;
#ifdef WITH_SHEAR
	double complex		*wremap;
	double *vx;
#endif
	double vx0;
	double vy0;
	double kappa_tau2;
	
	const int n_size1D[1] = {NY};
	
	DEBUG_START_FUNC;
	
	// In case NPARTICLE % NPROC /= 0
#ifdef MPI_SUPPORT
	if( param.particles_n % NPROC ) ERROR_HANDLER( ERROR_CRITICAL, "param.particles_n have to be a multiple of NPROC");
#endif
#ifdef WITH_2D
	ERROR_HANDLER( ERROR_CRITICAL, "Particles module incompatible with the 2D optimization");
#endif

	// init a real space grid
	x3D = (double *) malloc( (NX/NPROC+1) * (NY+1) * (NZ+1) * sizeof(double ));
	y3D = (double *) malloc( (NX/NPROC+1) * (NY+1) * (NZ+1) * sizeof(double ));
	z3D = (double *) malloc( (NX/NPROC+1) * (NY+1) * (NZ+1) * sizeof(double ));
	
	for(i = 0 ; i < NX/NPROC+1 ; i++) {
		for( j = 0 ; j < NY+1 ; j++) {
			for( k = 0 ; k < NZ+1 ; k++) {
				x3D[ k + (NZ+1) * j + (NZ+1) * (NY+1) * i ] = -param.lx / 2.0 + ((i + rank * (NX/NPROC)) * param.lx / ((double) NX));
				y3D[ k + (NZ+1) * j + (NZ+1) * (NY+1) * i ] = -param.ly / 2.0 + j * param.ly / ((double) NY);
				z3D[ k + (NZ+1) * j + (NZ+1) * (NY+1) * i ] = -param.lz / 2.0 + k * param.lz / ((double) NZ);
			}
		}
	}
	
	// init particle positions and velocity
	
	/*
	for(i = 0 ; i < 100/NPROC ; i++) {
		for(j = 0 ; j < 100 ; j++) {
			fldi.part[ j + 100 * i ].x = -param.lx / 2.0 + param.lx * (i+rank*100/NPROC) / 100;
			fldi.part[ j + 100 * i ].y = -param.ly / 2.0 + param.ly * j / 100;
			fldi.part[ j + 100 * i ].z = 0.0;
		}
	}
	*/
	/*
	fldi.part[0].x=0.0;
	fldi.part[0].y=0.0;
	fldi.part[0].z=0.0;
	*/
	
	kappa_tau2 = 2.0*param.omega*(2.0*param.omega-param.shear) * param.particles_stime * param.particles_stime + (param.particles_dg_ratio + 1.0) * (param.particles_dg_ratio + 1.0);

	// This is a non trivial equilibrium for the particles+gas system
	
	vx0=-param.particles_epsilon*param.particles_stime/kappa_tau2;
	vy0= param.particles_epsilon*(kappa_tau2-(param.particles_dg_ratio+1.0)) / ( 2.0*param.omega*kappa_tau2 );
	
	for(i = 0 ; i < param.particles_n/NPROC ; i++) {
		//fldi.part[i].x = 0.4-0.1*i;
		//fldi.part[i].y = 0.0;
		//fldi.part[i].z = 0.0;
		fldi.part[i].x = -param.lx/2.0 +param.lx*randm();
		fldi.part[i].y = -param.ly/2.0 +param.ly*randm();
		fldi.part[i].z = -param.lz/2.0 +param.lz*randm();
		
		fldi.part[i].vx = 0.5-randm();//vx0;
		fldi.part[i].vy = 0.5-randm();//vy0;
		fldi.part[i].vz = 0.5-randm();
		//fldi.part[i].mass = param.particles_mass;
		fldi.part[i].mass = param.particles_dg_ratio * param.lx * param.ly * param.lz / param.particles_n;
		fldi.part[i].stime = param.particles_stime;
	}

	// Init velocity vectors
	
	
	
#ifdef WITH_SHEAR
// Initialize 1D arrays for remaps
	wremap = (double complex *) fftw_malloc( sizeof(double complex) * (NY/2+1) * (NZ+1) );
	if (wremap == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for wremap allocation");

	vx = (double *) fftw_malloc( (NX/NPROC+1) * (NY+1) * (NZ+1) * sizeof(double ));
	if (vx == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for vx allocation");


// FFT plans (we use dummy arrays since we use the "guru" interface of fft3 in the code)
// The in place/ out of place will be set automatically at this stage

// The Following Fourier transforms takes an array of size ( NX+1, NY+1, NZ+1) but consider only the "included" array
// of size ( NX+1, NY, NZ+1) and transforms it in y, in an array of size (NX+1, NY/2+1, NZ+1). The j=NY plane is therefore 
// not modified by these calls

#ifdef _OPENMP
	fftw_plan_with_nthreads( 1 );
#endif

	fft_particle_forward = fftw_plan_many_dft_r2c(1, n_size1D, NZ+1,
											vx, NULL, NZ+1, 1,
											wremap, NULL, NZ+1, 1,
											FFT_PLANNING | FFTW_UNALIGNED );
											
	if (fft_particle_forward == NULL) ERROR_HANDLER( ERROR_CRITICAL, "FFTW 1D forward plan creation failed");
	
	fft_particle_backward = fftw_plan_many_dft_c2r(1, n_size1D, NZ+1,
											wremap, NULL, NZ+1, 1,
											vx,  NULL, NZ+1, 1,
											FFT_PLANNING | FFTW_UNALIGNED);
											
	if (fft_particle_backward == NULL) ERROR_HANDLER( ERROR_CRITICAL, "FFTW 1D backward plan creation failed");
	
	fftw_free(wremap);
	fftw_free(vx);
#endif	


	// All done
	
	DEBUG_END_FUNC;
	
	return;
}

void map_flow_forces(double *qi,
					 double *qo,
					 const double t,
					 const double tremap) {
					 
	
	int i,j,k;
	
	DEBUG_START_FUNC;
	
	// Sum the ghost zones in each direction
	
	// Periodize in the z direction k=0->k=NZ
	for( i = 0  ; i < NX/NPROC+1 ; i++) {
		for( j = 0 ; j < NY+1 ; j++) {
			qi[j * (NZ+1) + i * (NZ+1) * (NY+1)] += qi[NZ + j * (NZ+1) + i * (NZ+1) * (NY+1)];
		}
	}
	
	// Periodize in the y direction j=0->j=NY
	
	for( i = 0  ; i < NX/NPROC+1 ; i++) {
		for( k = 0 ; k < NZ ; k++) {
			qi[k + i * (NZ+1) * (NY+1)] += qi[k + NY * (NZ+1) + i * (NZ+1) * (NY+1)];
		}
	}

#ifdef WITH_SHEAR
	// Remap the flow
	remap_flow_inv(qi, t, tremap);
#endif

	// Periodize the output array in the x direction i=O->i=NX
	
	for( j = 0  ; j < NY ; j++) {
		for( k = 0 ; k < NZ ; k++) {
			qi[k + j * (NZ + 1)] += qi[k + j * (NZ+1) + NX * (NZ+1) * (NY+1)];
		}
	}

#ifdef _OPENMP
#pragma omp parallel for private(i,j,k) schedule(static)	
#endif
	for( i = 0 ; i < NX/NPROC ; i++) {
		for( j = 0  ; j < NY ; j++) {
			for( k = 0 ; k < NZ ; k++) {
				qo[k + j * (NZ + 2) + i * (NZ+2) * NY] = qi[k + j * (NZ+1) + i * (NZ+1) * (NY+1)];
			}
		}
	}
	
	DEBUG_END_FUNC;
	
	return;
}



/***********************************************************/
/** 
	Computes the flow velocity (minus shear if shear is present)
	in the traditional rectangular frame. It takes as an input the velocity
	field in a sheared frame (such as the one given by a call to gfft_c2r(fldi.qi))
	and produces an unsheared velocity field to which ghost zones are added in each direction
	
	The ghost zones are used by the interpolation routines later on
	
	We assume the input arrays have size (NX, NY, NZ+2) without mpi or (NY/NPROC, NX, NZ+2) with mpi
	
	CAUTION: With MPI, this routine uses temporary arrays wr4, wr5 and wr6
	
	@param qxi: input x velocity component
	@param qyi: input y velocity component
	@param qzi: input z velocity component
	@param qxo: output x velocity component
	@param qyo: output y velocity component
	@param qzo: output z velocity component
	@param t: current time (used by the remap routines)
	
*/
/***********************************************************/


void compute_flow_velocity(double *qxi,
						   double *qyi,
						   double *qzi,
						   double *qxo,
						   double *qyo,
						   double *qzo,
						   const double t,
						   const double tremap) {
	
	int i,j,k;
	int dest, source;
	
#ifdef MPI_SUPPORT
	MPI_Status status;
#endif
	DEBUG_START_FUNC;
	
#ifdef MPI_SUPPORT
	// Do a transposition of the array and store it in work arrays
	
	transpose_real(NY, NX, NZ+2, NPROC, qxi, wr4);
	transpose_real(NY, NX, NZ+2, NPROC, qyi, wr5);
	transpose_real(NY, NX, NZ+2, NPROC, qzi, wr6);
	
#ifdef _OPENMP
	#pragma omp parallel for private(i,j,k) schedule(static)	
#endif
	for( i = 0 ; i < NX/NPROC ; i++) {
		for( j = 0  ; j < NY ; j++) {
			for( k = 0 ; k < NZ ; k++) {
				qxo[k + j * (NZ+1) + i * (NZ+1) * (NY+1)] = wr4[k + j * (NZ + 2) + i * (NZ+2) * NY] / ((double)NTOTAL);
				qyo[k + j * (NZ+1) + i * (NZ+1) * (NY+1)] = wr5[k + j * (NZ + 2) + i * (NZ+2) * NY] / ((double)NTOTAL);
				qzo[k + j * (NZ+1) + i * (NZ+1) * (NY+1)] = wr6[k + j * (NZ + 2) + i * (NZ+2) * NY] / ((double)NTOTAL);
			}
		}
	}

#else	// NO MPI
	
	// Reshape the input array into a (NX+1)(NY+1)(NZ+1) array
#ifdef _OPENMP
	#pragma omp parallel for private(i,j,k) schedule(static)	
#endif
	for( i = 0 ; i < NX/NPROC ; i++) {
		for( j = 0  ; j < NY ; j++) {
			for( k = 0 ; k < NZ ; k++) {
				qxo[k + j * (NZ+1) + i * (NZ+1) * (NY+1)] = qxi[k + j * (NZ + 2) + i * (NZ+2) * NY] / ((double)NTOTAL);
				qyo[k + j * (NZ+1) + i * (NZ+1) * (NY+1)] = qyi[k + j * (NZ + 2) + i * (NZ+2) * NY] / ((double)NTOTAL);
				qzo[k + j * (NZ+1) + i * (NZ+1) * (NY+1)] = qzi[k + j * (NZ + 2) + i * (NZ+2) * NY] / ((double)NTOTAL);
			}
		}
	}
	
#endif

#ifdef MPI_SUPPORT
	// Transmit the data boundaries accross the processors to make interpolation easier...
	source = (rank+1) % NPROC;
	dest = rank-1;
	
	if(dest < 0) dest = NPROC-1;
	
	MPI_Sendrecv( qxo						  , (NY+1)*(NZ+1), MPI_DOUBLE, dest, 1, 
				  qxo + NX*(NY+1)*(NZ+1)/NPROC, (NY+1)*(NZ+1), MPI_DOUBLE, source, 1,
				  MPI_COMM_WORLD, &status );
				  
	MPI_Sendrecv( qyo						  , (NY+1)*(NZ+1), MPI_DOUBLE, dest, 1, 
				  qyo + NX*(NY+1)*(NZ+1)/NPROC, (NY+1)*(NZ+1), MPI_DOUBLE, source, 1,
				  MPI_COMM_WORLD, &status );
				  
	MPI_Sendrecv( qzo						  , (NY+1)*(NZ+1), MPI_DOUBLE, dest, 1, 
				  qzo + NX*(NY+1)*(NZ+1)/NPROC, (NY+1)*(NZ+1), MPI_DOUBLE, source, 1,
				  MPI_COMM_WORLD, &status );
	
#else
	// Periodize the output array in the x direction i=O->i=NX
	
	for( j = 0  ; j < NY ; j++) {
		for( k = 0 ; k < NZ ; k++) {
			qxo[k + j * (NZ+1) + NX * (NZ+1) * (NY+1)] = qxo[k + j * (NZ + 1)];
			qyo[k + j * (NZ+1) + NX * (NZ+1) * (NY+1)] = qyo[k + j * (NZ + 1)];
			qzo[k + j * (NZ+1) + NX * (NZ+1) * (NY+1)] = qzo[k + j * (NZ + 1)];
		}
	}
#endif
	
#ifdef WITH_SHEAR
	// Remap the flow
	
	remap_flow(qxo, t, tremap);
	remap_flow(qyo, t, tremap);
	remap_flow(qzo, t, tremap);
	
#endif

	// Periodize in the y direction j=0->j=NY
	
	for( i = 0  ; i < NX/NPROC+1 ; i++) {
		for( k = 0 ; k < NZ ; k++) {
			qxo[k + NY * (NZ+1) + i * (NZ+1) * (NY+1)] = qxo[k + i * (NZ+1) * (NY+1)];
			qyo[k + NY * (NZ+1) + i * (NZ+1) * (NY+1)] = qyo[k + i * (NZ+1) * (NY+1)];
			qzo[k + NY * (NZ+1) + i * (NZ+1) * (NY+1)] = qzo[k + i * (NZ+1) * (NY+1)];

		}
	}
	
	// Periodize in the z direction k=0->k=NZ
	for( i = 0  ; i < NX/NPROC+1 ; i++) {
		for( j = 0 ; j < NY+1 ; j++) {
			qxo[NZ + j * (NZ+1) + i * (NZ+1) * (NY+1)] = qxo[j * (NZ+1) + i * (NZ+1) * (NY+1)];
			qyo[NZ + j * (NZ+1) + i * (NZ+1) * (NY+1)] = qyo[j * (NZ+1) + i * (NZ+1) * (NY+1)];
			qzo[NZ + j * (NZ+1) + i * (NZ+1) * (NY+1)] = qzo[j * (NZ+1) + i * (NZ+1) * (NY+1)];
		}
	}
	
	// All done...
	
	DEBUG_END_FUNC;
	return;
}

#ifdef MPI_SUPPORT
int alltoall_schedule(int npes, int myrank, int step)
{
  int pe, n;

  if ((npes & 1) == 1)		/* skip self communication step for odd npes */
    step++;

  if (step == 0)		/* communicate with self */
    return myrank;

  if (myrank == npes - 1)	/* last processor communicates ordered */
    return step - 1;
  if (step == myrank + 1)
    return npes - 1;

  if ((npes & 1) == 1) {	/* re-insert self communication step at stall */
    if (myrank & 1) {
      if (step == (myrank + 1)/2)
	return myrank;
    } else {
      if (step == (npes + myrank + 1)/2)
	return myrank;
    }
  }

  n = npes + (npes & 1);	/* the next even number, starting from npes */
  pe = 2 * (step - 1) - myrank;	/* mirror along the (pe,npes-1) point */
  if (pe < 0) {
    pe += n - 1;
  } else if (pe >= n - 1) {
    pe -= n - 1;
  }
  return pe;
}
#endif


void build_indexArray(struct  Field fldi,
					  double  t,
						 int  *indexArray) {
						
	// Build an indexArray so that indexArray[3*i+j] is the the grid coordinate of particle i (x->j=0, y->j=1, z->j=2)
	// This routine assumes indexArray was already allocated
	int i,m,n,p;
	double x,y,z;
	
#ifdef _OPENMP
#pragma omp parallel for private(i, x, y, z, m, n, p) schedule(static)	
#endif
	for( i = 0 ; i < param.particles_n/NPROC ; i++) {
		// Compute particle position in the box
		x = fldi.part[i].x - param.lx * floor( fldi.part[i].x / param.lx + 0.5 );
#ifdef WITH_SHEAR
		y = fldi.part[i].y + (fldi.part[i].x - x) * param.shear * t;
#else
		y = fldi.part[i].y;
#endif
		y = y - param.ly * floor( y / param.ly + 0.5 );
		z = fldi.part[i].z - param.lz * floor( fldi.part[i].z / param.lz + 0.5 );

		// particle indices (these indices could be outside of the local box!!)
		m=(int) floor( (x/param.lx+0.5) * NX);
		n=(int) floor( (y/param.ly+0.5) * NY);
		p=(int) floor( (z/param.lz+0.5) * NZ);

		indexArray[3*i]   = m;
		indexArray[3*i+1] = n;
		indexArray[3*i+2] = p;
	}
	
	return;
}


	
	
// Map grid based scalar to a particle based scalar
void flow_to_particles(int *indexArray, struct Field fldi, double *vin, double *Varray) {
	int m,n,p;
	int i,j;
	
	double dx, dy, dz;
	double q000, q001, q010, q011, q100, q101, q110, q111;
	
#ifdef MPI_SUPPORT
	int nstep;
	int target;
	int nsend;
	int nrecv;
	
	MPI_Status	status;
		
	static int *indexSend = NULL;
	static int *indexRecv = NULL;
	static int *indexPos = NULL;
	static double *Vsend = NULL;
	static double *Vrecv = NULL;
	
	DEBUG_START_FUNC;

	// Allocate memory if it is the first time this function is called
	
	if(indexSend == NULL) {
		indexSend = (int *) malloc( 3 * sizeof(int) * param.particles_n/NPROC );
		if(indexSend == NULL) ERROR_HANDLER(ERROR_CRITICAL, "No memory for indexSend allocation");
		
		indexRecv = (int *) malloc( 3 * sizeof(int) * param.particles_n/NPROC );
		if(indexRecv == NULL) ERROR_HANDLER(ERROR_CRITICAL, "No memory for indexRecv allocation");
		
		indexPos = (int *) malloc(      sizeof(int) * param.particles_n/NPROC );
		if(indexPos == NULL) ERROR_HANDLER(ERROR_CRITICAL, "No memory for indexPos allocation");
		
		Vsend = (double *) malloc( 8 * sizeof(double) * param.particles_n/NPROC );
		if(Vsend == NULL) ERROR_HANDLER(ERROR_CRITICAL, "No memory for Vsend allocation");
		
		Vrecv = (double *) malloc( 8 * sizeof(double) * param.particles_n/NPROC );
		if(Vrecv == NULL) ERROR_HANDLER(ERROR_CRITICAL, "No memory for Vrecv allocation");
	}
	
	// Loop on each process
	for( nstep = 0; nstep < NPROC ; nstep++) {
		
		target = alltoall_schedule( NPROC, rank, nstep);	// Find the target processor
		
		nsend = 0;
		// Find the information we need from target processor
		for( i = 0 ; i < param.particles_n/NPROC ; i++) {
			if( (indexArray[ 3*i ] / (NX/NPROC) ) == target ) {
				//particle i belongs to the flow in processor target
				
				indexPos[ nsend ] = i; 
				indexSend[ 3*nsend ] = indexArray[ 3*i ];
				indexSend[ 3*nsend + 1] = indexArray[ 3*i + 1 ];
				indexSend[ 3*nsend + 2] = indexArray[ 3*i + 2 ];
		
				nsend++;
			}
		}
		
		// We're going to need nsend flow velocities points from processor target
			
		// Exchange the array sizes
			
		MPI_Sendrecv( &nsend, 1, MPI_INT, target, 1, 
					  &nrecv, 1, MPI_INT, target, 1,
					  MPI_COMM_WORLD, &status );
		
		// Exchange the indices
		MPI_Sendrecv( indexSend, 3*nsend, MPI_INT, target, 1,
					  indexRecv, 3*nrecv, MPI_INT, target, 1,
					  MPI_COMM_WORLD, &status );
		
		// From that point indexRecv contains the indices of the velocity field needed by the target.
			
		for( i = 0 ; i < nrecv ; i++) {
				m = indexRecv[ 3*i     ];
				n = indexRecv[ 3*i + 1 ];
				p = indexRecv[ 3*i + 2 ];
				
				m = m - (NX/NPROC)*rank;
				
				if( (m >= NX/NPROC) || (m < 0) ) ERROR_HANDLER(ERROR_CRITICAL, "You are asking for a flow information this process doesn't have");
				
				  
				Vsend[ 8*i    ] = vin[p   + (NZ+1)*n    + (NZ+1)*(NY+1)*m    ];
				Vsend[ 8*i + 1] = vin[p+1 + (NZ+1)*n    + (NZ+1)*(NY+1)*m    ];
				Vsend[ 8*i + 2] = vin[p   + (NZ+1)*(n+1)+ (NZ+1)*(NY+1)*m    ];
				Vsend[ 8*i + 3] = vin[p+1 + (NZ+1)*(n+1)+ (NZ+1)*(NY+1)*m    ];
				Vsend[ 8*i + 4] = vin[p   + (NZ+1)*n    + (NZ+1)*(NY+1)*(m+1)];
				Vsend[ 8*i + 5] = vin[p+1 + (NZ+1)*n    + (NZ+1)*(NY+1)*(m+1)];
				Vsend[ 8*i + 6] = vin[p   + (NZ+1)*(n+1)+ (NZ+1)*(NY+1)*(m+1)];
				Vsend[ 8*i + 7] = vin[p+1 + (NZ+1)*(n+1)+ (NZ+1)*(NY+1)*(m+1)];
				
				if((m==NX/NPROC-1)&&(rank==0)) {
					//printf("Vsend extreme is=%g although other is %g\n",Vsend[ 8*i + 7],Vsend[ 8*i + 3]);
				}
				
				//printf("Vsend[8*%d]=%g\n",i,Vsend[ 8*i    ]);
	
		}
		
		// Exchange back the velocity information
		
		MPI_Sendrecv( Vsend, 8*nrecv, MPI_DOUBLE, target, 1,
					  Vrecv, 8*nsend, MPI_DOUBLE, target, 1,
					  MPI_COMM_WORLD, &status);
					  
		// Compute the interpolated field in Varray
		
		for( j = 0 ; j < nsend ; j++) {
			i = indexPos[ j ];
			
			q000 = Vrecv[ 8*j   ];
			q001 = Vrecv[ 8*j+1 ];
			q010 = Vrecv[ 8*j+2 ];
			q011 = Vrecv[ 8*j+3 ];
			q100 = Vrecv[ 8*j+4 ];
			q101 = Vrecv[ 8*j+5 ];
			q110 = Vrecv[ 8*j+6 ];
			q111 = Vrecv[ 8*j+7 ];
#else
		// non MPI routine starts here
		
		DEBUG_START_FUNC;
		
#ifdef _OPENMP
#pragma omp parallel for private(i, m, n, p, dx, dy ,dz, q000, q001, q010, q011, q100, q101, q110, q111) schedule(static)	
#endif
		for( i = 0 ; i < param.particles_n ; i++) {
			m = indexArray[ 3*i     ];
			n = indexArray[ 3*i + 1 ];
			p = indexArray[ 3*i + 2 ];
		
			q000 = vin[p   + (NZ+1)*n    + (NZ+1)*(NY+1)*m    ];
			q001 = vin[p+1 + (NZ+1)*n    + (NZ+1)*(NY+1)*m    ];
			q010 = vin[p   + (NZ+1)*(n+1)+ (NZ+1)*(NY+1)*m    ];
			q011 = vin[p+1 + (NZ+1)*(n+1)+ (NZ+1)*(NY+1)*m    ];
			q100 = vin[p   + (NZ+1)*n    + (NZ+1)*(NY+1)*(m+1)];
			q101 = vin[p+1 + (NZ+1)*n    + (NZ+1)*(NY+1)*(m+1)];
			q110 = vin[p   + (NZ+1)*(n+1)+ (NZ+1)*(NY+1)*(m+1)];
			q111 = vin[p+1 + (NZ+1)*(n+1)+ (NZ+1)*(NY+1)*(m+1)];
#endif
	
					
			// Compute a linear interpolation of the flow at the particle location
			// NB: we add param.lx/y/z fldi.part[i].x/y/z so that the result returned by fmod is always positive.
			dx=fmod( fldi.part[i].x+param.lx , param.lx/((double)NX)) * NX / param.lx;
			dy=fmod( fldi.part[i].y+param.ly , param.ly/((double)NY)) * NY / param.ly;
			dz=fmod( fldi.part[i].z+param.lz , param.lz/((double)NZ)) * NZ / param.lz;

			Varray[i] =	  (1.0-dx) * (1.0-dy)*(1.0-dz) * q000
						+ (1.0-dx) * (1.0-dy)*(    dz) * q001
						+ (1.0-dx) * (    dy)*(1.0-dz) * q010
						+ (1.0-dx) * (    dy)*(    dz) * q011
						+ (    dx) * (1.0-dy)*(1.0-dz) * q100
						+ (    dx) * (1.0-dy)*(    dz) * q101
						+ (    dx) * (    dy)*(1.0-dz) * q110
						+ (    dx) * (    dy)*(    dz) * q111;

		}
#ifdef MPI_SUPPORT
	}
#endif
	
	
	DEBUG_END_FUNC;
	
	return;
}
			
// Map grid based scalar to a particle based scalar
void particles_to_flow(int *indexArray, struct Field fldi, double *Varray, double *vout) {
	int m,n,p;
	int i,j;
	
	double dx, dy, dz;
	double q000, q001, q010, q011, q100, q101, q110, q111;
	
#ifdef MPI_SUPPORT
	int nstep;
	int target;
	int nsend;
	int nrecv;
	
	MPI_Status	status;
		
	static int *indexSend = NULL;
	static int *indexRecv = NULL;
	static int *indexPos = NULL;
	static double *Vsend = NULL;
	static double *Vrecv = NULL;
	
	DEBUG_START_FUNC;
	
	// Allocate memory if it is the first time this function is called
	
	if(indexSend == NULL) {
		indexSend = (int *) malloc( 3 * sizeof(int) * param.particles_n/NPROC );
		if(indexSend == NULL) ERROR_HANDLER(ERROR_CRITICAL, "No memory for indexSend allocation");
		
		indexRecv = (int *) malloc( 3 * sizeof(int) * param.particles_n/NPROC );
		if(indexRecv == NULL) ERROR_HANDLER(ERROR_CRITICAL, "No memory for indexRecv allocation");
		
		indexPos = (int *) malloc(      sizeof(int) * param.particles_n/NPROC );
		if(indexPos == NULL) ERROR_HANDLER(ERROR_CRITICAL, "No memory for indexPos allocation");
		
		Vsend = (double *) malloc( 8 * sizeof(double) * param.particles_n/NPROC );
		if(Vsend == NULL) ERROR_HANDLER(ERROR_CRITICAL, "No memory for Vsend allocation");
		
		Vrecv = (double *) malloc( 8 * sizeof(double) * param.particles_n/NPROC );
		if(Vrecv == NULL) ERROR_HANDLER(ERROR_CRITICAL, "No memory for Vrecv allocation");
	}
	
	// Loop on each process
	for( nstep = 0; nstep < NPROC ; nstep++) {
		
		target = alltoall_schedule( NPROC, rank, nstep);	// Find the target processor
		
		nsend = 0;
		// Find the information we need to send to the target processor
		for( i = 0 ; i < param.particles_n/NPROC ; i++) {
			if( (indexArray[ 3*i ] / (NX/NPROC) ) == target ) {
				//particle i belongs to the flow in processor target
				
				indexPos[ nsend ] = i; 
				indexSend[ 3*nsend ] = indexArray[ 3*i ];
				indexSend[ 3*nsend + 1] = indexArray[ 3*i + 1 ];
				indexSend[ 3*nsend + 2] = indexArray[ 3*i + 2 ];
		
				// Build the velocity we have to send
				// Compute a linear interpolation of the flow at the particle location
				dx=fmod( fldi.part[i].x+param.lx/2 , param.lx/((double)NX)) * NX / param.lx;
				dy=fmod( fldi.part[i].y+param.ly/2 , param.ly/((double)NY)) * NY / param.ly;
				dz=fmod( fldi.part[i].z+param.lz/2 , param.lz/((double)NZ)) * NZ / param.lz;
			
				Vsend[ 8*nsend    ] = (1.0-dx) * (1.0-dy) * (1.0-dz) * Varray[i];
				Vsend[ 8*nsend + 1] = (1.0-dx) * (1.0-dy) * (    dz) * Varray[i];
				Vsend[ 8*nsend + 2] = (1.0-dx) * (    dy) * (1.0-dz) * Varray[i];
				Vsend[ 8*nsend + 3] = (1.0-dx) * (    dy) * (    dz) * Varray[i];
				Vsend[ 8*nsend + 4] = (    dx) * (1.0-dy) * (1.0-dz) * Varray[i];
				Vsend[ 8*nsend + 5] = (    dx) * (1.0-dy) * (    dz) * Varray[i];
				Vsend[ 8*nsend + 6] = (    dx) * (    dy) * (1.0-dz) * Varray[i];
				Vsend[ 8*nsend + 7] = (    dx) * (    dy) * (    dz) * Varray[i];

				nsend++;
			}
		}
		
		// We're going to send nsend flow velocities points to processor target
			
		// Exchange the array sizes
			
		MPI_Sendrecv( &nsend, 1, MPI_INT, target, 1, 
					  &nrecv, 1, MPI_INT, target, 1,
					  MPI_COMM_WORLD, &status );
		
		// Exchange the indices
		MPI_Sendrecv( indexSend, 3*nsend, MPI_INT, target, 1,
					  indexRecv, 3*nrecv, MPI_INT, target, 1,
					  MPI_COMM_WORLD, &status );
		
		// From that point indexRecv contains the indices of the velocity field we will receive from the target
			
				
		// Exchange back the velocity information
		
		MPI_Sendrecv( Vsend, 8*nsend, MPI_DOUBLE, target, 1,
					  Vrecv, 8*nrecv, MPI_DOUBLE, target, 1,
					  MPI_COMM_WORLD, &status);
		
		// Compute the interpolated field in Varray
		
		
		for( i = 0 ; i < nrecv ; i++) {
				m = indexRecv[ 3*i     ];
				n = indexRecv[ 3*i + 1 ];
				p = indexRecv[ 3*i + 2 ];
				
				m = m - (NX/NPROC)*rank;
				
				if( (m >= NX/NPROC) || (m < 0) ) ERROR_HANDLER(ERROR_CRITICAL, "You gave me flow information this process doesn't need");
				
				// Add to vout the result (caution, vout is not initialized!)
				
				vout[p   + (NZ+1)*n    + (NZ+1)*(NY+1)*m    ] += Vrecv[ 8*i     ];
				vout[p+1 + (NZ+1)*n    + (NZ+1)*(NY+1)*m    ] += Vrecv[ 8*i + 1 ];
				vout[p   + (NZ+1)*(n+1)+ (NZ+1)*(NY+1)*m    ] += Vrecv[ 8*i + 2 ];
				vout[p+1 + (NZ+1)*(n+1)+ (NZ+1)*(NY+1)*m    ] += Vrecv[ 8*i + 3 ];
				vout[p   + (NZ+1)*n    + (NZ+1)*(NY+1)*(m+1)] += Vrecv[ 8*i + 4 ];
				vout[p+1 + (NZ+1)*n    + (NZ+1)*(NY+1)*(m+1)] += Vrecv[ 8*i + 5 ];
				vout[p   + (NZ+1)*(n+1)+ (NZ+1)*(NY+1)*(m+1)] += Vrecv[ 8*i + 6 ];
				vout[p+1 + (NZ+1)*(n+1)+ (NZ+1)*(NY+1)*(m+1)] += Vrecv[ 8*i + 7 ];

	
		}

	}
#else			// no MPI

	DEBUG_START_FUNC;
	
// This loop cannot be parallelized with OpneMP due to race conditions appearing in the vout accesses.
// Usage of OMP ATOMIC is strongly inefficient in this particular context (tested).
	for( i = 0 ; i < param.particles_n ; i++) {
		m = indexArray[ 3*i     ];
		n = indexArray[ 3*i + 1 ];
		p = indexArray[ 3*i + 2 ];
			
		// Compute a linear interpolation of the flow at the particle location
		dx=fmod( fldi.part[i].x+param.lx , param.lx/((double)NX)) * NX / param.lx;
		dy=fmod( fldi.part[i].y+param.ly , param.ly/((double)NY)) * NY / param.ly;
		dz=fmod( fldi.part[i].z+param.lz , param.lz/((double)NZ)) * NZ / param.lz;
		
		if(dx>1.0 || dx < 0.0) printf("shit with dx=%g, x=%g\n",dx,fldi.part[i].x);
		if(dy>1.0 || dy < 0.0) printf("shit with dy=%g, y=%g\n",dy,fldi.part[i].y);
		if(dz>1.0 || dz < 0.0) printf("shit with dz=%g\n",dz);

		vout[p   + (NZ+1)*n    + (NZ+1)*(NY+1)*m    ] += (1.0-dx) * (1.0-dy)*(1.0-dz) * Varray[i];
		vout[p+1 + (NZ+1)*n    + (NZ+1)*(NY+1)*m    ] += (1.0-dx) * (1.0-dy)*(    dz) * Varray[i];
		vout[p   + (NZ+1)*(n+1)+ (NZ+1)*(NY+1)*m    ] += (1.0-dx) * (    dy)*(1.0-dz) * Varray[i];
		vout[p+1 + (NZ+1)*(n+1)+ (NZ+1)*(NY+1)*m    ] += (1.0-dx) * (    dy)*(    dz) * Varray[i];
		vout[p   + (NZ+1)*n    + (NZ+1)*(NY+1)*(m+1)] += (    dx) * (1.0-dy)*(1.0-dz) * Varray[i];
		vout[p+1 + (NZ+1)*n    + (NZ+1)*(NY+1)*(m+1)] += (    dx) * (1.0-dy)*(    dz) * Varray[i];
		vout[p   + (NZ+1)*(n+1)+ (NZ+1)*(NY+1)*(m+1)] += (    dx) * (    dy)*(1.0-dz) * Varray[i];
		vout[p+1 + (NZ+1)*(n+1)+ (NZ+1)*(NY+1)*(m+1)] += (    dx) * (    dy)*(    dz) * Varray[i];
	}
#endif

	DEBUG_END_FUNC;
	
	return;
}


/***********************************************************/
/** 
	Computes and adds the drag forces due to the flow
	on the particles using a linear interpolation at subgrid scale.
	
	The velocity field used as input is the one found in the sheared
	frame (such as the one given by a call to gfft_c2r(fldi.qi))
	
	We assume the input velocity arrays have size (NX, NY, NZ+2) 
	without mpi or (NY/NPROC, NX, NZ+2) with mpi
	
	@param dfldo: delta field strucutre (in this case dv_i) to which we have to add the drag force
	@param fldi: current state of the flow
	@param qx: x velocity component of the fluid
	@param qy: y velocity component of the fluid
	@param qz: z velocity component of the fluid
	@param t: current time
	@param dt: current timestep
*/
/***********************************************************/

				
void compute_drag_step(struct Field dfldo,
				   struct Field fldi,
				   double *qx,
				   double *qy,
				   double *qz,
				   const double t,
				   const double tremap,
				   const double dt) {
				   
	int i;
	double *vx;
	double *vy;
	double *vz;

	int *indexArray;
	double *VxArray;
	double *VyArray;
	double *VzArray;
	
	double fxt, fyt, fzt;
	double fxt1, fyt1, fzt1;
	
	vx = (double *) fftw_malloc( (NX/NPROC+1) * (NY+1) * (NZ+1) * sizeof(double ));
	if (vx == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for vx allocation");
	
	vy = (double *) fftw_malloc( (NX/NPROC+1) * (NY+1) * (NZ+1) * sizeof(double ));
	if (vy == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for vy allocation");
	
	vz = (double *) fftw_malloc( (NX/NPROC+1) * (NY+1) * (NZ+1) * sizeof(double ));
	if (vz == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for vz allocation");
	
	indexArray = (int *) malloc( param.particles_n/NPROC * 3 * sizeof(int) );	//3 spatial coordinates per particle
	VxArray = (double *) malloc( param.particles_n/NPROC * sizeof(double) ); // interpolated vx values at particle location
	VyArray = (double *) malloc( param.particles_n/NPROC * sizeof(double) ); // interpolated vy values at particle location
	VzArray = (double *) malloc( param.particles_n/NPROC * sizeof(double) ); // interpolated vz values at particle location
	
	// Create a velocity field we can use
	compute_flow_velocity(qx, qy, qz, vx, vy, vz, t, tremap);

	// Compute the index Array
	build_indexArray(fldi, t, indexArray);
	
	flow_to_particles(indexArray, fldi, vx, VxArray);
	flow_to_particles(indexArray, fldi, vy, VyArray);
	flow_to_particles(indexArray, fldi, vz, VzArray);
	
	fxt=0.0;
	fyt=0.0;
	fzt=0.0;
	
#ifdef _OPENMP
	#pragma omp parallel for private(i) schedule(static)	
#endif
	for( i = 0 ; i < param.particles_n/NPROC ; i++) {
		// Compute the forces in VjArray
		
		VxArray[i] = - (fldi.part[i].vx - VxArray[i]) / fldi.part[i].stime;
		VyArray[i] = - (fldi.part[i].vy - VyArray[i]) / fldi.part[i].stime;
		VzArray[i] = - (fldi.part[i].vz - VzArray[i]) / fldi.part[i].stime;
		
		//printf("m=%d, n=%d, p=%d, velocity on location is (%g, %g, %g)\n",m,n,p, vx[p + (NZ+1)*n + (NZ+1)*(NY+1)*m],vy[p + (NZ+1)*n + (NZ+1)*(NY+1)*m],vz[p + (NZ+1)*n + (NZ+1)*(NY+1)*m]);		
		dfldo.part[i].vx += VxArray[i];
		dfldo.part[i].vy += VyArray[i];
		dfldo.part[i].vz += VzArray[i];
	}

	// initialize to 0 vx, vy, vz
	
	for(i=0 ; i < (NX/NPROC+1) * (NY+1) * (NZ+1) ; i++) {
		vx[i]=0.0;
		vy[i]=0.0;
		vz[i]=0.0;
	}
	
	// Map back forces fo the grid.
	particles_to_flow(indexArray, fldi, VxArray, vx);
	particles_to_flow(indexArray, fldi, VyArray, vy);
	particles_to_flow(indexArray, fldi, VzArray, vz);
	

	// remove ghost zone and remap
	map_flow_forces(vx,wr4,t,tremap);
	map_flow_forces(vy,wr5,t,tremap);
	map_flow_forces(vz,wr6,t,tremap);
	
	
	// fourier transform
	gfft_r2c(wr4);
	gfft_r2c(wr5);
	gfft_r2c(wr6);
	
		// Add to the deviations (including antialiasing)
	
#ifdef _OPENMP
	#pragma omp parallel for private(i) schedule(static)	
#endif
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		dfldo.vx[i] -= param.particles_dg_ratio * NX * NY * NZ * mask[i] * w4[i] / param.particles_n;
		dfldo.vy[i] -= param.particles_dg_ratio * NX * NY * NZ * mask[i] * w5[i] / param.particles_n;
		dfldo.vz[i] -= param.particles_dg_ratio * NX * NY * NZ * mask[i] * w6[i] / param.particles_n;
	}


	fftw_free(vx);
	fftw_free(vy);
	fftw_free(vz);	
	free(indexArray);
	free(VxArray);
	free(VyArray);
	free(VzArray);

	return;
		
}

/***********************************************************/
/** 
	This routine computes the forces and displacements
	of the particles. It is supposed to be called from the mainloop.
	As an optimization, it takes as an argument the (real) velocity
	field components of the flow. This is useful since this routine
	can be called from timestep(...) where these fields are already computed
	for the navier-stokes equation.
	
	The velocity field used as input is the one found in the sheared
	frame (such as the one given by a call to gfft_c2r(fldi.qi))
	
	We assume the input velocity arrays have size (NX, NY, NZ+2) 
	without mpi or (NY/NPROC, NX, NZ+2) with mpi

	@param dfldo: delta field strucutre (in this case dv_i) to which we have to add the drag force
	@param fldi: current state of the flow
	@param qx: x velocity component of the fluid
	@param qy: y velocity component of the fluid
	@param qz: z velocity component of the fluid
	@param t: current time
	@param dt: current timestep


*/
/***********************************************************/
void particle_step(struct Field dfldo,
				   struct Field fldi,
				   double *qx,
				   double *qy,
				   double *qz,
				   const double t,
				   const double tremap,
				   const double dt) {
			
	int i;
		
	DEBUG_START_FUNC;
	
	// Timestep for the particle dynamics
	// Unless we include an explicit back-reaction, we should only modify 
	// dfldo.part, and nothing else!!
	
	// Update the position according to the velocity field
	
#ifdef _OPENMP
	#pragma omp parallel for private(i) schedule(static)	
#endif
	for( i = 0 ; i < param.particles_n/NPROC ; i++) {
		dfldo.part[i].vx =   2.0 * param.omega * fldi.part[i].vy - param.particles_epsilon;
#ifdef WITH_SHEAR
		dfldo.part[i].vy = - (2.0 * param.omega - param.shear) * fldi.part[i].vx;
#else
		dfldo.part[i].vy = - 2.0 * param.omega * fldi.part[i].vx;
#endif
		dfldo.part[i].vz = 0.0;
		
		dfldo.part[i].x = fldi.part[i].vx;
#ifdef WITH_SHEAR
		dfldo.part[i].y = fldi.part[i].vy - param.shear * fldi.part[i].x;
#else
		dfldo.part[i].y = fldi.part[i].vy;
#endif
		
		dfldo.part[i].z = fldi.part[i].vz;
	}
		
	compute_drag_step(dfldo, fldi, qx, qy, qz, t, tremap, dt);
	
	DEBUG_END_FUNC;
}
	
/***********************************************************/
/** 
	Move the particles back in the simulation box if they moved out.
	This routine is intended to be called after the full explicit
	step of mainloop, like the implict step for the flow.
	
	@param fldi: current state to be modified
	@param t: current time
	@param dt: timestep (unused here, just for consistancy with the other flow implicitstep).
*/
/***********************************************************/
void particle_implicit_step(struct Field fldi,
						    const double t,
							const double dt) {
	double x0;
	int i;
	
	DEBUG_START_FUNC;
	// Implicit step for particles, actually just a routine to keep the particles in the simulation "box"
	
	for(i = 0 ; i < param.particles_n/NPROC ; i++) {
		x0 = fldi.part[i].x;
		fldi.part[i].x = fldi.part[i].x - param.lx * floor( fldi.part[i].x / param.lx + 0.5 );
		
#ifdef WITH_SHEAR
		fldi.part[i].y = fldi.part[i].y + (x0 - fldi.part[i].x) * param.shear * t;
#endif

		fldi.part[i].y = fldi.part[i].y - param.ly * floor( fldi.part[i].y / param.ly + 0.5 );

		fldi.part[i].z = fldi.part[i].z - param.lz * floor( fldi.part[i].z / param.lz + 0.5 );

		//printf("t=%g, Particle %d: x=%g, y=%g, z=%g\n",t+dt, i, fldi.part[i].vx, fldi.part[i].vy, fldi.part[i].vz);
	}
	
	DEBUG_END_FUNC;
}


/***********************************************************/
/** 
	Clean the memory and the allocation due to fft routines.
*/
/***********************************************************/


void finish_particles() {

	DEBUG_START_FUNC;
	
	free( x3D );
	free( y3D );
	free( z3D );
	
#ifdef WITH_SHEAR
	
	fftw_destroy_plan(fft_particle_forward);
	fftw_destroy_plan(fft_particle_backward);	
#endif


	DEBUG_END_FUNC;
	return;
}

#endif