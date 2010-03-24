// routines for particle dynamic

#include <stdlib.h>

#include "common.h"
#include "gfft.h"
#include "debug.h"
#include "shear.h"

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
					const double t) {
					
	int i,j,k;
	double tvelocity;
	double tremap;
	complex double wexp;
	complex double phase;
	double complex		*wremap;						/** 1D arrays used by the remap methods */
	
	DEBUG_START_FUNC;
	
	tremap = time_shift(t);
	tvelocity = fmod(t, 2.0 * param.ly / (param.shear * param.lx));
	
	wremap = (double complex *) fftw_malloc( sizeof(double complex) * (NY/2+1) * (NZ+1) );
	if (wremap == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for wremap allocation");
	
	for( i = 0 ; i < NX+1 ; i++) {
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


#endif

/***********************************************************/
/** 
	Output the particle positions and mass to a VTK file
	using a POLYDATA format. The particle position can then be 
	easily represented in paraview using glyphs.
	
	@param n file number to be created
	@param t Current time of the simulation.
	
*/
/***********************************************************/

void output_particles(const int n, double t) {

	FILE *ht = NULL;
	char  filename[50];
	int num_remain_field;
	int array_size, i;
	float q0;
	
	DEBUG_START_FUNC;

	sprintf(filename,"data/p%04i.vtk",n);

	if(rank==0) {
		ht=fopen(filename,"w");
	
		fprintf(ht, "# vtk DataFile Version 2.0\n");
		fprintf(ht, "t= %015.15e Snoopy Code v5.0\n",t);
		fprintf(ht, "BINARY\n");
		fprintf(ht, "DATASET POLYDATA\n");
		fprintf(ht, "POINTS %d float\n",NPARTICLES);
	}
	// Write Particle position...
	
	for(i = 0 ; i < NPARTICLES ; i++) {
		q0 = big_endian( (float) fld.part[i].x);
		fwrite( &q0, sizeof(float), 1, ht);
		
		q0 = big_endian( (float) fld.part[i].y);
		fwrite( &q0, sizeof(float), 1, ht);
		
		q0 = big_endian( (float) fld.part[i].z);
		fwrite( &q0, sizeof(float), 1, ht);
	}
	
	if(rank==0) {
		fprintf(ht, "POINT_DATA %d\n",NPARTICLES);
		fprintf(ht, "SCALARS %s float\n","mass");
		fprintf(ht, "LOOKUP_TABLE default\n");
	}
		
	for(i = 0 ; i < NPARTICLES ; i++) {
		q0 = big_endian( (float) fld.part[i].mass);
		fwrite( &q0, sizeof(float), 1, ht);
	}
	
	fclose(ht);
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

void write_vtk_particles(FILE * ht, const double t) {
	int i,j,k;
	int m, n, p;
	
	float q0;
	
	for(i = 0 ; i < 2*NTOTAL_COMPLEX ; i++) {
		wr1[i] = 0.0;
	}
	
	for(i = 0 ; i < NPARTICLES ; i++) {
		m = floor( (fld.part[i].x / param.lx + 0.5) * NX );
		n = floor( (fld.part[i].y / param.ly + 0.5) * NY );
		p = floor( (fld.part[i].z / param.lz + 0.5) * NZ );
		
		wr1[ p + n * (NZ + 2) + m * (NZ + 2) * NY]++;
	}
	
	for( k = 0 ; k < NZ; k++) {
		for( j = 0; j < NY; j++) {	
			for( i = 0; i < NX; i++) {
				q0 = big_endian( (float) wr1[k + j * (NZ + 2) + i * NY * (NZ + 2) ] );
				fwrite(&q0, sizeof(float), 1, ht);
			}
		}
	}

	return;
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

void init_particles() {
	int i,j,k;
#ifdef WITH_SHEAR
	double complex		*wremap;
	double *vx;
#endif
	
	const int n_size1D[1] = {NY};
	
	DEBUG_START_FUNC;
	
	// init a real space grid
	x3D = (double *) malloc( (NX+1) * (NY+1) * (NZ+1) * sizeof(double ));
	y3D = (double *) malloc( (NX+1) * (NY+1) * (NZ+1) * sizeof(double ));
	z3D = (double *) malloc( (NX+1) * (NY+1) * (NZ+1) * sizeof(double ));
	
	for(i = 0 ; i < NX+1 ; i++) {
		for( j = 0 ; j < NY+1 ; j++) {
			for( k = 0 ; k < NZ+1 ; k++) {
				x3D[ k + (NZ+1) * j + (NZ+1) * (NY+1) * i ] = -param.lx / 2.0 + i * param.lx / ((double) NX);
				y3D[ k + (NZ+1) * j + (NZ+1) * (NY+1) * i ] = -param.ly / 2.0 + j * param.ly / ((double) NY);
				z3D[ k + (NZ+1) * j + (NZ+1) * (NY+1) * i ] = -param.lz / 2.0 + k * param.lz / ((double) NZ);
			}
		}
	}
	
	// init particle positions and velocity
	
	for(i = 0 ; i < 100 ; i++) {
		for(j = 0 ; j < 100 ; j++) {
			fld.part[ j + 100 * i ].x = -param.lx / 2.0 + param.lx * i / 100;
			fld.part[ j + 100 * i ].y = -param.ly / 2.0 + param.ly * j / 100;
			fld.part[ j + 100 * i ].z = 0.0;
		}
	}
			
	for(i = 0 ; i < NPARTICLES ; i++) {
		//fld.part[i].x = 0.4-0.1*i;
		//fld.part[i].y = 0.0;
		//fld.part[i].z = 0.0;
		//fld.part[i].x = -param.lx/2.0 +param.lx*randm();
		//fld.part[i].y = -param.ly/2.0 +param.ly*randm();
		//fld.part[i].z = 0.0;
		
		fld.part[i].vx = 0.0;//0.1*randm()-0.5;
		fld.part[i].vy = 0.0;//*randm()-0.5;
		fld.part[i].vz = 0.0;
		fld.part[i].mass = 1.0;
	}
	
	// Init velocity vectors
	
	
	
#ifdef WITH_SHEAR
// Initialize 1D arrays for remaps
	wremap = (double complex *) fftw_malloc( sizeof(double complex) * (NY/2+1) * (NZ+1) );
	if (wremap == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for wremap allocation");

	vx = (double *) fftw_malloc( (NX+1) * (NY+1) * (NZ+1) * sizeof(double ));
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

/***********************************************************/
/** 
	Computes the flow velocity (minus shear if shear is present)
	in the traditional rectangular frame. It takes as an input the velocity
	field in a sheared frame (such as the one given by a call to gfft_c2r(fld.qi))
	and produces an unsheared velocity field to which ghost zones are added in each direction
	
	The ghost zones are used by the interpolation routines later one
	
	We assume the input arrays have size (NX, NY, NZ+2) without mpi or (NY/NPROC, NX, NZ+2) with mpi
	
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
						   const double t) {
	
	int i,j,k;
	
	DEBUG_START_FUNC;
	
	// Reshape the input array into a (NX+1)(NY+1)(NZ+1) array
#ifdef _OPENMP
	#pragma omp parallel for private(i,j,k) schedule(static)	
#endif
	for( i = 0 ; i < NX ; i++) {
		for( j = 0  ; j < NY ; j++) {
			for( k = 0 ; k < NZ ; k++) {
				qxo[k + j * (NZ+1) + i * (NZ+1) * (NY+1)] = qxi[k + j * (NZ + 2) + i * (NZ+2) * NY] / ((double)NTOTAL);
				qyo[k + j * (NZ+1) + i * (NZ+1) * (NY+1)] = qyi[k + j * (NZ + 2) + i * (NZ+2) * NY] / ((double)NTOTAL);
				qzo[k + j * (NZ+1) + i * (NZ+1) * (NY+1)] = qzi[k + j * (NZ + 2) + i * (NZ+2) * NY] / ((double)NTOTAL);
			}
		}
	}

	// Periodize the output array in the x direction i=O->i=NX
	
	for( j = 0  ; j < NY ; j++) {
		for( k = 0 ; k < NZ ; k++) {
			qxo[k + j * (NZ+1) + NX * (NZ+1) * (NY+1)] = qxo[k + j * (NZ + 1)];
			qyo[k + j * (NZ+1) + NX * (NZ+1) * (NY+1)] = qyo[k + j * (NZ + 1)];
			qzo[k + j * (NZ+1) + NX * (NZ+1) * (NY+1)] = qzo[k + j * (NZ + 1)];
		}
	}
	
#ifdef WITH_SHEAR
	// Remap the flow
	
	remap_flow(qxo, t);
	remap_flow(qyo, t);
	remap_flow(qzo, t);
	
#endif

	// Periodize in the y direction j=0->j=NY
	
	for( i = 0  ; i < NX+1 ; i++) {
		for( k = 0 ; k < NZ ; k++) {
			qxo[k + NY * (NZ+1) + i * (NZ+1) * (NY+1)] = qxo[k + i * (NZ+1) * (NY+1)];
			qyo[k + NY * (NZ+1) + i * (NZ+1) * (NY+1)] = qyo[k + i * (NZ+1) * (NY+1)];
			qzo[k + NY * (NZ+1) + i * (NZ+1) * (NY+1)] = qzo[k + i * (NZ+1) * (NY+1)];

		}
	}
	
	// Periodize in the z direction k=0->k=NZ
	for( i = 0  ; i < NX+1 ; i++) {
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

/***********************************************************/
/** 
	Computes and adds the drag forces due to the flow
	on the particles using a linear interpolation at subgrid scale.
	
	The velocity field used as input is the one found in the sheared
	frame (such as the one given by a call to gfft_c2r(fld.qi))
	
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
				   const double dt) {
				   
	int i;
	int m,n,p;
	double *vx;
	double *vy;
	double *vz;
	double x, y, z;
	double dx, dy, dz;
	double q000, q001, q010, q011, q100, q101, q110, q111;
	double partvx, partvy, partvz;
	
	vx = (double *) fftw_malloc( (NX+1) * (NY+1) * (NZ+1) * sizeof(double ));
	if (vx == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for vx allocation");
	
	vy = (double *) fftw_malloc( (NX+1) * (NY+1) * (NZ+1) * sizeof(double ));
	if (vy == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for vy allocation");
	
	vz = (double *) fftw_malloc( (NX+1) * (NY+1) * (NZ+1) * sizeof(double ));
	if (vz == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for vz allocation");
	
	// Create a velocity field we can use
	compute_flow_velocity(qx, qy, qz, vx, vy, vz, t);
	
#ifdef _OPENMP
	#pragma omp parallel for private(i, x, y, z, m, n, p, dx, dy ,dz, q000, q001, q010, q011, q100, q101, q110, q111, partvx, partvy, partvz) schedule(static)	
#endif
	for( i = 0 ; i < NPARTICLES ; i++) {
		// Add drag forces
		// Compute particle position in the box
		x = fldi.part[i].x - param.lx * floor( fld.part[i].x / param.lx + 0.5 );
#ifdef WITH_SHEAR
		y = fldi.part[i].y + (fld.part[i].x - x) * param.shear * t;
#endif
		y = y - param.ly * floor( y / param.ly + 0.5 );
		z = fldi.part[i].z - param.lz * floor( fld.part[i].z / param.lz + 0.5 );

		// particle indices
		m=(int) floor( (x/param.lx+0.5) * NX);
		n=(int) floor( (y/param.ly+0.5) * NY);
		p=(int) floor( (z/param.lz+0.5) * NZ);
		
		// Compute a linear interpolation of the flow at the particle location
		dx = (x - x3D[ p + n*(NZ+1) + m*(NZ+1)*(NY+1) ])/param.lx * ((double)NX);
		dy = (y - y3D[ p + n*(NZ+1) + m*(NZ+1)*(NY+1) ])/param.ly * ((double)NY);
		dz = (z - z3D[ p + n*(NZ+1) + m*(NZ+1)*(NY+1) ])/param.lz * ((double)NZ);
		
		////////////// VX interpolation
		q000 = vx[p   + (NZ+1)*n    + (NZ+1)*(NY+1)*m    ];
		q001 = vx[p+1 + (NZ+1)*n    + (NZ+1)*(NY+1)*m    ];
		q010 = vx[p   + (NZ+1)*(n+1)+ (NZ+1)*(NY+1)*m    ];
		q011 = vx[p+1 + (NZ+1)*(n+1)+ (NZ+1)*(NY+1)*m    ];
		q100 = vx[p   + (NZ+1)*n    + (NZ+1)*(NY+1)*(m+1)];
		q101 = vx[p+1 + (NZ+1)*n    + (NZ+1)*(NY+1)*(m+1)];
		q110 = vx[p   + (NZ+1)*(n+1)+ (NZ+1)*(NY+1)*(m+1)];
		q111 = vx[p+1 + (NZ+1)*(n+1)+ (NZ+1)*(NY+1)*(m+1)];
		
		partvx = (1.0-dx) * (1.0-dy)*(1.0-dz) * q000
			   + (1.0-dx) * (1.0-dy)*(    dz) * q001
			   + (1.0-dx) * (    dy)*(1.0-dz) * q010
			   + (1.0-dx) * (    dy)*(    dz) * q011
			   + (    dx) * (1.0-dy)*(1.0-dz) * q100
			   + (    dx) * (1.0-dy)*(    dz) * q101
			   + (    dx) * (    dy)*(1.0-dz) * q110
			   + (    dx) * (    dy)*(    dz) * q111;
			  
		////////////// VY interpolation
		
		q000 = vy[p   + (NZ+1)*n    + (NZ+1)*(NY+1)*m    ];
		q001 = vy[p+1 + (NZ+1)*n    + (NZ+1)*(NY+1)*m    ];
		q010 = vy[p   + (NZ+1)*(n+1)+ (NZ+1)*(NY+1)*m    ];
		q011 = vy[p+1 + (NZ+1)*(n+1)+ (NZ+1)*(NY+1)*m    ];
		q100 = vy[p   + (NZ+1)*n    + (NZ+1)*(NY+1)*(m+1)];
		q101 = vy[p+1 + (NZ+1)*n    + (NZ+1)*(NY+1)*(m+1)];
		q110 = vy[p   + (NZ+1)*(n+1)+ (NZ+1)*(NY+1)*(m+1)];
		q111 = vy[p+1 + (NZ+1)*(n+1)+ (NZ+1)*(NY+1)*(m+1)];
		
		partvy = (1.0-dx) * (1.0-dy)*(1.0-dz) * q000
			   + (1.0-dx) * (1.0-dy)*(    dz) * q001
			   + (1.0-dx) * (    dy)*(1.0-dz) * q010
			   + (1.0-dx) * (    dy)*(    dz) * q011
			   + (    dx) * (1.0-dy)*(1.0-dz) * q100
			   + (    dx) * (1.0-dy)*(    dz) * q101
			   + (    dx) * (    dy)*(1.0-dz) * q110
			   + (    dx) * (    dy)*(    dz) * q111;
			   
		////////////// VZ interpolation
		
		q000 = vz[p   + (NZ+1)*n    + (NZ+1)*(NY+1)*m    ];
		q001 = vz[p+1 + (NZ+1)*n    + (NZ+1)*(NY+1)*m    ];
		q010 = vz[p   + (NZ+1)*(n+1)+ (NZ+1)*(NY+1)*m    ];
		q011 = vz[p+1 + (NZ+1)*(n+1)+ (NZ+1)*(NY+1)*m    ];
		q100 = vz[p   + (NZ+1)*n    + (NZ+1)*(NY+1)*(m+1)];
		q101 = vz[p+1 + (NZ+1)*n    + (NZ+1)*(NY+1)*(m+1)];
		q110 = vz[p   + (NZ+1)*(n+1)+ (NZ+1)*(NY+1)*(m+1)];
		q111 = vz[p+1 + (NZ+1)*(n+1)+ (NZ+1)*(NY+1)*(m+1)];
		
		partvz = (1.0-dx) * (1.0-dy)*(1.0-dz) * q000
			   + (1.0-dx) * (1.0-dy)*(    dz) * q001
			   + (1.0-dx) * (    dy)*(1.0-dz) * q010
			   + (1.0-dx) * (    dy)*(    dz) * q011
			   + (    dx) * (1.0-dy)*(1.0-dz) * q100
			   + (    dx) * (1.0-dy)*(    dz) * q101
			   + (    dx) * (    dy)*(1.0-dz) * q110
			   + (    dx) * (    dy)*(    dz) * q111;
		
		
		if(p>NZ||p<0) ERROR_HANDLER( ERROR_CRITICAL, "Error with p");
		if(n>NY||n<0) ERROR_HANDLER( ERROR_CRITICAL, "Error with n");
		if(m>NX||m<0) ERROR_HANDLER( ERROR_CRITICAL, "Error with m");
		// Apply forces
		
		//printf("m=%d, n=%d, p=%d, velocity on location is (%g, %g, %g)\n",m,n,p, vx[p + (NZ+1)*n + (NZ+1)*(NY+1)*m],vy[p + (NZ+1)*n + (NZ+1)*(NY+1)*m],vz[p + (NZ+1)*n + (NZ+1)*(NY+1)*m]);
		dfldo.part[i].vx += - (fldi.part[i].vx - partvx);
		dfldo.part[i].vy += - (fldi.part[i].vy - partvy);
		dfldo.part[i].vz += - (fldi.part[i].vz - partvz);
	}

	fftw_free(vx);
	fftw_free(vy);
	fftw_free(vz);
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
	frame (such as the one given by a call to gfft_c2r(fld.qi))
	
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
	for( i = 0 ; i < NPARTICLES ; i++) {
		dfldo.part[i].vx =   2.0 * param.omega * fldi.part[i].vy;
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
		
	compute_drag_step(dfldo, fldi, qx, qy, qz,t, dt);
	
	DEBUG_END_FUNC;
}
	
void particle_implicit_step(struct Field fldi,
						    const double t,
							const double dt) {
	double x0;
	int i;
	
	DEBUG_START_FUNC;
	// Implicit step for particles, actually just a routine to keep the particles in the simulation "box"
	
	for(i = 0 ; i < NPARTICLES ; i++) {
		x0 = fld.part[i].x;
		fld.part[i].x = fld.part[i].x - param.lx * floor( fld.part[i].x / param.lx + 0.5 );
		
#ifdef WITH_SHEAR
		fld.part[i].y = fld.part[i].y + (x0 - fld.part[i].x) * param.shear * t;
#endif

		fld.part[i].y = fld.part[i].y - param.ly * floor( fld.part[i].y / param.ly + 0.5 );
		fld.part[i].z = fld.part[i].z - param.lz * floor( fld.part[i].z / param.lz + 0.5 );

		//printf("t=%g, Particle %d: x=%g, y=%g, z=%g\n",t+dt, i, fld.part[i].vx, fld.part[i].vy, fld.part[i].vz);
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