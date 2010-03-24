// A structure for particles.

#ifndef __PARTICLES_H
#define __PARTICLES_H

void output_particles(const int n, double t);
void write_vtk_particles(FILE * ht, const double t);

void init_particles();
void particle_step(struct Field dfldo,
				   struct Field fldi,
				   double *vx,
				   double *vy,
				   double *vz,
				   const double t,
				   const double dt);

void particle_implicit_step(struct Field fldi,
						    const double t,
							const double dt);
							
void finish_particles();
#endif
