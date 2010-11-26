#ifndef __PARTICLES_H
#define __PARTICLES_H

#ifdef WITH_PARTICLES
void output_partvart(struct Particle *part, double t);
void output_particles(struct Field fldi, const int n, double t);
void write_vtk_particles(struct Field fldi, FILE * ht, const double t);
void read_particle_dump(FILE *ht, struct Particle *part);
void write_particle_dump(FILE *ht, struct Particle *part);

void init_particles(struct Field fldi);
void particle_step(struct Field dfldo,
				   struct Field fldi,
				   double *vx,
				   double *vy,
				   double *vz,
				   const double t,
				   const double tremap,
				   const double dt);

void particle_implicit_step(struct Field fldi,
						    const double t,
							const double dt);
							
void finish_particles();
#endif
#endif