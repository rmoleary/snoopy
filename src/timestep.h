void timestep( struct Field dfldo,
			   struct Field fldi,
			   double complex *po,
			   const double t,
			   const double dt);
			   			   
void implicitstep(
			   struct Field fldi,
			   const double t,
			   const double dt );