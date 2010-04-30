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
    along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
*/


#ifndef _GVARS_
#define _GVARS_

#define		NX				96			/**< X Dimension in real space. Must be multiples of NPROC when using MPI.*/
#define		NY				96			/**< Y Dimension in real space. Must be multiples of NPROC when using MPI.*/
#define		NZ				96			/**< Z Dimension in real space. */

//#define		MHD						/**< Uncomment to activate MHD*/

#define		BOUSSINESQ				/**< Uncomment to activate Boussinesq */
#define		VERTSTRAT				/**< Vertical stratification. Otherwise, Boussinesq stratification is in X */
//#define		N2PROFILE				/**< Assumes N2 depends on spatial coordinates. The profile has to be set up by hand in common.c/init_N2_profile(). This feature in in alpha version. */

//#define		WITH_ROTATION			/**< Enable a Coriolis force, assuming the rotation axis is z. The rotation speed is set up in snoopy.cfg */

//#define		WITH_SHEAR				/**< Uncomment to activate mean SHEAR */
//#define		TIME_DEPENDANT_SHEAR	/**< Enable Time dependant shear */

#define			WITH_LINEAR_TIDE		/**< Enable linear tide evolution following Goodman & Oh (1997) */

//#define		ANELASTIC				/**< Anelastic approximation in the x direction, assuming an exponential density profile. This work only for the moment with pure hydro (no MHD/BOUSSINESQ).*/

//#define		WITH_2D					/**< 2D (x,y) reduction. Enforce NZ=1 and optimise the Fourier transform routines accordingly. Using this option should lead to a speedup>4 compared to the old NZ=2 method.*/

//#define		SGS						/**< Chollet-Lesieur Subgrid scale model for Hydro. This is in alpha version (see timestep.c).*/

//#define		WITH_PARTICLES			/**< Enable Lagrangian test particles tracing. */

//#define		BOUNDARY_C				/**< Enforce NON-PERIODIC boundary conditions to the flow using Fourier space symmetries. The precise boundary condition has to be set up in boundary.c/boundary_c(); This feature is in alpha version. */

//#define		FORCING					/**< Uncomment to use internal forcing of the velocity field (see forcing in timestep.c) */

#define		FFT_PLANNING	FFTW_MEASURE  /**< can be either FFTW_ESTIMATE, FFTW_MEASURE, FFTW_PATIENT or FFTW_EXHAUSTIVE (see fftw3 doc). Measure leads to longer initialisation of fft routines */



#endif