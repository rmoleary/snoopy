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
#include <string.h>
#include "common.h"
#include "debug.h"
#include "libconfig/libconfig.h"

#define CONFIG_FILENAME		"snoopy.cfg"

void read_config() {
	// Read the config file and initialize everyting
	config_t	config;		// Initialize the structure
	config_setting_t * setting;	// a setting structure
	long tmp_v;
	int i,n;
	
	const char * configname;
	
	const char * temp_string;
	
	
	DEBUG_START_FUNC;
	
	if(rank==0) {
		config_init(&config);
	
		if(!config_read_file(&config, CONFIG_FILENAME)) {
			MPI_Printf("Error reading configuration file in line %d: %s\n", config_error_line(&config), config_error_text(&config));
			ERROR_HANDLER(ERROR_CRITICAL, "Failed to read the configuration file");
		}

		if(config_lookup_string(&config,"configname",&configname)) {
			MPI_Printf("Using config file: %s.\n",configname);
		}
		// read physics parameters-------------------------------------------------------------------------------
		if(!config_lookup_float(&config, "physics.boxsize.[0]",&param.lx)) {
			param.lx = 1.0;
		}
		if(!config_lookup_float(&config, "physics.boxsize.[1]",&param.ly)) {
		param.ly = 1.0;
		}
		if(!config_lookup_float(&config, "physics.boxsize.[2]",&param.lz)) {
			param.lz = 1.0;
		}
		
		if(!config_lookup_float(&config, "physics.reynolds",&param.reynolds)) {
			param.reynolds = 1.0;
		}	
	
		if(!config_lookup_float(&config, "physics.reynolds_magnetic",&param.reynolds_m)) {
			param.reynolds_m = 1.0;
		}
	
		if(!config_lookup_float(&config, "physics.reynolds_thermic",&param.reynolds_th)) {
			param.reynolds_th = 1.0;
		}
	
		if(!config_lookup_float(&config, "physics.reynolds_Braginskii",&param.reynolds_B)) {
			param.reynolds_B = 1.0;
		}
		
		if(!config_lookup_float(&config, "physics.brunt_vaissala_squared",&param.N2)) {
			param.N2 = 0.0;
		}
	
		if(!config_lookup_float(&config, "physics.omega",&param.omega)) {
			param.omega = 0.0;
		}
#ifndef WITH_ROTATION
		// Omega should be forced to zero in order to be fool-proof
		param.omega = 0.0;
#endif
		if(!config_lookup_float(&config, "physics.shear",&param.shear)) {
			param.shear = 0.0;
		}
#ifndef WITH_SHEAR
		// same for the shear
		param.shear = 0.0;
#endif	
		if(!config_lookup_float(&config, "physics.omega_shear",&param.omega_shear)) {
			param.omega_shear = 0.0;
		}
	
		// Particles parameters-------------------------------------------------------------------------------------
		if(!config_lookup_int(&config, "particles.n",&tmp_v)) {
			param.particles_n = 1000;
		}
		else {
			param.particles_n = (int) tmp_v;
		}
		
		if(!config_lookup_float(&config, "particles.mass",&param.particles_mass)) {
			param.particles_mass = 1.0;
		}
		
		if(!config_lookup_float(&config, "particles.stime",&param.particles_stime)) {
			param.particles_stime = 1.0;
		}
		
		if(!config_lookup_float(&config, "particles.dg_ratio",&param.particles_dg_ratio)) {
			param.particles_dg_ratio = 0.01;
		}
		
		if(!config_lookup_float(&config, "particles.epsilon",&param.particles_epsilon)) {
			param.particles_epsilon = 0.1;
		}
		
		// Code parameters-------------------------------------------------------------------------------------
	
		if(!config_lookup_float(&config, "code.cfl",&param.cfl)) {
			param.cfl = 1.5;
		}
		if(!config_lookup_float(&config, "code.safety_source",&param.safety_source)) {
			param.safety_source = 0.2;
		}
		if(!config_lookup_float(&config, "code.t_initial",&param.t_initial)) {
			param.t_initial = 0.0;
		}
		if(!config_lookup_float(&config, "code.t_final",&param.t_final)) {
			param.t_final = 1.0;
		}
		if(!config_lookup_float(&config, "code.max_t_elapsed",&param.max_t_elapsed)) {
			param.max_t_elapsed = 1e30;
		}
		if(!config_lookup_int(&config, "code.interface_check",&tmp_v)) {
			param.interface_check = 5;
		}
		else {
			param.interface_check = (int) tmp_v;
		}
		if(!config_lookup_bool(&config, "code.interface_output_file",&param.interface_output_file)) {
			param.interface_output_file = 0;
		}
		if(!config_lookup_bool(&config, "code.force_symmetries",&param.force_symmetries)) {
			param.force_symmetries = 0;
		}
		if(!config_lookup_int(&config, "code.symmetries_step",&tmp_v)) {
			param.symmetries_step = 20;
		}
		else {
			param.symmetries_step = (int) tmp_v;
		}
		if(!config_lookup_bool(&config, "code.antialiasing",&param.antialiasing)) {
			param.antialiasing = 1;
		}
		if(!config_lookup_bool(&config, "code.restart",&param.restart)) {
			param.restart = 0;
		}

		// Output parameters-------------------------------------------------------------------------------------
		if(!config_lookup_float(&config, "output.timevar_step",&param.toutput_time)) {
			param.toutput_time = 1.0;
		}
		if(!config_lookup_float(&config, "output.snapshot_step",&param.toutput_flow)) {
			param.toutput_flow = 1.0;
		}
		if(!config_lookup_float(&config, "output.dump_step",&param.toutput_dump)) {
			param.toutput_dump = 1.0;
		}
		if(!config_lookup_bool(&config, "output.vorticity",&param.output_vorticity)) {
			param.output_vorticity = 0;
		}
		// find which parameters are requested in the timevar file
		setting = config_lookup(&config, "output.timevar_vars");
		
		if(setting == NULL) {
			ERROR_HANDLER(ERROR_WARNING, "You did not provide any variable in timevar outputs");
		}
		else {
			param.timevar_vars.length = config_setting_length( setting );
		
			// Allocate output_vars
			param.timevar_vars.name = malloc( param.timevar_vars.length * sizeof(char*) );
		
			for(i = 0 ; i < param.timevar_vars.length ; i++) {
				temp_string = config_setting_get_string_elem( setting, i);
			
				// Allocate the string
				param.timevar_vars.name[i] = malloc( sizeof(char) * (strlen(temp_string) + 1));
			
				// Copy the string in the right location
				strcpy(param.timevar_vars.name[i], temp_string);
			}
		}
		
		// find which parameters are requested in the snapshot files
		setting = config_lookup(&config, "output.snapshot_vars");
		
		if(setting == NULL) {
			ERROR_HANDLER(ERROR_WARNING, "You did not provide any variable in snapshots outputs");
		}
		else {
			param.snapshot_vars.length = config_setting_length( setting );
		
			// Allocate output_vars
			param.snapshot_vars.name = malloc( param.snapshot_vars.length * sizeof(char*) );
		
			for(i = 0 ; i < param.snapshot_vars.length ; i++) {
				temp_string = config_setting_get_string_elem( setting, i);
			
				// Allocate the string
				param.snapshot_vars.name[i] = malloc( sizeof(char) * (strlen(temp_string) + 1));
			
				// Copy the string in the right location
				strcpy(param.snapshot_vars.name[i], temp_string);
			}
		}

		
		// Initial conditions parameters-------------------------------------------------------------------------
		if(!config_lookup_bool(&config, "init.vortex.enable",&param.init_vortex)) {
			param.init_vortex = 0;
		}
		if(!config_lookup_float(&config, "init.vortex.a",&param.vortex_a)) {
			param.vortex_a = 1.0;
		}
		if(!config_lookup_float(&config, "init.vortex.b",&param.vortex_b)) {
			param.vortex_b = 2.0;
		}
		if(!config_lookup_bool(&config, "init.spatial_structure",&param.init_spatial_structure)) {
			param.init_spatial_structure = 0;
		}
		if(!config_lookup_bool(&config, "init.large_scale_noise.enable",&param.init_large_scale_noise)) {
			param.init_large_scale_noise = 0;
		}	
		if(!config_lookup_float(&config, "init.large_scale_noise.amplitude",&param.per_amplitude_large)) {
			param.per_amplitude_large = 0.0;
		}	
		if(!config_lookup_float(&config, "init.large_scale_noise.cut_length",&param.noise_cut_length)) {
			param.noise_cut_length = 0.0;
		}
		if(!config_lookup_bool(&config, "init.large_scale_2D_noise.enable",&param.init_large_scale_2D_noise)) {
			param.init_large_scale_2D_noise = 0;
		}	
		if(!config_lookup_float(&config, "init.large_scale_2D_noise.amplitude",&param.per_amplitude_large_2D)) {
			param.per_amplitude_large_2D = 0.0;
		}	
		if(!config_lookup_float(&config, "init.large_scale_2D_noise.cut_length",&param.noise_cut_length_2D)) {
			param.noise_cut_length_2D = 0.0;
		}
		if(!config_lookup_bool(&config, "init.white_noise.enable",&param.init_white_noise)) {
			param.init_white_noise = 0;
		}
		if(!config_lookup_float(&config, "init.white_noise.amplitude",&param.per_amplitude_noise)) {
			param.per_amplitude_noise = 0.0;
		}
		if(!config_lookup_bool(&config, "init.mean_field.enable",&param.init_mean_field)) {
			param.init_mean_field = 0;
		}
		if(!config_lookup_float(&config, "init.mean_field.bx0",&param.bx0)) {
			param.bx0 = 0.0;
		}
		if(!config_lookup_float(&config, "init.mean_field.by0",&param.by0)) {
			param.by0 = 0.0;
		}
		if(!config_lookup_float(&config, "init.mean_field.bz0",&param.bz0)) {
			param.bz0 = 0.0;
		}
		if(!config_lookup_bool(&config, "init.dump",&param.init_dump)) {
			param.init_dump = 0;
		}
		if(!config_lookup_bool(&config, "init.bench",&param.init_bench)) {
			param.init_bench = 0;
		}
		config_destroy(&config);
	}
#ifdef MPI_SUPPORT
	MPI_Bcast( &param, sizeof(struct Parameters), MPI_CHAR, 0, MPI_COMM_WORLD);
	
	// Copy varname structures properly (Broadcast does not work because of the allocation structure we use)
	if(rank !=0 ) {
		// Allocate the name list
		param.timevar_vars.name = malloc( param.timevar_vars.length * sizeof(char*) );
	}
	
	// Next, allocate each name and copy it
	for(i = 0 ; i < param.timevar_vars.length ; i++) {
		if(rank==0) n = strlen(param.timevar_vars.name[i]);
		
		// Broadcast the string length and allocate it
		MPI_Bcast( &n, 1, MPI_INT, 0, MPI_COMM_WORLD);
		
		if(rank != 0) param.timevar_vars.name[i] = malloc( sizeof(char) * (n + 1));
		
		// Broadcast the string itself
		MPI_Bcast( param.timevar_vars.name[i], n+1, MPI_CHAR, 0, MPI_COMM_WORLD);
		
	}
	
	// Same for snapshot varname structure
	if(rank !=0 ) {
		// Allocate the name list
		param.snapshot_vars.name = malloc( param.snapshot_vars.length * sizeof(char*) );
	}
	
	// Next, allocate each name and copy it

	for(i = 0 ; i < param.snapshot_vars.length ; i++) {
		if(rank==0) n = strlen(param.snapshot_vars.name[i]);
		
		// Broadcast the string length and allocate it
		MPI_Bcast( &n, 1, MPI_INT, 0, MPI_COMM_WORLD);
		
		if(rank != 0) param.snapshot_vars.name[i] = malloc( sizeof(char) * (n + 1));
		
		// Broadcast the string itself
		MPI_Bcast( param.snapshot_vars.name[i], n+1, MPI_CHAR, 0, MPI_COMM_WORLD);
	}

#endif
		
	DEBUG_END_FUNC;
	
	return;
}
	
