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


#include "common.h"
#include "debug.h"
#include "libconfig/libconfig.h"

#define CONFIG_FILENAME		"snoopy.cfg"

read_config() {
	// Read the config file and initialize everyting
	config_t	config;		// Initialize the structure
	long tmp_v;
	
	const char * configname;
	
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
	
		if(!config_lookup_float(&config, "physics.brunt_vaissala_squared",&param.N2)) {
			param.N2 = 0.0;
		}
	
		if(!config_lookup_float(&config, "physics.omega",&param.omega)) {
			param.omega = 0.0;
		}
	
		if(!config_lookup_float(&config, "physics.shear",&param.shear)) {
			param.shear = 0.0;
		}
	
		if(!config_lookup_float(&config, "physics.omega_shear",&param.omega_shear)) {
			param.omega_shear = 0.0;
		}
		
		if(!config_lookup_float(&config, "physics.anelastic_lambda",&param.anelastic_lambda)) {
			param.anelastic_lambda = 1.0;
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
		if(!config_lookup_bool(&config, "output.pressure",&param.output_pressure)) {
			param.output_pressure = 0;
		}
		if(!config_lookup_bool(&config, "output.vorticity",&param.output_vorticity)) {
			param.output_vorticity = 0;
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
#endif
	
	DEBUG_END_FUNC;
	
	return;
}
	
