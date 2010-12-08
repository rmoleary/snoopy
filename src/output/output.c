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

#include "../common.h"
#include "../gfft.h"
#include "../shear.h"
#include "../debug.h"
#include "../particles.h"
#include "output_vtk.h"
#include "output_dump.h"
#include "output_timevar.h"
#include "output_spectrum.h"
#include "output_common.h"


int noutput_flow;										/**< Next snapshot output number */
double lastoutput_time;								/**< Time when the last timevar output was done */
double lastoutput_flow;								/**< Time when the las snapshot output was done */
double lastoutput_dump;								/**< Time when the last dump output was done */

double output_timer;

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
	
	output_timer=0.0;
	
	init_output_vtk();

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
		init_timevar();
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
	
	finish_output_vtk();
	
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

void output(struct Field fldi, const double t) {
	
	DEBUG_START_FUNC;
	
	// Very rough output function
	if( (t-lastoutput_time)>=param.toutput_time) {
		output_timer = output_timer - get_c_time();
		
		output_timevar(fldi,t);
		output1Dspectrum(fldi,t);
#ifdef WITH_PARTICLES
		output_partvar(fldi.part,t);
#endif
		lastoutput_time = lastoutput_time + param.toutput_time;
		
		output_timer = output_timer + get_c_time();
	}
	
	if( (t-lastoutput_flow)>=param.toutput_flow) {
		output_timer = output_timer - get_c_time();
		
		output_vtk(fldi,noutput_flow,t);
#ifdef WITH_PARTICLES
		output_particles(fldi,noutput_flow,t);
#endif
		noutput_flow++;
		lastoutput_flow = lastoutput_flow + param.toutput_flow;
		
		output_timer = output_timer + get_c_time();
	}
		
	if( (t-lastoutput_dump)>=param.toutput_dump) {
		output_timer = output_timer - get_c_time();
		lastoutput_dump=lastoutput_dump+param.toutput_dump;
		output_dump(fldi,t);
		output_timer = output_timer + get_c_time();
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

void output_immediate(struct Field fldi, const double t) {
	// Very rough output function
	// Immediate output
	output_timevar(fldi,t);
	output_vtk(fldi, noutput_flow,t);
	noutput_flow++;
	output_dump(fldi,t);
	return;
}

/**************************************************************************************/
/** 
	Immediatly output a dump file, regardless of the output
	parameters.
	@param t Current time of the simulation.
*/
/**************************************************************************************/

void dump_immediate(struct Field fldi, const double t) {
	output_dump(fldi,t);
	return;
}

double read_output_timer() {
	return(output_timer);
}

