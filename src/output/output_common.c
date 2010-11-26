#include <stdlib.h>
#include "../common.h"

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
