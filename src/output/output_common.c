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

/**********************************************************
     Check that the file size is equal to the expected size
	 return 1 if size is different
	 This check is useful for filesystems implementing quotas which do not
	 return an error when the file quota is reached (leading to impossible restarts).
	 @param filename	name of the file one has to check
	 @param filesize	expected file size
*/
/**********************************************************/
	 
int check_file_size(char *filename, long int filesize) {
	FILE* ht;
	int rvalue;
	
	if( !file_exist(filename) ) {
		MPI_Printf("check_file_size: File %s not found in current directory.\n",filename);
		return(-1);
	}
	
	if(rank==0) {
		ht=fopen(filename,"r");
		fseek(ht, 0L, SEEK_END);
		if( filesize != ftell(ht) ) 
			rvalue = ftell(ht);
		else
			rvalue = 0;
		fclose(ht);
	}

#ifdef MPI_SUPPORT
	MPI_Bcast(&rvalue, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

	return(rvalue);
	
}
