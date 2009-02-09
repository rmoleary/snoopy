/* Error.h: Error handling, useful for debugging purposes */

#define		ERROR_HANDLER(ERROR_TYPE, ERROR_MESSAGE)		error_h(ERROR_TYPE, ERROR_MESSAGE, __func__ , __LINE__, __FILE__)
#define		ERROR_WARNING		1
#define		ERROR_NONCRITICAL	2
#define		ERROR_CRITICAL		3

void error_h (const int ErrorType,
			  const char ErrorMessage[], 
			  const char ErrorRoutine[], 
			  const int line, 
			  const char Filename[] );
