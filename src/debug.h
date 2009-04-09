#ifdef DEBUG

void D_show_all(struct Field fldi);
void D_show_field(PRECISION complex * field);
void debug_start_f(const char ErrorRoutine[], const int line, const char ErrorFile[]);
void debug_end_f(const char ErrorRoutine[], const int line, const char ErrorFile[]);

#define		DEBUG_START_FUNC		debug_start_f(__func__, __LINE__, __FILE__)
#define		DEBUG_END_FUNC			debug_end_f(__func__, __LINE__, __FILE__)

#else

#define		DEBUG_START_FUNC
#define		DEBUG_END_FUNC

#endif
