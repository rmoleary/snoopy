void output(const double t);
void init_output();
void finish_output();
void read_dump(   struct Field fldo,
				  double *t);
void output_status(FILE * iostream);
void output_immediate(const double t);
void dump_immediate(const double t);
void clear_timevar();