
void init_transpose();
void finish_transpose();
double read_transpose_timer();
void transpose_real(const int nxin, const int nyin, const int nzin, const int nproc, double *qin, double *qout);
void transpose_complex_XY(double complex *qin, double complex *qout);
void transpose_complex_YX(double complex *qin, double complex *qout);
void transpose_complex_ZY(double complex *qin, double complex *qout);
void transpose_complex_YZ(double complex *qin, double complex *qout);