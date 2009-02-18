
void init_transpose();
void transpose_real(const int nxin, const int nyin, const int nzin, const int nproc, PRECISION *qin, PRECISION *qout);
void transpose_complex_XY(PRECISION complex *qin, PRECISION complex *qout);
void transpose_complex_YX(PRECISION complex *qin, PRECISION complex *qout);
void transpose_complex_ZY(PRECISION complex *qin, PRECISION complex *qout);
void transpose_complex_YZ(PRECISION complex *qin, PRECISION complex *qout);