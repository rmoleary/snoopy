extern double time_1D, time_2D, time_transpose;

void init_gfft() ;
void finish_gfft();

void gfft_r2c_t(double *wrin);
void gfft_c2r_t(double complex *win);
void gfft_r2c(double *wrin);
void gfft_c2r(double complex *win);