extern double time_1D, time_2D, time_transpose;

void init_gfft() ;

void gfft_r2c_t(PRECISION *wrin);
void gfft_c2r_t(PRECISION complex *win);
void gfft_r2c(PRECISION *wrin);
void gfft_c2r(PRECISION complex *win);