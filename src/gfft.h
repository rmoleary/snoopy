extern double time_1D, time_2D, time_transpose;

void init_gfft() ;

void gfft_r2c_t(PRECISION *wrin);
void gfft_c2r_t(PRECISION complex *win);
void gfft_r2c_tt(PRECISION *wrin);
void gfft_c2r_tt(PRECISION complex *win);
void gfft_r2c(PRECISION *wrin);
void gfft_c2r(PRECISION complex *win);
void gfft3_r2c_t(PRECISION *wrin1, PRECISION *wrin2, PRECISION *wrin3);
void gfft3_c2r_t(PRECISION complex *win1, PRECISION complex *win2, PRECISION complex *win3);