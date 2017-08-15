#include "fftw3.h"

typedef struct{
  char path[100];
  char d_flat[100];
  char f_flat[100];
  char f_name[100];
  char f_dark[100];
  char f_badp[100];
  char hwp[5];
  char filt[2];
  int n_dith;
  int f_arr[2];
  int band_ndx[2];
  int add_back;
  int fft_coor;
  int seq;
  int ndx;
  int st;
}config;

int darksub(config tconfig, fitsobj *dark, int ndx, int st);
int doskies(config tconfig, int ndx, int st);
int skysub(config tconfig, fitsobj *sky, int ndx, int st);
int flatdivide(config tconfig, fitsobj *flat, int ndx, int st);
int badpixrem(config tconfig, fitsobj *badpix, int ndx, int st);
int fftcorr(config tconfig, int ndx, int st);
int tile_images(config tconfig, int ndx, int st);
int darksub_all(config tconfig);
int doskies_all(config tconfig);
int skysub_all(config tconfig);
int flatdivide_all(config tconfig);
int badpixrem_all(config tconfig);
int fftcorr_all(config tconfig);
int tile_all(config tconfig);
int badpixfun(float *image, float *badpixim, int naxis1, int naxis2); 
int fftcorrfun(float* image, int naxis1, int naxis2, fftw_plan pfor,
	       fftw_plan pback, fftw_complex*, fftw_complex*);
