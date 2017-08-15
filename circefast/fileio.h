

typedef struct {
  float *data;
  float *data2;
  char f_name[100];
  int nramps;
  int ngroups;
  int naxis1;
  int naxis2;
  int status;
} fitsobj;

int fits_info(fitsobj*, char*);
int read_single_fits(char *f_name, fitsobj* obj);
int read_fits(char *f_name, fitsobj* obj);
int read_fits_omp(char *f_name, fitsobj* obj);
int read_diff_fits(char *f_name, fitsobj* obj);
int read_diff_framp(char *f_name, fitsobj* obj);
int write_single_fits(char *f_name, fitsobj *obj);
int write_fits(char *f_name, fitsobj *obj);

//#include "fileio.h"

/*
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

int darksub_all(config tconfig);

int doskies_all(config tconfig);

int skysub_all(config tconfig);

int flatdivide_all(config tconfig);
*/
