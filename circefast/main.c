#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "omp.h"
#include "gsl/gsl_sort.h"
#include "gsl/gsl_statistics.h"

#include "fileio.h"
#include "engine.h"

//#define NUM_THREADS 4


int main(int argc, char* argv[]){

  double time;

  printf("yo\n");
  time = omp_get_wtime();
  printf("yo2\n");
  config tconfig;
  strcpy(tconfig.path, "/scratch1/CIRCE/2017-05-23/");
  strcpy(tconfig.f_dark, "darks/dark_band_0341_1706_02_005_04869.fits");
  strcpy(tconfig.f_badp, "darks/badpix_band_0341_1706_02_005_04869.fits");
  strcpy(tconfig.f_flat, "darks/flat_34_36.fits");
  //strcpy(tconfig.f_flat, "darks/sky_flat_lin.fits");
  strcpy(tconfig.f_name, "CIRCE2017-05-23-%04d.fits");
  //tconfig.n_dith = 5;
  //tconfig.f_arr[0] = 79;
  //tconfig.f_arr[1] = 84;
  tconfig.seq = atoi(argv[1]);
  tconfig.n_dith = 5;
  tconfig.st = atoi(argv[2]);
  tconfig.band_ndx[0] = 341;
  tconfig.band_ndx[1] = 1707;

  //omp_set_nested(1);
  
  printf("***seqs %d\n", tconfig.seq);

  darksub_all(tconfig);
  doskies_all(tconfig);
  skysub_all(tconfig);
  fftcorr_all(tconfig);	 
  flatdivide_all(tconfig);
  badpixrem_all(tconfig);
  tile_all(tconfig);

  time -= omp_get_wtime();
  time = - time;
  printf("%.2f ms\n", time*1000);

  return(0);
}
