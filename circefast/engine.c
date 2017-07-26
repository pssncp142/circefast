#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "omp.h"
#include "gsl/gsl_sort.h"
#include "gsl/gsl_statistics.h"
#include "complex.h"
#include "fftw3.h"

#include "fileio.h"
#include "engine.h"

#define NUM_THREADS 5

int darksub(config tconfig, fitsobj *dark, int ndx, int st){  

  fitsobj *data = malloc(sizeof(fitsobj));
  fitsobj *proc = malloc(sizeof(fitsobj));
  char f_name_tmp[100], f_name[100];
  int npixels;
  int i;
  char f_fmt[100] = "darksub/IPA_seq%d_dith_%d.fits";
  double time;

  strcpy(f_name_tmp, tconfig.path);
  strcat(f_name_tmp, tconfig.f_name);
  sprintf(f_name, f_name_tmp, ndx);

  read_diff_fits(f_name, data);
  sprintf(f_name, f_fmt, tconfig.seq, ndx-st+1);  

  proc->naxis1 = data->naxis1;
  proc->naxis2 = data->naxis2;
  proc->nramps = data->nramps;
  proc->ngroups = data->ngroups;

  npixels = proc->naxis1*proc->naxis2*(proc->ngroups-1)*proc->nramps;
  proc->data = malloc(sizeof(float)*npixels);

  time = omp_get_wtime();

  #pragma omp parallel num_threads(1)
  {
    float val;
    float *data_img2 = &data->data2[0];
    float *data_img1 = &data->data[0];
    float *dark_img = &dark->data[0];
    #pragma omp for 
    for (i=0; i<npixels; i++){
      val = data_img2[i]-data_img1[i]-dark_img[i];
      val = 0.9701104*val+
      1.106431e-5*val*val+
      -7.3981e-10*val*val*val+
      2.490898e-14*val*val*val*val;
      proc->data[i] = val;
    }

  }
  

  time -= omp_get_wtime();
  time = - time;
  printf("Thread %d darksub calc in %.2f ms\n", 
	 omp_get_thread_num(), time*1000);
  
  sprintf(f_name, f_fmt, tconfig.seq, ndx-st+1);  
  write_fits(f_name, proc);

  free(data->data);
  free(data->data2);
  free(proc->data);
  free(data);
  free(proc);
  
  return(0);
}

int doskies(config tconfig, int ndx, int st){  

  fitsobj *data = malloc(sizeof(fitsobj));
  fitsobj *proc = malloc(sizeof(fitsobj));
  char  f_name[100];
  int npixels, nimages;
  int i;
  char f_fmt_in[100] = "darksub/IPA_seq%d_dith_%d.fits";
  char f_fmt_out[100] = "skies/IPA_seq%d_dith_%d.fits";
  double time;

  sprintf(f_name, f_fmt_in, tconfig.seq, ndx-st+1);
  read_fits(f_name, data);

  proc->naxis1 = data->naxis1;
  proc->naxis2 = data->naxis2;
  proc->nramps = data->nramps;
  proc->ngroups = data->ngroups;

  npixels = proc->naxis1*proc->naxis2*(proc->ngroups-1)*proc->nramps;
  proc->data = malloc(sizeof(float)*npixels);

  time = omp_get_wtime();

  //printf("here\n.");

  nimages = proc->nramps*(proc->ngroups-1);

  int pix_img = npixels/nimages;

#pragma omp parallel num_threads(1)
  {
    double *tmp_list = malloc(sizeof(double)*nimages);
    int j;
#pragma omp for schedule(dynamic, 2048)
    for (i=0; i<pix_img; i++){
      for (j=0; j<nimages; j++){
	tmp_list[j] = (double)data->data[i+j*pix_img];
      }
      gsl_sort(tmp_list, 1, nimages);
      proc->data[i] = (float) 
	gsl_stats_median_from_sorted_data(tmp_list, 1, nimages);
    }
  }


  time -= omp_get_wtime();
  time = - time;
  printf("Thread %d dosky calc %.2f ms\n", omp_get_thread_num(), time*1000);
  
  sprintf(f_name, f_fmt_out, tconfig.seq, ndx-st+1);
  write_single_fits(f_name, proc);

  free(data->data);
  free(proc->data);
  free(data);
  free(proc);

  return(0);
}

int skysub(config tconfig, fitsobj *sky, int ndx, int st){  

  fitsobj *data = malloc(sizeof(fitsobj));
  fitsobj *proc = malloc(sizeof(fitsobj));
  char f_name[100];
  int npixels, nimages;
  int pix_img;
  int i, j;
  char f_fmt_in[100] = "darksub/IPA_seq%d_dith_%d.fits";
  char f_fmt_out[100] = "skysub/IPA_seq%d_dith_%d.fits";
  double time;

  sprintf(f_name, f_fmt_in, tconfig.seq, ndx-st+1);  
  read_fits(f_name, data);

  proc->naxis1 = data->naxis1;
  proc->naxis2 = data->naxis2;
  proc->nramps = data->nramps;
  proc->ngroups = data->ngroups;

  npixels = proc->naxis1*proc->naxis2*(proc->ngroups-1)*proc->nramps;
  proc->data = malloc(sizeof(float)*npixels);

  
  pix_img = proc->naxis1*proc->naxis2;
  nimages = (proc->ngroups-1)*proc->nramps;
  
  time = omp_get_wtime();

  #pragma omp parallel private(j) num_threads(1)
  {
    #pragma omp for 
    for (i=0; i<pix_img; i++){
      for (j=0; j<nimages; j++)
	proc->data[i+j*pix_img] = data->data[i+j*pix_img]-sky->data[i];
    }
  }
  
  printf("Thread %d skysub calc in %.2f ms\n", 
	 omp_get_thread_num(), (omp_get_wtime()-time)*1000);

  sprintf(f_name, f_fmt_out, tconfig.seq, ndx-st+1);  
  write_fits(f_name, proc);


  free(data->data);
  free(proc->data);
  free(data);
  free(proc);

  return(0);
}

int flatdivide(config tconfig, fitsobj *flat, int ndx, int st){  

  fitsobj *data = malloc(sizeof(fitsobj));
  fitsobj *proc = malloc(sizeof(fitsobj));
  char f_name[100];
  int npixels, nimages;
  int pix_img;
  int i, j;
  char f_fmt_in[100] = "fftcorr/IPA_seq%d_dith_%d.fits";
  char f_fmt_out[100] = "flat_divide/IPA_seq%d_dith_%d.fits";
  double time;

  sprintf(f_name, f_fmt_in, tconfig.seq, ndx-st+1);  
  read_fits(f_name, data);

  proc->naxis1 = data->naxis1;
  proc->naxis2 = data->naxis2;
  proc->nramps = data->nramps;
  proc->ngroups = data->ngroups;

  npixels = proc->naxis1*proc->naxis2*(proc->ngroups-1)*proc->nramps;
  proc->data = malloc(sizeof(float)*npixels);
  
  pix_img = proc->naxis1*proc->naxis2;
  nimages = (proc->ngroups-1)*proc->nramps;
  
  time = omp_get_wtime();

  for (j=0; j<nimages; j++)
    for (i=0; i<pix_img; i++)
      proc->data[i+j*pix_img] = data->data[i+j*pix_img]/
	flat->data[tconfig.band_ndx[0]*2048+i];

  printf("Thread %d flat_divide calc in %.2f ms\n", 
	 omp_get_thread_num(), (omp_get_wtime()-time)*1000);

  sprintf(f_name, f_fmt_out, tconfig.seq, ndx-st+1);  
  write_fits(f_name, proc);

  free(data->data);
  free(proc->data);
  free(data);
  free(proc);

  return(0);
}

int badpixrem(config tconfig, fitsobj *badpix, int ndx, int st){  

  fitsobj *data = malloc(sizeof(fitsobj));
  fitsobj *proc = malloc(sizeof(fitsobj));
  char f_name[100];
  int npixels, nimages;
  int pix_img;
  int i;
  char f_fmt_in[100] = "flat_divide/IPA_seq%d_dith_%d.fits";
  char f_fmt_out[100] = "badpixrem/IPA_seq%d_dith_%d.fits";
  double time;

  sprintf(f_name, f_fmt_in, tconfig.seq, ndx-st+1);  
  read_fits(f_name, data);

  proc->naxis1 = data->naxis1;
  proc->naxis2 = data->naxis2;
  proc->nramps = data->nramps;
  proc->ngroups = data->ngroups;

  npixels = proc->naxis1*proc->naxis2*(proc->ngroups-1)*proc->nramps;
  proc->data = malloc(sizeof(float)*npixels);
  
  pix_img = proc->naxis1*proc->naxis2;
  nimages = (proc->ngroups-1)*proc->nramps;
  
  time = omp_get_wtime();

#pragma omp parallel num_threads(1)
  {
    //float *badpixim = malloc(sizeof(float)*pix_img);
    //float *image = malloc(sizeof(float)*pix_img);
    //float *data_ref = &data->data[0];
    float *proc_ref = &proc->data[0];
    float *image = &data->data[0];
    float *badpixim = &badpix->data[0];

    #pragma omp for
    for (i=0; i<nimages; i++){
      //memcpy(image, &data_ref[i*pix_img], sizeof(float)*pix_img);
      //memcpy(badpixim, &badpix->data[0], sizeof(float)*pix_img);
      badpixfun(&image[i*pix_img], &badpixim[0], data->naxis1, data->naxis2);
      memcpy(&proc_ref[i*pix_img], &image[i*pix_img], sizeof(float)*pix_img);
    }

    //free(badpixim);
    //free(image);
  }
  
  printf("Thread %d badpix calc in %.2f ms\n", 
	 omp_get_thread_num(), (omp_get_wtime()-time)*1000);

  sprintf(f_name, f_fmt_out, tconfig.seq, ndx-st+1);  
  write_fits(f_name, proc);

  free(data->data);
  free(proc->data);
  free(data);
  free(proc);

  return(0);
}

int fftcorr(config tconfig, int ndx, int st){  

  fitsobj *data = malloc(sizeof(fitsobj));
  fitsobj *proc = malloc(sizeof(fitsobj));
  char f_name[100];
  int npixels, nimages;
  int pix_img;
  int i;
  char f_fmt_in[100] = "skysub/IPA_seq%d_dith_%d.fits";
  char f_fmt_out[100] = "fftcorr/IPA_seq%d_dith_%d.fits";
  double time;

  sprintf(f_name, f_fmt_in, tconfig.seq, ndx-st+1);  
  read_fits(f_name, data);

  proc->naxis1 = data->naxis1;
  proc->naxis2 = data->naxis2;
  proc->nramps = data->nramps;
  proc->ngroups = data->ngroups;

  npixels = proc->naxis1*proc->naxis2*(proc->ngroups-1)*proc->nramps;
  proc->data = malloc(sizeof(float)*npixels);
  
  pix_img = proc->naxis1*proc->naxis2;
  nimages = (proc->ngroups-1)*proc->nramps;
  
  time = omp_get_wtime();

  #pragma omp parallel  num_threads(1)
  {
    int npix_32 = data->naxis1*data->naxis2/32;
    float *image = malloc(sizeof(float)*pix_img);
    fftw_complex *ff2dfor = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*npix_32);
    fftw_complex *ff2dback = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*npix_32);
    fftw_plan pfor, pback;
    #pragma omp critical
    {
    pfor = fftw_plan_dft_2d(data->naxis2, 64, ff2dfor, ff2dfor, FFTW_FORWARD, FFTW_ESTIMATE);
    pback = fftw_plan_dft_2d(data->naxis2, 64, ff2dfor, ff2dback, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_free(ff2dfor);
    fftw_free(ff2dback);
    }

    #pragma omp for
    for (i=0; i<nimages; i++){
      memcpy(image, &data->data[i*pix_img], sizeof(float)*pix_img);
      fftcorrfun(image, data->naxis1, data->naxis2, pfor, pback);
      memcpy(&proc->data[i*pix_img], &image[0], sizeof(float)*pix_img);
      //fftcorrfun(&proc->data[i*pix_img], data->naxis1, data->naxis2, pfor, pback);
    }

    #pragma omp critical
    {
    free(image);
    fftw_destroy_plan(pfor);
    fftw_destroy_plan(pback);
    }
  }
  
  printf("Thread %d fftcorr calc in %.2f ms\n", 
	 omp_get_thread_num(), (omp_get_wtime()-time)*1000);

  sprintf(f_name, f_fmt_out, tconfig.seq, ndx-st+1);  
  write_fits(f_name, proc);//

  free(data->data);
  free(proc->data);
  free(data);
  free(proc);

  return(0);
}

int tile_images(config tconfig, int ndx, int st){  

  fitsobj *data = malloc(sizeof(fitsobj));
  fitsobj *proc = malloc(sizeof(fitsobj));
  char f_name[100];
  int nimages;
  int pix_img;
  int i, j, k;
  char f_fmt_in[100] = "badpixrem/IPA_seq%d_dith_%d.fits";
  char f_fmt_out[100] = "tile/Image_%04d.fits";
  double time;
  //int bin = 1;
  int crop[2] = {1, 30};

  sprintf(f_name, f_fmt_in, tconfig.seq, ndx-st+1);  
  read_fits(f_name, data);
  
  pix_img = (data->naxis1)*data->naxis2;
  nimages = (data->ngroups-1)*data->nramps;

  proc->data = malloc(sizeof(float)*pix_img);
  proc->naxis1 = (crop[1]-crop[0])*64;
  proc->naxis2 = data->naxis2;
  
  time = omp_get_wtime();

  #pragma omp parallel num_threads(1)
  {
    #pragma omp for 
    for (i=0; i<nimages; i++){
      sprintf(f_name, f_fmt_out, (tconfig.seq-1)*tconfig.n_dith*nimages+(ndx-st)*nimages+i+1);  
      for (j=0; j<proc->naxis2;j++) 
	for (k=crop[0]*64; k<proc->naxis1+crop[0]*64;k++) 
	  proc->data[j*proc->naxis1+k-crop[0]*64] = data->data[i*pix_img+j*data->naxis1+k]; 
      for (j=0; j<proc->naxis2;j++) 
	for (k=64*5; k<64*6;k++) 
	  proc->data[j*proc->naxis1+k-crop[0]*64] = 0; 
      write_single_fits(f_name, proc);
    }
  }
  
  printf("Thread %d tile calc in %.2f ms\n", 
	 omp_get_thread_num(), (omp_get_wtime()-time)*1000);


  free(proc->data);
  free(data->data);
  free(data);
  free(proc);

  return(0);
}


int darksub_all(config tconfig){

  int i;   double time;
  printf("%s\n", tconfig.f_dark);

#pragma omp parallel num_threads(NUM_THREADS)
  {

    fitsobj *dark = malloc(sizeof(fitsobj));
    read_fits(tconfig.f_dark, dark);

    #pragma omp master
    {
      printf("**Darksub with %d threads\n", omp_get_num_threads());
      time = omp_get_wtime();
    }

    #pragma omp for
    for (i=tconfig.st; i<tconfig.st+tconfig.n_dith; i++){
      darksub(tconfig, dark, i, tconfig.st);
    }

    #pragma omp barrier

    #pragma omp master
    printf("Darksub %.2f ms\n", (omp_get_wtime()-time)*1000);

    free(dark->data);
    free(dark);

    }

  return(0);
}

int doskies_all(config tconfig){

  int i; double time;
  
#pragma omp parallel private(time) num_threads(NUM_THREADS)
  {
#pragma omp master
    printf("**Dosky with %d threads\n", omp_get_num_threads());
    time = omp_get_wtime();

#pragma omp for
    for (i=tconfig.st; i<tconfig.st+tconfig.n_dith; i++){
      doskies(tconfig, i, tconfig.st);
    }

    #pragma omp master
    printf("dosky %.2f ms\n", (omp_get_wtime()-time)*1000);
  }

  
  return(0);
}

int skysub_all(config tconfig){

  int i,j; double time;
  float *tmp_skies;
  double tmp_arr[tconfig.n_dith];
  char f_fmt_in[100] = "skies/IPA_seq%d_dith_%d.fits";
  fitsobj * sky = malloc(sizeof(fitsobj));
  int npixels;

  char f_name[100];

#pragma omp parallel private(tmp_arr) private(f_name) num_threads(NUM_THREADS)
  {
    fitsobj *skies = malloc(sizeof(fitsobj));

    //id = omp_get_thread_num();

    #pragma omp master
    printf("**Skysub with %d threads\n", omp_get_num_threads());
    time = omp_get_wtime();

    #pragma omp single
    {
      sprintf(f_name, f_fmt_in, tconfig.seq, 1);
      read_single_fits(f_name, sky);
      npixels = sky->naxis1*sky->naxis2;
      tmp_skies = malloc(sizeof(float)*npixels*tconfig.n_dith);
      sky->data = malloc(sizeof(float)*npixels); 
    }

    #pragma omp for private(j)
    for (i=tconfig.st; i<tconfig.st+tconfig.n_dith; i++){
      sprintf(f_name, f_fmt_in, tconfig.seq, i-tconfig.st+1);
      read_single_fits(f_name, skies);
      for (j=0; j<npixels; j++)
	tmp_skies[j+(i-tconfig.st)*npixels] = skies->data[j];
    }


    #pragma omp barrier
    free(skies->data);
    free(skies);


    #pragma omp for private(j)
    for (i=0; i<npixels; i++){
      for (j=0; j<tconfig.n_dith; j++)
	tmp_arr[j] = (double) tmp_skies[i+j*npixels];
      gsl_sort(tmp_arr, 1, tconfig.n_dith);
      sky->data[i] = (float) gsl_stats_median_from_sorted_data(tmp_arr, 1, tconfig.n_dith);
    }

    #pragma omp for
    for (i=tconfig.st; i<tconfig.st+tconfig.n_dith; i++){
      skysub(tconfig, sky, i, tconfig.st);
    }

    
    #pragma omp barrier

    #pragma omp master
    printf("Skysub %.2f ms\n", (omp_get_wtime()-time)*1000);
  }

  fitsobj *flat = malloc(sizeof(fitsobj));
  read_single_fits(tconfig.f_flat, flat);
  for (i=0; i<npixels;i++)
    sky->data[i] /= flat->data[tconfig.band_ndx[0]*2048+i];
  sprintf(f_name, "skies/IPA_deq%d_chk.fits", tconfig.seq);
  write_single_fits(f_name, sky);  
  
  free(flat->data);
  free(sky->data);
  free(sky);
  free(flat);
  return(0);
}

int flatdivide_all(config tconfig){

  int i; double time;
  
#pragma omp parallel num_threads(NUM_THREADS)
  {
#pragma omp master
    printf("**flatdivide with %d threads\n", omp_get_num_threads());
    time = omp_get_wtime();

  fitsobj *flat = malloc(sizeof(fitsobj));
  read_single_fits(tconfig.f_flat, flat);

#pragma omp for
    for (i=tconfig.st; i<tconfig.st+tconfig.n_dith; i++){
      flatdivide(tconfig, flat, i, tconfig.st);
    }

    free(flat->data); 
    free(flat);
 
#pragma omp master
    printf("flatdivide %.2f ms\n", (omp_get_wtime()-time)*1000);
  }

  
  return(0);
}

int badpixrem_all(config tconfig){

  int i; double time;
  
#pragma omp parallel num_threads(NUM_THREADS)
  {
    fitsobj *badpix = malloc(sizeof(fitsobj));
    read_single_fits(tconfig.f_badp, badpix);
    #pragma omp single
    {
      printf("**badpix with %d threads\n", omp_get_num_threads());
      time = omp_get_wtime();
    }


  #pragma omp for
    for (i=tconfig.st; i<tconfig.st+tconfig.n_dith; i++){
      badpixrem(tconfig, badpix, i, tconfig.st);
    }

  #pragma omp barrier
  
  #pragma omp master
    printf("badpix %.2f ms\n", (omp_get_wtime()-time)*1000);

    free(badpix->data);
    free(badpix);

  }

  
  return(0);
}

int fftcorr_all(config tconfig){

  int i; double time;
  
#pragma omp parallel private(time) num_threads(NUM_THREADS)
  {
    #pragma omp master
    printf("**fftcorr with %d threads\n", omp_get_num_threads());
    time = omp_get_wtime();

  #pragma omp for
  for (i=tconfig.st; i<tconfig.st+tconfig.n_dith; i++){
    fftcorr(tconfig, i, tconfig.st);
  }

    #pragma omp master
    printf("fftcorr %.2f ms\n", (omp_get_wtime()-time)*1000);
  }

  
  return(0);
}

int tile_all(config tconfig){

  int i; double time;
  
#pragma omp parallel private(time) num_threads(NUM_THREADS)
  {
    #pragma omp master
    printf("**tile with %d threads\n", omp_get_num_threads());
    time = omp_get_wtime();

  #pragma omp for
  for (i=tconfig.st; i<tconfig.st+tconfig.n_dith; i++){
    tile_images(tconfig, i, tconfig.st);
  }

    #pragma omp master
    printf("tile %.2f ms\n", (omp_get_wtime()-time)*1000);
  }

  
  return(0);
}


int badpixfun(float *image, float *badpixim, int naxis1, int naxis2){

  int i, j, k, l;
  double tmp_list[9];
  double med, med2, std, std2;
  float* im_proc = malloc(sizeof(float)*naxis1*naxis2);
  double* tmp_im = malloc(sizeof(double)*naxis1*naxis2);
  double tmp_arr[2048];

  long npix = naxis1*naxis2;

  for (j=0; j<naxis2; j++){
    for (i=0; i<naxis1; i++){
      if(!isfinite(image[naxis1*j+i]))
	image[naxis1*j+i] = 0;
    }
  }

  for (i=0; i<npix; i++) tmp_im[i] = image[i];
  gsl_sort(tmp_im, 1, npix);
  med2 = gsl_stats_median_from_sorted_data(&tmp_im[npix/20], 1, 18*npix/20);
  i = 3;
  while(1){
    std2 = gsl_stats_sd_m(&tmp_im[i*npix/20], 1, (20-2*i)*npix/20, med2);
    i++;
    if (isfinite(std2)) break;
    if (i > 8) break;
  }
  
  free(tmp_im);

  for (j=1; j<naxis2-1; j++){
    for (i=1; i<naxis1-1; i++){
      if ((badpixim[naxis1*j+i] == 1) | ((image[naxis1*j+i]-med2) < -3*std2)){
	for (l=j-1; l<j+2; l++){
	  for (k=i-1; k<i+2; k++){
	    tmp_list[(l-j+1)*3+k-i+1] = (double) image[naxis1*l+k];
	  } 
	}
	tmp_list[4] = tmp_list[8];
      
	gsl_sort(tmp_list, 1, 8);
	med = (tmp_list[3]+tmp_list[4])*.5;
	image[naxis1*j+i] = (float) med;
      }
    }
  }

  
  for (j=1; j<naxis2-1; j++){
    for (i=1; i<naxis1-1; i++){
      for (l=j-1; l<j+2; l++){
	for (k=i-1; k<i+2; k++){
	  tmp_list[(l-j+1)*3+k-i+1] = (double) image[naxis1*l+k];
	} 
      }
      tmp_list[4] = tmp_list[8];
      
      gsl_sort(tmp_list, 1, 8);
      med = (tmp_list[3]+tmp_list[4])*.5;
      std = gsl_stats_sd_m(&tmp_list[1], 1, 6, 1);

      if (abs(image[naxis1*j+i]-med) > 5*std)
	im_proc[naxis1*j+i] = (float) med;
      else
	im_proc[naxis1*j+i] = image[naxis1*j+i];
      
    }
  }

  memcpy(image, im_proc, sizeof(float)*naxis1*naxis2);   

  for (j = 0; j<naxis2; j++ ){
    for(i=0; i<naxis1; i++)
      tmp_arr[i] = image[naxis1*j+i];
    
    gsl_sort(tmp_arr, 1, 2048);
    med = gsl_stats_median_from_sorted_data(tmp_arr, 1, 2048);
    for(i=0; i<naxis1; i++)
      image[naxis1*j+i] -= med;
    for(i=1472; i<1472+64; i++)
      image[naxis1*j+i] = 0;
    for(i=320; i<320+64; i++)
      image[naxis1*j+i] = 0;

  }

  free(im_proc);
  return(0);

}


int fftcorrfun(float* image, int naxis1, int naxis2, fftw_plan pfor, fftw_plan pback){

  int npix = naxis1*naxis2;
  int npix_32 = npix /32;
  fftw_complex *ff2dfor = 
    (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*npix_32);
  fftw_complex *ff2dback = 
    (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*npix_32);
  float *tmp_im = malloc(sizeof(float)*npix);
  float *orig = malloc(sizeof(float)*npix);
  double *tmp_im2 = malloc(sizeof(double)*npix);
  int i, j, k;
  int stride;
  double tmp_arr[2048];
  double mean, std;
  double med;

  memcpy(orig, image, sizeof(float)*npix);
  
  
  for (k=0; k<32;k++){
    stride = k*64;
    for (j=0; j<naxis2; j++ ){
      for (i=stride; i<64+stride; i++){
	ff2dfor[64*j+i-stride] = (double) image[naxis1*j+i]+0*I;
      }
    }

    fftw_execute_dft(pfor, ff2dfor, ff2dfor);

    for (j=0; j<naxis2; j++ ){
      for (i=stride; i<64+stride; i++){
	if (i != stride)
	  tmp_im[naxis1*j+i] = (float) (cabs(ff2dfor[64*j+i-stride]));
	else
	  tmp_im[naxis1*j+i] = 0;
      }
    }
        
    for (j=1; j<naxis2-1; j++ ){
      for (i=stride; i<64+stride; i++){
	image[naxis1*j+i] = (tmp_im[naxis1*(j+1)+i] + tmp_im[naxis1*(j-1)+i]
			   -2*tmp_im[naxis1*j+i]);
	tmp_im2[naxis1*j+i] = image[naxis1*j+i];
      }
    }

    mean =0; std=0;

    for(i=stride; i<stride+64; i++){
      for (j=0; j< naxis2; j++){
    	  ff2dback[64*j+i-stride] = ff2dfor[64*j+i-stride];
      }
    }


    for(i=stride+2; i<stride+62; i++){
      mean = gsl_stats_mean(&tmp_im2[naxis1+i], naxis1, naxis2-2); 
      std = gsl_stats_sd(&tmp_im2[naxis1+i], naxis1, naxis2-2); 
      for (j=2; j< naxis2-2; j++){
	image[naxis1*j+i] -= mean;
	image[naxis1*j+i] /= std;     
	if (abs(image[naxis1*j+i]) > 3){
	  ff2dback[64*j+i-stride] = mean+0.*0.5*ff2dfor[64*j+i-stride]/cabs(ff2dfor[64*j+i-stride])* 
	    (cabs(ff2dfor[64*(j+1)+i-stride])+cabs(ff2dfor[64*(j-1)+i-stride]));
	}	
      }
    }

    fftw_execute_dft(pback, ff2dback, ff2dback);

    for (j=0; j<naxis2; j++ ){
      for (i=stride; i<64+stride; i++){
	image[naxis1*j+i] = (float) creal(ff2dback[64*j+i-stride])/npix_32 ;
      }
    }
  }


  //for (i=0; i< npix; i++) image[i] = -image[i]+orig[i];

  free(orig);
  free(tmp_im);
  free(tmp_im2);			  
  fftw_free(ff2dfor);
  fftw_free(ff2dback);

  return(0);
}
