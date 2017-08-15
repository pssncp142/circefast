#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "fitsio.h"
#include "omp.h"

#include "fileio.h"

#define VERBOSE 1

//print fits object details.
int fits_info(fitsobj* fitso, char* func){

  char error[30];
  printf("Nramps      : %4d, Ngroups     : %4d\n", 
	 fitso->nramps, fitso->ngroups);
  printf("Naxis1      : %4d, Naxis2      : %4d\n", 
	 fitso->naxis1, fitso->naxis2);
  fits_get_errstatus(fitso->status, error);

  if (fitso->status > 0){
    printf("Exiting with status : %4d %s\n", fitso->status, error);
    exit(0);
  }

  return(0);
}

//initialize fits object structure
int init_fits_object(fitsfile* fptr, fitsobj* obj, char* f_name, int opt, 
		     char* func){

  int hdutype;
  char card[FLEN_CARD];
  int id = omp_get_thread_num();

  #pragma omp critical
  {
  if(VERBOSE) printf("%s[#%d] => Init fits object for : %s\n", 
		     func, id, f_name);

  strcpy(obj->f_name, f_name);
  
  if (opt == 0){
  //get ngroups and nramps information from PrimaryHDU
  fits_movabs_hdu(fptr, 1, &hdutype, &obj->status);
  fits_read_key(fptr, TINT, "NGROUPS", &obj->ngroups, card, &obj->status);
  fits_read_key(fptr, TINT, "NRAMPS", &obj->nramps, card, &obj->status);
  } else if (opt == 1) { // to read single fits
    obj->ngroups = 2;
    obj->nramps = 1;
  } else if (opt == 2) { // to read single fits
    obj->ngroups = 2;
    fits_movabs_hdu(fptr, 1, &hdutype, &obj->status);
    fits_read_key(fptr, TINT, "NRAMPS", &obj->nramps, card, &obj->status);
  }


  if ((opt == 0) || (opt == 2))
    fits_movabs_hdu(fptr, 2, &hdutype, &obj->status);    
  else
    fits_movabs_hdu(fptr, 1, &hdutype, &obj->status);
  
  //get the size of the image
  fits_read_key(fptr, TINT, "NAXIS1", &obj->naxis1, card, &obj->status);
  fits_read_key(fptr, TINT, "NAXIS2", &obj->naxis2, card, &obj->status);

  //Accordingly set the data size
  obj->data = malloc(sizeof(float)*obj->naxis1*obj->naxis2*obj->nramps*
		     (obj->ngroups-1));  

  fits_info(obj, &func[0]);

  }
  return(1);
}

//read fits into a fits object
int read_fits(char *f_name, fitsobj* obj){

  char func[100] = "read_fits";
  int hdutype, npixels, anynull;
  int nullval = 0;
  fitsfile *fptr;
  int i,j;
  int hdu_num;

  obj->status = 0;
  fits_open_file(&fptr, f_name, READONLY, &obj->status); // get the pointer
  init_fits_object(fptr, obj, f_name, 0, func);

  npixels = obj->naxis1*obj->naxis2;

  for (i=0; i<obj->nramps; i++){
    for (j=0; j<obj->ngroups-1; j++){ 
      hdu_num = 2+j+i*(obj->ngroups-1);
      fits_movabs_hdu(fptr, hdu_num, &hdutype, &obj->status);
      fits_read_img(fptr, TFLOAT, 1, npixels, &nullval, 
		    &obj->data[npixels*(i*(obj->ngroups-1)+j)],
		    &anynull, &obj->status);
    }
  }
  
  fits_close_file(fptr, &obj->status);
  return(0);
}

//read fits into a fits object. 
//Need to fix nramps, nroups
int read_single_fits(char *f_name, fitsobj* obj){

  char func[100] = "read_single_fits";
  int npixels, anynull;
  int nullval = 0;
  fitsfile *fptr;

  obj->status = 0;
  fits_open_file(&fptr, f_name, READONLY, &obj->status); // get the pointer
  init_fits_object(fptr, obj, f_name, 1, func);

  npixels = obj->naxis1*obj->naxis2;

  fits_read_img(fptr, TFLOAT, 1, npixels, &nullval, 
		&obj->data[0], &anynull, &obj->status);

  fits_close_file(fptr, &obj->status);
  return(0);
}

//read raw files in data2, data. Resulting images data2-data.
int read_diff_fits(char *f_name, fitsobj* obj){

  char func[100] = "read_diff_fits";
  int hdutype, npixels, anynull;
  int nullval = 0;
  fitsfile *fptr;
  int i,j,k;
  int hdu_num;
  int st_pos;
  
  float* tmp_img1, *tmp_img2, *data_ptr;

  obj->status = 0;
  fits_open_file(&fptr, f_name, READONLY, &obj->status); // get the pointer
  init_fits_object(fptr, obj, f_name, 0, func);

  npixels = obj->naxis1*obj->naxis2;

  tmp_img1 = malloc(sizeof(float)*obj->naxis1*obj->naxis2);
  tmp_img2 = malloc(sizeof(float)*obj->naxis1*obj->naxis2);
  data_ptr = &obj->data[0];

  for (i=0; i<obj->nramps; i++){
    for (j=0; j<obj->ngroups-1; j++){
      hdu_num = 2+j+i*(obj->ngroups);
      fits_movabs_hdu(fptr, hdu_num, &hdutype, &obj->status);
      fits_read_img(fptr, TFLOAT, 1, npixels, &nullval, 
		    &tmp_img1[0],
		    &anynull, &obj->status);
      fits_movabs_hdu(fptr, hdu_num+1, &hdutype, &obj->status);
      fits_read_img(fptr, TFLOAT, 1, npixels, &nullval, 
		    &tmp_img2[0],
		    &anynull, &obj->status);
      st_pos = ((obj->ngroups-1)*i+j)*npixels;
      for (k=0; k<npixels; k++) 
	data_ptr[st_pos+k] = tmp_img2[k] - tmp_img1[k];
 
    }
  }

  free(tmp_img1);
  free(tmp_img2);
  fits_close_file(fptr, &obj->status);
  return(0);
}

//read raw files in data2, data. Resulting images data2-data.
// Attempt to read final - first diff image
int read_diff_fits_framp(char *f_name, fitsobj* obj){

  char func[100] = "read_diff_fits";
  int hdutype, npixels, anynull;
  char card[FLEN_CARD];
  int nullval = 0;
  fitsfile *fptr;
  int i,j,k;
  int hdu_num;
  int st_pos;
  
  float* tmp_img1, *tmp_img2, *data_ptr;

  obj->status = 0;
  fits_open_file(&fptr, f_name, READONLY, &obj->status); // get the pointer
  init_fits_object(fptr, obj, f_name, 2, func);

  npixels = obj->naxis1*obj->naxis2;

  tmp_img1 = malloc(sizeof(float)*obj->naxis1*obj->naxis2);
  tmp_img2 = malloc(sizeof(float)*obj->naxis1*obj->naxis2);
  data_ptr = &obj->data[0];

  fits_movabs_hdu(fptr, 1, &hdutype, &obj->status);
  fits_read_key(fptr, TINT, "NGROUPS", &obj->ngroups, card, &obj->status);

  for (i=0; i<obj->nramps; i++){
      hdu_num = 2+i*(obj->ngroups);
      fits_movabs_hdu(fptr, hdu_num, &hdutype, &obj->status);
      fits_read_img(fptr, TFLOAT, 1, npixels, &nullval, 
		    &tmp_img1[0],
		    &anynull, &obj->status);
      fits_movabs_hdu(fptr, hdu_num+obj->ngroups-1, &hdutype, &obj->status);
      fits_read_img(fptr, TFLOAT, 1, npixels, &nullval, 
		    &tmp_img2[0],
		    &anynull, &obj->status);
      st_pos = i*npixels;
      for (k=0; k<npixels; k++) 
	data_ptr[st_pos+k] = tmp_img2[k] - tmp_img1[k];
 
  }

  obj->ngroups = 2;

  free(tmp_img1);
  free(tmp_img2);
  fits_close_file(fptr, &obj->status);
  return(0);
}

int write_single_fits(char *f_name, fitsobj *obj){

  fitsfile *fptr;
  int status;
  long naxes[2] = {obj->naxis1, obj->naxis2};
  long npixels = naxes[0]*naxes[1];

  status= 0;
  remove(f_name);

  fits_create_file(&fptr, f_name, &status);
  fits_create_img(fptr, FLOAT_IMG, 2, &naxes[0],  &status);
  fits_write_img(fptr, TFLOAT, 1, npixels, 
		 &obj->data[0], &status);
  
  fits_close_file(fptr, &status);
  return(0);
}


int write_fits(char *f_name, fitsobj *obj){

  fitsfile *fptr;
  int status;
  long naxes[2] = {obj->naxis1, obj->naxis2};
  long npixels = naxes[0]*naxes[1];
  int i, j;

  status= 0;
  remove(f_name);

  //printf("%s\n", f_name);

  fits_create_file(&fptr, f_name, &status);
  fits_create_img(fptr, FLOAT_IMG, 0, &naxes[0],  &status);
  fits_update_key(fptr, TINT, "NRAMPS", &obj->nramps, " ", &status);
  fits_update_key(fptr, TINT, "NGROUPS", &obj->ngroups, " ", &status);

  for (i=0; i<obj->nramps; i++){
    for (j=0; j<obj->ngroups-1; j++){
      fits_create_img(fptr, FLOAT_IMG, 2, &naxes[0],  &status);
      //printf("stat:%d\n", status);
      fits_write_img(fptr, TFLOAT, 1, npixels, 
		     &obj->data[npixels*(i*(obj->ngroups-1)+j)], &status);
      //printf("stat2:%d\n", status);
    }
  }
  
  fits_close_file(fptr, &status);
  return(0);
}
