#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "fitsio.h"
#include "omp.h"

#include "fileio.h"

int fits_info(fitsobj* fitso, char* func){

  #pragma omp critical
  {
  int id = omp_get_thread_num();
  printf("-> Thread %d, I/O func %s\n", id, func);
  printf("Filename    : %-30s\n", fitso->f_name);
  printf("Nramps      : %4d, Ngroups     : %4d\n", 
	 fitso->nramps, fitso->ngroups);
  printf("Naxis1      : %4d, Naxis2      : %4d\n", 
	 fitso->naxis1, fitso->naxis2);
  printf("Status      : %4d\n", fitso->status);
  }

  return(0);
}


//read fits into a fits object
int read_fits(char *f_name, fitsobj* obj){

  char func[100] = "read_fits";
  int status, hdutype, npixels, anynull;
  int nullval = 0;
  fitsfile *fptr;
  char card[FLEN_CARD];
  int i,j;
  int hdu_num;

  obj->status = 0;

  fits_open_file(&fptr, f_name, READONLY, &obj->status);
  
  fits_movabs_hdu(fptr, 1, &hdutype, &obj->status);
  fits_read_key(fptr, TINT, "NGROUPS", &obj->ngroups, card, &obj->status);
  fits_read_key(fptr, TINT, "NRAMPS", &obj->nramps, card, &obj->status);

  fits_movabs_hdu(fptr, 2, &hdutype, &obj->status);
  fits_read_key(fptr, TINT, "NAXIS1", &obj->naxis1, card, &obj->status);
  fits_read_key(fptr, TINT, "NAXIS2", &obj->naxis2, card, &obj->status);

  obj->data = malloc(sizeof(float)*obj->naxis1*obj->naxis2*obj->nramps*
		     (obj->ngroups-1));
  npixels = obj->naxis1*obj->naxis2;

  //obj->status = status;
  strcpy(obj->f_name, f_name);
  fits_info(obj, &func[0]);
  if (obj->status > 0) exit(1);

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
  int status, hdutype, npixels, anynull;
  int nullval = 0;
  fitsfile *fptr;
  char card[FLEN_CARD];

  status = 0;

  fits_open_file(&fptr, f_name, READONLY, &status);
  
  fits_movabs_hdu(fptr, 1, &hdutype, &status);

  fits_read_key(fptr, TINT, "NAXIS1", &obj->naxis1, card, &status);
  fits_read_key(fptr, TINT, "NAXIS2", &obj->naxis2, card, &status);

  obj->data = malloc(sizeof(float)*obj->naxis1*obj->naxis2);
  npixels = obj->naxis1*obj->naxis2;

  fits_read_img(fptr, TFLOAT, 1, npixels, &nullval, 
		&obj->data[0], &anynull, &status);

  obj->status = status;
  strcpy(obj->f_name, f_name);
  fits_info(obj, &func[0]);
  if (obj->status > 0) exit(1);

  fits_close_file(fptr, &status);
  return(0);
}

//read raw files in data2, data. Resulting images data2-data.
int read_diff_fits(char *f_name, fitsobj* obj){

  char func[100] = "read_diff_fits";
  int status, hdutype, npixels, anynull;
  int nullval = 0;
  fitsfile *fptr;
  char card[FLEN_CARD];
  int i,j;
  int hdu_num;

  status = 0;

  fits_open_file(&fptr, f_name, READONLY, &status);
  
  fits_movabs_hdu(fptr, 1, &hdutype, &status);
  fits_read_key(fptr, TINT, "NGROUPS", &obj->ngroups, card, &status);
  fits_read_key(fptr, TINT, "NRAMPS", &obj->nramps, card, &status);
  
  fits_movabs_hdu(fptr, 2, &hdutype, &status);

  fits_read_key(fptr, TINT, "NAXIS1", &obj->naxis1, card, &status);
  fits_read_key(fptr, TINT, "NAXIS2", &obj->naxis2, card, &status);

  obj->data = malloc(sizeof(float)*obj->naxis1*obj->naxis2*obj->nramps*
		     (obj->ngroups-1));
  obj->data2 = malloc(sizeof(float)*obj->naxis1*obj->naxis2*obj->nramps*
		     (obj->ngroups-1));

  npixels = obj->naxis1*obj->naxis2;

  obj->status = status;
  strcpy(obj->f_name, f_name);
  fits_info(obj, &func[0]);
  if (obj->status > 0) exit(1);

  for (i=0; i<obj->nramps; i++){
    for (j=0; j<obj->ngroups-1; j++){
      hdu_num = 2+j+i*(obj->ngroups);
      fits_movabs_hdu(fptr, hdu_num, &hdutype, &status);
      fits_read_img(fptr, TFLOAT, 1, npixels, &nullval, 
		    &obj->data[npixels*(i*(obj->ngroups-1)+j)],
		    &anynull, &status);
      fits_movabs_hdu(fptr, hdu_num+1, &hdutype, &status);
      fits_read_img(fptr, TFLOAT, 1, npixels, &nullval, 
		    &obj->data2[npixels*(i*(obj->ngroups-1)+j)],
		    &anynull, &status);
 
    }
  }

  //free(tmp);
  fits_close_file(fptr, &status);
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
