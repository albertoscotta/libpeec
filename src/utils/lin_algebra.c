/**
 * \file
 * \author Scotta Alberto
 *
\verbatim
 |¯¯| |    |¯¯\ \  /
 |--| |    |--<  \/
 |  | |___ |__/  /   
\endverbatim
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "lin_algebra.h"

SPARSE_MAT* create_sparse_mat(int row_num, int col_num) {
  SPARSE_MAT* mat;

  mat = malloc(sizeof(SPARSE_MAT));
  if(!mat) {
    #ifdef DEBUG
    fprintf(stderr, "Error: create_sparse_mat(), 'mat', malloc() failed.\n");
    #endif 
    return NULL;
  }
  mat->row_num = row_num;
  mat->col_num = col_num;
  mat->nz_max = 1;
  mat->nz = 0;
  mat->i = malloc(sizeof(int));
  mat->j = malloc(sizeof(int));
  mat->value = malloc(sizeof(double));
  if(!mat->i || !mat->j || !mat->value) {
    #ifdef DEBUG
    fprintf(stderr, "Error: create_spare_mat(), malloc() failed.\n");
    #endif
    destroy_sparse_mat(mat);
    return NULL;
  }

  return mat;
}

void destroy_sparse_mat(SPARSE_MAT* mat) {

  if(!mat)
    return;

  if(!mat->i)
    free(mat->i);
  if(!mat->j)
    free(mat->j);
  if(!mat->value)
    free(mat->value);
  free(mat);

  return;
}


int sparse_mat_add(SPARSE_MAT* mat, int row, int col, double value) {

  int* mat_i, * mat_j;
  double* mat_value;

  if(!mat) {
    #ifdef DEBUG
    fprintf(stderr, "Error: sparse_mat_add(), NULL pointer detected in ");
    fprintf(stderr, "input arguments.\n");
    #endif
    return 1;
  }

  if(((row+1) > mat->row_num) || ((col+1) > mat->col_num)) {
    #ifdef DEBUG
    fprintf(stderr, "Error: sparse_mat_add(), element indices exceed ");
    fprintf(stderr, "matrix dimensions.\n");
    #endif
    return 1;
  }

  if(mat->nz == mat->nz_max) {
    mat_i = realloc(mat->i, 2*mat->nz_max*sizeof(int));
    if(mat_i)
      mat->i = mat_i;
    else {
      #ifdef DEBUG
      fprintf(stderr, "Error: sparse_mat_add(), 'mat_i', realloc() failed.\n.");
      #endif
      return 1;
    }
    mat_j = realloc(mat->j, 2*mat->nz_max*sizeof(int));
    if(mat_j)
      mat->j = mat_j;
    else {
      #ifdef DEBUG
      fprintf(stderr, "Error: sparse_mat_add(), 'mat_j', realloc() failed.\n");
      #endif
      return 1;
    }
    mat_value = realloc(mat->value, 2*mat->nz_max*sizeof(double));
    if(value)
      mat->value = mat_value;
    else {
      #ifdef DEBUG
      fprintf(stderr, "Error: sparse_mat_add(), 'mat_value', ");
      fprintf(stderr, "realloc() failed.\n");
      #endif
      return 1;
    }
    mat->nz_max *= 2;
  }

  mat->value[mat->nz] = value;
  mat->i[mat->nz] = row;
  mat->j[mat->nz] = col;
  mat->nz++;

  return 0;
}

DENSE_MAT* create_dense_mat(int row_num, int col_num) {
  DENSE_MAT* mat;

  mat = malloc(sizeof(DENSE_MAT));
  if(!mat) {
    #ifdef DEBUG
    fprintf(stderr, "Error: create_dense_mat(), 'mat', malloc() failed.\n");
    #endif
    goto failure;
  }

  mat->row_num = row_num;
  mat->col_num = col_num;
  mat->data = malloc(row_num*col_num*sizeof(double));

  if(!mat->data) {
    #ifdef DEBUG
    fprintf(stderr, "Error: create_dense_mat(), 'data', malloc() failed.\n");
    #endif
    goto failure;
  }
  mat->data = malloc(row_num*col_num*sizeof(double));

  return mat;

  failure:
  destroy_dense_mat(mat);

  return NULL;
}

void destroy_dense_mat(DENSE_MAT* mat) {

  if(!mat)
    return;

  if(mat->data)
    free(mat->data);
  free(mat);

  return;
}

void vector_init(double* a) {

  for(int i = 0; i < 3; i ++)
    a[i] = 0.0;

  return;
}

void vector_copy(double* a, const double* b) {

  for(int i = 0; i < 3; i ++)
    a[i] = b[i];

  return;
}

void vector_add(double* res, const double* a, const double* b) {

  for(int i = 0; i < 3; i ++)
    res[i] = a[i] + b[i];

  return;
}

void vector_sub(double* res, const double* a, const double* b) {

  for(int i = 0; i < 3; i ++)
    res[i] = a[i] - b[i];

  return;
}

void vector_scale(double* a, const double x) {

  for(int i = 0; i < 3; i ++)
    a[i] *= x;

  return;
}

double dot_product(const double* a, const double* b) {
  double res;

  res = 0.0;
  for(int i = 0; i < 3; i ++)
    res += a[i]*b[i];

  return res;
}

void cross_product(double* res, const double* a, const double* b) {

  res[0] = a[1]*b[2]-a[2]*b[1];
  res[1] = a[2]*b[0]-a[0]*b[2];
  res[2] = a[0]*b[1]-a[1]*b[0];

  return;
}


void mat_product(int row_num, int col_num, int prod_dim, double(* res)[col_num],
  const double(* a)[prod_dim], const double(* b)[col_num]) {

  for(int i = 0; i < row_num; i ++)
    for(int j = 0; j < col_num; j ++) {
      res[i][j] = 0.0;
      for(int k = 0; k < prod_dim; k ++)
        res[i][j] += a[i][k]*b[k][j];
    }

  return;
}


void mat_transpose(int row_trans, int col_trans, double(* trans_mat)[col_trans],
  const double(* mat)[row_trans]) {

  for(int i = 0; i < row_trans; i ++)
    for(int j = 0; j < col_trans; j ++)
      trans_mat[i][j] = mat[j][i];

  return;
}

void get_rotation_matrix(double p0[3], double pz[3], double px[3],
  double(* rot_mat)[3]) {
  double uz[3], ux[3], uxr[3];
  double phi, theta, psi;
  double cosphi, sinphi;
  double costheta, sintheta;
  double cospsi, sinpsi;
  double eta;

  vector_sub(uz, pz, p0);
  vector_scale(uz, 1.0/sqrt(dot_product(uz, uz)));

  vector_sub(ux, px, p0);
  vector_scale(ux, 1.0/sqrt(dot_product(ux, ux)));

  phi = atan2(uz[0], -uz[1]),
  sinphi = sin(phi);
  cosphi = cos(phi);

  eta = -uz[0]*sinphi+uz[1]*cosphi;
  costheta = uz[2];
  theta = acos(costheta);
  if(eta > 0.0)
    theta *= -1.0;
  sintheta = sin(theta);

  uxr[0] = cosphi;
  uxr[1] = sinphi;
  uxr[2] = 0.0;
  eta = -ux[0]*costheta*sinphi + ux[1]*costheta*cosphi + ux[2]*sintheta;
  cospsi = dot_product(ux, uxr);
  psi = acos(cospsi);
  if(eta < 0.0)
    psi *= -1.0;
  sinpsi = sin(psi);

  rot_mat[0][0] = cospsi*cosphi - costheta*sinphi*sinpsi;
  rot_mat[1][0] = -sinpsi*cosphi - costheta*sinphi*cospsi;
  rot_mat[2][0] = sintheta*sinphi;

  rot_mat[0][1] = cospsi*sinphi + costheta*cosphi*sinpsi;
  rot_mat[1][1] = -sinpsi*sinphi + costheta*cosphi*cospsi;
  rot_mat[2][1] = -sintheta*cosphi;

  rot_mat[0][2] = sinpsi*sintheta;
  rot_mat[1][2] = cospsi*sintheta;
  rot_mat[2][2] = costheta;

  return;
}
