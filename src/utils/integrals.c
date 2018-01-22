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
#include "../mesh.h"

#define TOL 1e-6

int compute_face_integrals(double(* node)[3], int obs_pnt_num,
  double(* obs_pnt)[3], double* i1, double(* in)[3]) {
  double threshold;
  double p0[3], px[3], pz[3];
  double e1[3], e2[3], e3[3];
  double rot_mat[3][3], rot_mat_transposed[3][3];
  double node_trans[3][3], node_transrot[3][3];
  double obs_pnt_trans[3], obs_pnt_transrot[3];
  double edge_vector[3][3], edge_length[3];
  double m[3][3];
  double f2, f3, beta;
  double iu, iv, iua, iva;
  double r_minus, r_plus, r0;
  double t_minus, t_plus, t0;
  double s_minus, s_plus;
  double w0;
  double t0_vect[3];
  double p, pint[3];

  /* local-to-global rotation matrix */
  vector_copy(p0, node[0]); /* origin of local coordinate system */
  vector_copy(px, node[1]); /* point on local x-axis */
  vector_sub(e1, node[1], node[0]);
  vector_sub(e2, node[2], node[0]);

  cross_product(e3, e1, e2);

  threshold = TOL*sqrt(sqrt(dot_product(e3, e3)));

  vector_add(pz, p0, e3); /* point on local z-axis */

  get_rotation_matrix(p0, pz, px, rot_mat);
  mat_transpose(3, 3, rot_mat_transposed, rot_mat);

  /* translate nodes */
  for(int i = 0; i < 3; i ++)
    vector_sub(node_trans[i], node[i], p0);
  /* rotate nodes */
  mat_product(3, 3, 3, node_transrot, node_trans, rot_mat_transposed);

  /* vector of edges */
  for(int i = 0; i < 3; i ++) { /* loop over edges */
    vector_sub(edge_vector[i], node_transrot[(i+1)%3], node_transrot[i]);
    edge_length[i] = sqrt(dot_product(edge_vector[i], edge_vector[i]));
    m[i][0] = edge_vector[i][1]/edge_length[i];
    m[i][1] = -edge_vector[i][0]/edge_length[i];
    m[i][2] = 0.0;
  }

  for(int i = 0; i < obs_pnt_num;  i ++) {
    i1[i] = 0.0;
    iua = 0.0;
    iva = 0.0;

    /* translate observation point */
    vector_sub(obs_pnt_trans, obs_pnt[i], p0);
    /* rotate */
    mat_product(1, 3, 3, (double(*)[3])obs_pnt_transrot,
      (double(*)[3])obs_pnt_trans, rot_mat_transposed);

    w0 = obs_pnt_transrot[2];

    for(int j = 0; j < 3; j ++) { /* loop over edges */
      /* compute p of intersection */
      p = ((obs_pnt_transrot[0] - node_transrot[j][0])*edge_vector[j][0] +
           (obs_pnt_transrot[1] - node_transrot[j][1])*edge_vector[j][1])/
           (edge_length[j]*edge_length[j]);

      /* intersection point */
      pint[0] = node_transrot[j][0] + edge_vector[j][0]*p;
      pint[1] = node_transrot[j][1] + edge_vector[j][1]*p;
      pint[2] = 0.0;

      s_minus = -p*edge_length[j];
      s_plus = (1.0-p)*edge_length[j];

      t_minus =  sqrt((obs_pnt_transrot[0] - node_transrot[j][0])*
        (obs_pnt_transrot[0] - node_transrot[j][0])+
        (obs_pnt_transrot[1] - node_transrot[j][1])*
        (obs_pnt_transrot[1] - node_transrot[j][1]));  
      t_plus = sqrt((obs_pnt_transrot[0] - node_transrot[(j+1)%3][0])*
        (obs_pnt_transrot[0] - node_transrot[(j+1)%3][0])+
        (obs_pnt_transrot[1] - node_transrot[(j+1)%3][1])*
        (obs_pnt_transrot[1] - node_transrot[(j+1)%3][1]));

      vector_sub(t0_vect, pint, obs_pnt_transrot);
      t0 = dot_product(t0_vect, m[j]);

      r0 = sqrt(t0*t0 + w0*w0);
      r_minus = sqrt(t_minus*t_minus + w0*w0);
      r_plus = sqrt(t_plus*t_plus + w0*w0);

      if(fabs(r0) < threshold) { /* point on edge or its extension */
      /* case treated apart because it leads to indeterminate forms */

        if(in != NULL)
          f3 = s_plus*r_plus - s_minus*r_minus;

      } else { /* standard point */
        f2 = log((r_plus + s_plus)/(r_minus + s_minus));

        beta = atan(t0*s_plus/(r0*r0 + fabs(w0)*r_plus)) -
          atan(t0*s_minus/(r0*r0 + fabs(w0)*r_minus));

        i1[i] += t0*f2-fabs(w0)*beta;

        if(in != NULL)
          f3 = s_plus*r_plus - s_minus*r_minus + f2*r0*r0;
      }

      if(in != NULL) {
        iua += 0.5*m[j][0]*f3;
        iva += 0.5*m[j][1]*f3;
      }
    }

    /* integral n/r */
    if(in != NULL) {
      iu = obs_pnt_transrot[0]*i1[i] + iua;
      iv = obs_pnt_transrot[1]*i1[i] + iva;
      in[i][0] = i1[i] - iu/node_transrot[1][0] +
        iv/node_transrot[2][1]*(node_transrot[2][0]/node_transrot[1][0]-1.0);
      in[i][1] = (iu - node_transrot[2][0]*iv/node_transrot[2][1])/
        node_transrot[1][0];
      in[i][2] = iv/node_transrot[2][1];
    }
  }

  return 0;
}

#ifdef MACKENZIE

int compute_integrals(const MESH* mesh, int face_num,
  double(* integrals)[face_num][3]) {
  int status;
  double node[3][3];
  double edge_midpoint[3];
  double(* obs_pnt)[3];
  double(* is)[face_num];

  obs_pnt = malloc(face_num*3*sizeof(double));
  if(!obs_pnt) {
    #ifdef DEBUG
    fprintf(stderr, "Error: compute_integrals(), 'obs_pnt', ");
    fprintf(stderr, " malloc() failed.\n");
    #endif
    if(obs_pnt)
      free(obs_pnt);
    return 1;
  }

  #pragma omp parallel for
  for(int i = 0; i < face_num; i ++)
    vector_copy(obs_pnt[i], mesh->face[i].center);

  status = 0;
  #pragma omp parallel for private(edge_midpoint, node, is)
  for(int i = 0; i < face_num; i ++)
    if(!status) {
      if(!(is = malloc(face_num*4*sizeof(double)))) {
        #pragma omp critical
        status = 1;
      } else {
        for(int j = 0; j < 4; j ++) {
          if(j < 2) {
            vector_add(edge_midpoint,
              (double *)&(mesh->node[mesh->face[i].node[j+1]]),
              (double *)&(mesh->node[mesh->face[i].node[0]]));
            vector_scale(edge_midpoint, 1/2.0);

            vector_copy(node[0], (double *)&(mesh->node[mesh->face[i].node[0]]));
            vector_copy(node[1], mesh->face[i].center);
            vector_copy(node[2], edge_midpoint);
          } else {
            vector_add(edge_midpoint,
              (double *)&(mesh->node[mesh->face[i].node[1]]),
              (double *)&(mesh->node[mesh->face[i].node[2]]));
            vector_scale(edge_midpoint, 1.0/2.0);

            vector_copy(node[0], (double *)&(mesh->node[mesh->face[i].node[0]]));
            vector_copy(node[1], edge_midpoint);
            vector_copy(node[2],
              (double *)&(mesh->node[mesh->face[i].node[j-1]]));
          }
          compute_face_integrals(node, face_num, obs_pnt, is[j], NULL);
        }
        for(int j = 0; j < face_num; j ++) {
          integrals[i][j][0] = is[2][j] + is[3][j] - is[0][j] - is[1][j];
          integrals[i][j][1] = is[3][j] + is[0][j];
          integrals[i][j][2] = is[2][j] + is[1][j];
        }
        free(is);
      }
    }

  if(status) {
    #ifdef DEBUG
    fprintf(stderr, "Error: compute_integrals(), 'is', malloc() failed.\n");
    #endif
  }

  free(obs_pnt);

  return status;
}

#else

int compute_integrals(const MESH* mesh, int face_num,
  double(* integrals)[face_num][3]) {
  int status;
  double node[3][3];
  double(* obs_pnt)[3];
  double* i1;

  if(!mesh || !mesh->node || !mesh->face || !integrals) {
    #ifdef DEBUG
    fprintf(stderr, "Error: compute_integrals(), NULL pointer detected in ");
    fprintf(stderr, "input arguments.\n");
    #endif
    return 1;
  }

  obs_pnt = malloc(face_num*3*sizeof(double));
  if(!obs_pnt) {
    #ifdef DEBUG
    fprintf(stderr, "Error: compute_integrals(), 'obs_pnt', ");
    fprintf(stderr, "malloc() failed.\n");
    #endif
    if(obs_pnt)
      free(obs_pnt);
    return 1;
  }

  #pragma omp parallel for
  for(int i = 0; i < face_num; i ++)
    vector_copy(obs_pnt[i], mesh->face[i].center);

  status = 0;
  #pragma omp parallel for private(i1, node)
  for(int i = 0; i < face_num; i ++)
    if(!status) {
      if(!(i1 = malloc(face_num*sizeof(double)))) {
        #pragma omp critical
        status = 1;
      } else {
        for(int j = 0; j < 3; j ++) /* loop over nodes */
          vector_copy(node[j], (double *)&(mesh->node[mesh->face[i].node[j]]));
        compute_face_integrals(node, face_num, obs_pnt, i1, integrals[i]);

        free(i1);
      }
    }

  if(status) {
    #ifdef DEBUG
    fprintf(stderr, "Error: compute_integrals(), 'i1', malloc() failed.\n");
    #endif
  }
  
  free(obs_pnt);

  return status;
}

#endif
