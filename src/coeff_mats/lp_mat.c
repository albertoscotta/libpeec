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

#include <stdio.h>
#include "../mesh.h"
#include "../utils/lin_algebra.h"
#include "../constants.h"

#ifdef MACKENZIE

#ifdef RECT_APPROX

int fill_lp_mat(const MESH* mesh, int face_num,
  const double(* integrals)[face_num][3], const DENSE_MAT* lp_mat) {

  if(!(mesh && mesh->node && mesh->iedge && mesh->face && integrals &&
    lp_mat && lp_mat->data)) {
      #ifdef DEBUG
      fprintf(stderr, "Error: fill_lp_mat(), NULL pointer detected in ");
      fprintf(stderr, "input arguments.\n");
      #endif
      return 1;
  }

  double(* lp_mat_data)[mesh->iedge_num];
  double edge_one_vector[3], edge_two_vector[3];
  double midpoint_center_vector[3];  
  double iedge_vector[3];
  double face_normal_vector[3];
  double iedge_normal_unit_vector[3];
  double a[3];

  lp_mat_data = (double(*)[mesh->iedge_num])lp_mat->data;

  #pragma omp parallel for
  for(int i = 0; i < mesh->iedge_num; i ++)
    for(int j = 0; j < mesh->iedge_num; j ++)
      lp_mat_data[i][j] = 0.0;

  #pragma omp parallel for private(edge_one_vector, edge_two_vector, \
    iedge_vector, face_normal_vector, iedge_normal_unit_vector, a, \
    midpoint_center_vector)
  for(int i = 0; i < face_num; i ++)
    for(int j = 0; j < 3; j ++) { /* loop over interior edges */
      if(mesh->face[i].iedge[j] != -1) {

        vector_sub(edge_one_vector,
          (double *)&(mesh->node[mesh->face[i].node[(2-j)/2]]),
          (double *)&(mesh->node[mesh->face[i].node[j]]));
        vector_sub(edge_two_vector,
          (double *)&(mesh->node[mesh->face[i].node[(5-j)/2]]),
          (double *)&(mesh->node[mesh->face[i].node[j]]));

        vector_sub(iedge_vector,
          (double *)&(mesh->node[mesh->iedge[mesh->face[i].iedge[j]].node[1]]),
          (double *)&(mesh->node[mesh->iedge[mesh->face[i].iedge[j]].node[0]]));

        cross_product(face_normal_vector, edge_one_vector, edge_two_vector);
        cross_product(iedge_normal_unit_vector, iedge_vector,
          face_normal_vector);

        vector_scale(iedge_normal_unit_vector, 1.0/(2.0*mesh->face[i].area*
          mesh->iedge[mesh->face[i].iedge[j]].length));

        if(mesh->iedge[mesh->face[i].iedge[j]].face[1] == i)
          vector_scale(iedge_normal_unit_vector, -1.0);

        for(int k = 0; k < face_num; k ++) {

          vector_copy(a, iedge_normal_unit_vector);
          vector_scale(a, MU_0/(4.0*M_PI*
            mesh->iedge[mesh->face[i].iedge[j]].length)*integrals[i][k][j]);

          for(int l = 0; l < 3; l ++)
            if(mesh->face[k].iedge[l] != -1) {

              vector_sub(midpoint_center_vector,
                mesh->iedge[mesh->face[k].iedge[l]].midpoint,
                mesh->face[k].center);

              if(mesh->iedge[mesh->face[k].iedge[l]].face[1] == k)
                vector_scale(midpoint_center_vector, -1.0);

              #pragma omp atomic
              lp_mat_data[mesh->face[k].iedge[l]][mesh->face[i].iedge[j]] +=
                dot_product(a, midpoint_center_vector);
            }
        }
      }
    }

  return 0;
}

#else

int fill_lp_mat(const MESH* mesh, int face_num,
  const double(* integrals)[face_num][3], const DENSE_MAT* lp_mat) {

  if(!(mesh && mesh->node && mesh->iedge && mesh->face && integrals &&
    lp_mat && lp_mat->data)) {
      #ifdef DEBUG
      fprintf(stderr, "Error: fill_lp_mat(), NULL pointer detected in ");
      fprintf(stderr, "input arguments.\n");
      #endif
      return 1;
  }

  double(* lp_mat_data)[mesh->iedge_num];
  double edge_one_vector[3], edge_two_vector[3];
  double midpoint_center_unit_vector[2][3], midpoint_center_vector_length[2];  
  double iedge_vector[3];
  double face_normal_vector[3];
  double iedge_normal_unit_vector[3];
  double a[2][3];  

  lp_mat_data = (double(*)[mesh->iedge_num])lp_mat->data;

  #pragma omp parallel for
  for(int i = 0; i < mesh->iedge_num; i ++)
    for(int j = 0; j < mesh->iedge_num; j ++)
      lp_mat_data[i][j] = 0.0;

  #pragma omp parallel for private(edge_one_vector, edge_two_vector, \
    iedge_vector, face_normal_vector, iedge_normal_unit_vector, a, \
    midpoint_center_unit_vector, midpoint_center_vector_length)
  for(int i = 0; i < face_num; i ++)
    for(int j = 0; j < 3; j ++) /* loop over interior edges */
      if(mesh->face[i].iedge[j] != -1) {

        vector_sub(edge_one_vector,
          (double *)&(mesh->node[mesh->face[i].node[(2-j)/2]]),
          (double *)&(mesh->node[mesh->face[i].node[j]]));
        vector_sub(edge_two_vector,
          (double *)&(mesh->node[mesh->face[i].node[(5-j)/2]]),
          (double *)&(mesh->node[mesh->face[i].node[j]]));

        vector_sub(iedge_vector,
          (double *)&(mesh->node[mesh->iedge[mesh->face[i].iedge[j]].node[1]]),
          (double *)&(mesh->node[mesh->iedge[mesh->face[i].iedge[j]].node[0]]));

        cross_product(face_normal_vector, edge_one_vector, edge_two_vector);
        cross_product(iedge_normal_unit_vector, iedge_vector,
          face_normal_vector);

        vector_scale(iedge_normal_unit_vector, 1.0/(2.0*mesh->face[i].area*
          mesh->iedge[mesh->face[i].iedge[j]].length));

        if(mesh->iedge[mesh->face[i].iedge[j]].face[1] == i)
          vector_scale(iedge_normal_unit_vector, -1.0);

        for(int k = 0; k < mesh->iedge_num; k ++) {
          for(int l = 0; l < 2; l ++) {

          vector_copy(a[l], iedge_normal_unit_vector);
          vector_scale(a[l], MU_0/(4.0*M_PI*
            mesh->iedge[mesh->face[i].iedge[j]].length)*integrals[i][
            mesh->iedge[k].face[l]][j]);
  
            if(l)
              vector_sub(midpoint_center_unit_vector[l],
                mesh->face[mesh->iedge[k].face[l]].center,
                mesh->iedge[k].midpoint);
            else
              vector_sub(midpoint_center_unit_vector[l],
                mesh->iedge[k].midpoint,
                mesh->face[mesh->iedge[k].face[l]].center);

            midpoint_center_vector_length[l] = sqrt(dot_product(
              midpoint_center_unit_vector[l], midpoint_center_unit_vector[l]));
  
            vector_scale(midpoint_center_unit_vector[l],
              1.0/midpoint_center_vector_length[l]);
          }

          #pragma omp atomic 
          lp_mat_data[k][mesh->face[i].iedge[j]] +=
            (dot_product(a[0], midpoint_center_unit_vector[0]) +
             dot_product(a[1], midpoint_center_unit_vector[1]))*
            (midpoint_center_vector_length[0] +
             midpoint_center_vector_length[1])/2.0;
        }
      }

  return 0;
}

#endif

#else

#ifdef RECT_APPROX

int fill_lp_mat(const MESH* mesh, int face_num,
  const double(* integrals)[face_num][3], const DENSE_MAT* lp_mat) {

  if(!(mesh && mesh->node && mesh->iedge && mesh->face && integrals &&
    lp_mat && lp_mat->data)) {
      #ifdef DEBUG
      fprintf(stderr, "Error: fill_lp_mat(), NULL pointer detected in ");
      fprintf(stderr, "input arguments.\n");
      #endif
      return 1;
  }

  double(* lp_mat_data)[mesh->iedge_num];
  double edge_one_vector[3], edge_two_vector[3];
  double midpoint_center_vector[3];
  double a_one[3], a_two[3], a[3];

  lp_mat_data = (double(*)[mesh->iedge_num])lp_mat->data;

  #pragma omp parallel for
  for(int i = 0; i < mesh->iedge_num; i ++)
    for(int j = 0; j < mesh->iedge_num; j ++)
      lp_mat_data[i][j] = 0.0;

  #pragma omp parallel for private(edge_one_vector, edge_two_vector, a_one, \
    a_two,a, midpoint_center_vector)
  for(int i = 0; i < face_num; i ++)
    for(int j = 0; j < 3; j ++) { /* loop over interior edges */
      if(mesh->face[i].iedge[j] != -1) {

        vector_sub(edge_one_vector,
          (double *)&(mesh->node[mesh->face[i].node[(2-j)/2]]),
          (double *)&(mesh->node[mesh->face[i].node[j]]));
        vector_sub(edge_two_vector,
          (double *)&(mesh->node[mesh->face[i].node[(5-j)/2]]),
          (double *)&(mesh->node[mesh->face[i].node[j]]));

        if(mesh->iedge[mesh->face[i].iedge[j]].face[1] == i) {
          vector_scale(edge_one_vector, -1.0);
          vector_scale(edge_two_vector, -1.0);
        }

        for(int k = 0; k < face_num; k ++) {

          vector_copy(a_one, edge_one_vector);
          vector_scale(a_one, MU_0/(8.0*M_PI*mesh->face[i].area)*
            integrals[i][k][(2-j)/2]);
          vector_copy(a_two, edge_two_vector);
          vector_scale(a_two, MU_0/(8.0*M_PI*mesh->face[i].area)*
            integrals[i][k][(5-j)/2]);
          vector_add(a, a_one, a_two);

          for(int l = 0; l < 3; l ++)
            if(mesh->face[k].iedge[l] != -1) {

              vector_sub(midpoint_center_vector,
                mesh->iedge[mesh->face[k].iedge[l]].midpoint,
                mesh->face[k].center);

              if(mesh->iedge[mesh->face[k].iedge[l]].face[1] == k)
                vector_scale(midpoint_center_vector, -1.0);

              #pragma omp atomic
              lp_mat_data[mesh->face[k].iedge[l]][mesh->face[i].iedge[j]] +=
                dot_product(a, midpoint_center_vector);
            }
        }
      }
    }

  return 0;
}


#else

int fill_lp_mat(const MESH* mesh, int face_num,
  const double(* integrals)[face_num][3], const DENSE_MAT* lp_mat) {

  if(!(mesh && mesh->node && mesh->iedge && mesh->face && integrals &&
    lp_mat && lp_mat->data)) {
      #ifdef DEBUG
      fprintf(stderr, "Error: fill_lp_mat(), NULL pointer detected in ");
      fprintf(stderr, "input arguments.\n");
      #endif
      return 1;
  }

  double(* lp_mat_data)[mesh->iedge_num];
  double edge_one_vector[3], edge_two_vector[3];
  double midpoint_center_unit_vector[2][3], midpoint_center_vector_length[2];
  double a_one[3], a_two[3], a[2][3];

  lp_mat_data = (double(*)[mesh->iedge_num])lp_mat->data;

  #pragma omp parallel for
  for(int i = 0; i < mesh->iedge_num; i ++)
    for(int j = 0; j < mesh->iedge_num; j ++)
      lp_mat_data[i][j] = 0.0;

  #pragma omp parallel for private(edge_one_vector, edge_two_vector, a_one, \
    a_two, a, midpoint_center_unit_vector, midpoint_center_vector_length)
  for(int i = 0; i < face_num; i ++)
    for(int j = 0; j < 3; j ++) { /* loop over interior edges */
      if(mesh->face[i].iedge[j] != -1) {

        vector_sub(edge_one_vector,
          (double *)&(mesh->node[mesh->face[i].node[(2-j)/2]]),
          (double *)&(mesh->node[mesh->face[i].node[j]]));
        vector_sub(edge_two_vector,
          (double *)&(mesh->node[mesh->face[i].node[(5-j)/2]]),
          (double *)&(mesh->node[mesh->face[i].node[j]]));

        if(mesh->iedge[mesh->face[i].iedge[j]].face[1] == i) {
          vector_scale(edge_one_vector, -1.0);
          vector_scale(edge_two_vector, -1.0);
        }

        for(int k = 0; k < mesh->iedge_num; k++) {
          for(int l = 0; l < 2; l ++) {

            vector_copy(a_one, edge_one_vector);
            vector_scale(a_one, MU_0/(8.0*M_PI*mesh->face[i].area)*
              integrals[i][mesh->iedge[k].face[l]][(2-j)/2]);
            vector_copy(a_two, edge_two_vector);
            vector_scale(a_two, MU_0/(8.0*M_PI*mesh->face[i].area)*
              integrals[i][mesh->iedge[k].face[l]][(5-j)/2]);
            vector_add(a[l], a_one, a_two);
  
            if(l)
              vector_sub(midpoint_center_unit_vector[l],
                mesh->face[mesh->iedge[k].face[l]].center,
                mesh->iedge[k].midpoint);
            else
              vector_sub(midpoint_center_unit_vector[l],
                mesh->iedge[k].midpoint,
                mesh->face[mesh->iedge[k].face[l]].center);
  
            midpoint_center_vector_length[l] = sqrt(dot_product(
              midpoint_center_unit_vector[l], midpoint_center_unit_vector[l]));
  
            vector_scale(midpoint_center_unit_vector[l],
              1.0/midpoint_center_vector_length[l]);
          }

          #pragma omp atomic 
          lp_mat_data[k][mesh->face[i].iedge[j]] +=
            (dot_product(a[0], midpoint_center_unit_vector[0]) +
             dot_product(a[1], midpoint_center_unit_vector[1]))*
            (midpoint_center_vector_length[0] +
             midpoint_center_vector_length[1])/2.0;
        }
      }
    }

  return 0;
}

#endif

#endif
