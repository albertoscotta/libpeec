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
#include <math.h>
#include "../constants.h"
#include "../mesh.h"
#include "../utils/lin_algebra.h"

#ifdef MACKENZIE

int fill_rp_mat(CONSTS cs, const MESH* mesh, SPARSE_MAT* rp_mat) {
  double midpoint_center_vector[3];
  double iedge_vector[3];
  double face_normal_vector[3];
  double self_rp;

  if(!mesh && mesh->node && mesh->iedge && mesh->face && rp_mat &&
    rp_mat->i && rp_mat->j && rp_mat->value) {
    #ifdef DEBUG
    fprintf(stderr, "Error: fill_rp_mat(), NULL pointer detected in ");
    fprintf(stderr, "input arguments.\n");
    #endif
    return 1;
  }

  #pragma omp parallel for private(midpoint_center_vector, iedge_vector, \
    face_normal_vector, self_rp)
  for(int i = 0; i < mesh->face_num; i ++)
    for(int k = 0; k < 3; k ++)
      if(mesh->face[i].iedge[k] != -1) {
        vector_sub(midpoint_center_vector,
          mesh->iedge[mesh->face[i].iedge[k]].midpoint, mesh->face[i].center);

        vector_sub(iedge_vector,
          (double *)&(mesh->node[mesh->iedge[mesh->face[i].iedge[k]].node[1]]),
          (double *)&(mesh->node[mesh->iedge[mesh->face[i].iedge[k]].node[0]]));

        cross_product(face_normal_vector, midpoint_center_vector, iedge_vector);

        self_rp = sqrt(dot_product(face_normal_vector, face_normal_vector));
        self_rp /= cs.delta*cs.sigma*dot_product(iedge_vector, iedge_vector);

        #pragma omp critical
        sparse_mat_add(rp_mat, mesh->face[i].iedge[k], mesh->face[i].iedge[k],
          self_rp);
      }

  return 0;
}

#else

int fill_rp_mat(CONSTS cs, const MESH* mesh, SPARSE_MAT* rp_mat) {
  double obs_edge_vector[3], other_edge_vector[3];
  double midpoint_center_vector[3][3];
  double self_rp, mutual_rp;
  double mutual_rp1, mutual_rp2;

  if(!mesh && mesh->node && mesh->iedge && mesh->face && rp_mat &&
    rp_mat->i && rp_mat->j && rp_mat->value) {
    #ifdef DEBUG
    fprintf(stderr, "Error: fill_rp_mat(), NULL pointer detected in ");
    fprintf(stderr, "input arguments.\n");
    #endif
    return 1;
  }

  #pragma omp parallel for private(midpoint_center_vector, self_rp, \
    obs_edge_vector, other_edge_vector, mutual_rp1, mutual_rp2, mutual_rp)
  for(int i = 0; i < mesh->face_num; i ++) {
    for(int k = 0; k < 3; k ++)
      if(mesh->face[i].iedge[k] != -1) {
        vector_sub(midpoint_center_vector[k],
          mesh->iedge[mesh->face[i].iedge[k]].midpoint, mesh->face[i].center);
        if(mesh->iedge[mesh->face[i].iedge[k]].face[1] == i)
          vector_scale(midpoint_center_vector[k], -1.0);
      }

    for(int j = 0; j < 3; j ++) /* src */
      if(mesh->face[i].iedge[j] != -1) {
        for(int k = 0; k < 3; k ++) /* obs */
          if(mesh->face[i].iedge[k] != -1) {
            if(k == j) {
              self_rp = dot_product(midpoint_center_vector[k],
                midpoint_center_vector[k]);
              self_rp *= 5.0/(4.0*cs.delta*cs.sigma*mesh->face[i].area);

              #pragma omp critical
              sparse_mat_add(rp_mat, mesh->face[i].iedge[k],
                mesh->face[i].iedge[j], self_rp);
            } else {
              vector_sub(obs_edge_vector,
                (double *)&(mesh->node[mesh->face[i].node[3-(j+k)]]),
                (double *)&(mesh->node[mesh->face[i].node[j]]));
              vector_sub(other_edge_vector,
                (double *)&(mesh->node[mesh->face[i].node[k]]),
                (double *)&(mesh->node[mesh->face[i].node[j]]));

              if(mesh->iedge[mesh->face[i].iedge[j]].face[1] == i) {
                vector_scale(obs_edge_vector, -1.0);
                vector_scale(other_edge_vector, -1.0);
              }
  
              mutual_rp1 = 5.0/(24.0*cs.delta*cs.sigma*mesh->face[i].area)*
                dot_product(obs_edge_vector, midpoint_center_vector[k]);
              mutual_rp2 = 1.0/(12.0*cs.delta*cs.sigma*mesh->face[i].area)*
                dot_product(other_edge_vector, midpoint_center_vector[k]);

              mutual_rp = mutual_rp1 + mutual_rp2;

              #pragma omp critical
              sparse_mat_add(rp_mat, mesh->face[i].iedge[k],
                mesh->face[i].iedge[j], mutual_rp);
        }
      }
    }
  }

  return 0;
}

#endif
