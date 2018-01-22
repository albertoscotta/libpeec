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

int fill_pp_mat(double epsilon_r, const MESH* mesh, int face_num,
  const double(* integrals)[face_num][3], DENSE_MAT* pp_mat) {
  double epsilon;
  double i1;
  double(* pp_mat_data)[face_num];

  if(!(pp_mat && pp_mat->data && integrals && mesh && mesh->face)) {
    #ifdef DEBUG
    fprintf(stderr, "Error: fill_pp_mat(), detected NULL pointer in input ");
    fprintf(stderr, "arguments.\n");
    #endif
    return 1;
  }

  epsilon = epsilon_r*EPSILON_0;
  pp_mat_data = (double(*)[face_num])pp_mat->data;

  #pragma omp parallel for private(i1)
  for(int i = 0; i < face_num; i ++)
    for(int j = 0; j < face_num; j ++) {
      i1 = 0.0;
      for(int k = 0; k < 3; k ++)
        i1 += integrals[j][i][k];
      #ifdef MACKENZIE
      i1 /= 2.0;
      #endif
      pp_mat_data[i][j] = i1/(4.0*M_PI*epsilon*mesh->face[j].area);
    }

  return 0;
}
