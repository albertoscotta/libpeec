#include <stdio.h>
#include <peec/mesh.h>
#include <peec/constants.h>
#include <peec/coeff_mats.h>

int main() {
  int status = 0; /* 1, if something goes wrong */
  MESH* mesh = NULL; /* mesh and mesh properties */
  COEFF_MATS* coeff_mats = NULL; /* Rp, Pp, Lp matrices */
  CONSTS cs; /* container for constants */
  int node_num = 8;
  int face_num = 8;
  double node[][3] = {{0.0, 0.0, 0.0},
                      {0.1, 0.0, 0.0},
                      {0.1, 0.3, 0.0},
                      {0.0, 0.3, 0.0},
                      {0.1, 0.15, 0.0},
                      {0.0, 0.15, 0.0},
                      {0.05, 0.075, 0.0},
                      {0.05, 0.225, 0.0}};
  int face2node[][3] = {{0, 1, 6},
                        {4, 5, 6},
                        {2, 3, 7},
                        {4, 7, 5},
                        {0, 6, 5},
                        {1, 4, 6},
                        {2, 7, 4},
                        {3, 5, 7}};
  double epsilon_r = 1.0;
  double sigma = 59.6e6;
  double delta = 1e-3;
  double temp;

  /* build mesh */
  mesh = build_mesh(node_num, face_num, node, face2node);
  if(!mesh) {
    status = 1; /* update status on error */
    printf("Error: couldn't build mesh.\n");
    goto failure; /* go to clean up section */
  }

  /* wrap up constants */
  set_constants(epsilon_r, sigma, delta, &cs);
  /* build Rp, Pp, Lp matrices */ 
  coeff_mats = build_coeff_mats(cs, mesh);
  if(!coeff_mats) {
    status = 1; /* update status on error */
    printf("Error: couldn't build coefficient matrices.\n");
    goto failure; /* go to clean up section */
  }

   /* print results */

  /* Rp */
  printf("Rp mat [%ix%i]:\n", coeff_mats->rp->row_num, coeff_mats->rp->col_num);
  for(int i = 0; i < coeff_mats->rp->row_num; i ++) {
    for(int j = 0; j < coeff_mats->rp->col_num; j ++) {
      temp = 0.0;
      for(int k = 0; k < coeff_mats->rp->nz; k ++)
        if(coeff_mats->rp->i[k] == i && coeff_mats->rp->j[k] == j)
          temp += coeff_mats->rp->value[k];
      printf("%8.1e", temp);
    }
    printf("\n");
  }

  /* Pp */
  printf("Pp mat [%ix%i]:\n", coeff_mats->pp->row_num, coeff_mats->pp->col_num);
  for(int i = 0; i < coeff_mats->pp->row_num; i ++) {
    for(int j = 0; j < coeff_mats->pp->col_num; j ++)
      printf("%8.1e", coeff_mats->pp->data[i*coeff_mats->pp->col_num + j]);
    printf("\n");
  }

  /* Lp */
  printf("Lp mat [%ix%i]:\n", coeff_mats->lp->row_num, coeff_mats->lp->col_num);
  for(int i = 0; i < coeff_mats->lp->row_num; i ++) {
    for(int j = 0; j < coeff_mats->lp->col_num; j ++)
      printf("%8.1e", coeff_mats->lp->data[i*(coeff_mats->lp->col_num) + j]);
    printf("\n");
  }

  failure:
  destroy_coeff_mats(coeff_mats);
  destroy_mesh(mesh);

  return status;
  
}
