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
#include "coeff_mats.h"
#include "constants.h"
#include "mesh.h"
#include "utils/lin_algebra.h"
#include "utils/integrals.h"
#include "coeff_mats/rp_mat.h"
#include "coeff_mats/pp_mat.h"
#include "coeff_mats/lp_mat.h"

COEFF_MATS* build_coeff_mats(CONSTS cs, const MESH* mesh) {

  if(!mesh) {
    #ifdef DEBUG
    fprintf(stderr, "Error: build_coeff_mats(), 'mesh' NULL pointer.\n");
    #endif
    return NULL;
  }

  COEFF_MATS* coeff_mats;
  double(* integrals)[mesh->face_num][3];

  coeff_mats = NULL;
  integrals = NULL;

  coeff_mats = malloc(sizeof(COEFF_MATS));
  if(!coeff_mats) {
    #ifdef DEBUG
    fprintf(stderr, "Error: build_coeff_mats(), 'coeff_mats', ");
    fprintf(stderr, "malloc() failed.\n");
    #endif
    goto failure;
  }

  coeff_mats->rp = create_sparse_mat(mesh->iedge_num, mesh->iedge_num);
  coeff_mats->pp = create_dense_mat(mesh->face_num, mesh->face_num);
  coeff_mats->lp = create_dense_mat(mesh->iedge_num, mesh->iedge_num);

  if(fill_rp_mat(cs, mesh, coeff_mats->rp))
    goto failure;

  integrals = malloc(mesh->face_num*mesh->face_num*3*sizeof(double));

  compute_integrals(mesh, mesh->face_num, integrals);

  if(fill_pp_mat(cs.epsilon_r, mesh, mesh->face_num, integrals, coeff_mats->pp))
    goto failure;
  if(fill_lp_mat(mesh, mesh->face_num, integrals, coeff_mats->lp))
    goto failure;

  free(integrals);

  return coeff_mats;

  failure:
  if(integrals)
    free(integrals);
  destroy_coeff_mats(coeff_mats);

  return NULL;
}

void destroy_coeff_mats(COEFF_MATS* coeff_mats) {

  if(!coeff_mats)
    return;

  if(coeff_mats->rp)
    destroy_sparse_mat(coeff_mats->rp);
  if(coeff_mats->pp)
    destroy_dense_mat(coeff_mats->pp);
  if(coeff_mats->lp)
    destroy_dense_mat(coeff_mats->lp);
  free(coeff_mats);

  return;
}
