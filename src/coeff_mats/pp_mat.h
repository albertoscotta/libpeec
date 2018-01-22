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

#ifndef _PP_MAT
#define _PP_MAT

#include "../constants.h"
#include "../mesh.h"
#include "../utils/lin_algebra.h"

/* Function prototypes */

/**
 * Fills the partial coefficient of potential matrix.
 *
 * \param[in] epsilon_r Relative dielectric constant.
 * \param[in] mesh Mesh and mesh properties.
 * \param[in] face_num Number of faces.
 * \param[in] integrals Integrals from \c compute_integrals().
 * \param[out] pp_mat Pointer to a \c DENSE_MAT instance to be filled with
 *             partial coefficients of potential.
 * \return \c 0 on success, \c 1 on failure.
 */
int fill_pp_mat(double epsilon_r, const MESH* mesh, int face_num,
  const double(* integrals)[face_num][3], DENSE_MAT* pp_mat);

#endif  
