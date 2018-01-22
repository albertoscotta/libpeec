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

#ifndef _LP_MAT
#define _LP_MAT

#include "../mesh.h"
#include "../utils/lin_algebra.h"

/* Function prototypes */

/**
 * Fills the partial inductance matrix.
 *
 * \param[in] mesh Mesh and mesh properties.
 * \param[in] face_num Number.
 * \param[in] integrals Integrals from \c compute_integrals(). 
 * \param[out] pp_mat Pointer to a \c DENSE_MAT instance to be filled with
 *             partial inductances.
 * \return \c 0 on success, \c 1 on failure.
 */
int fill_lp_mat(const MESH* mesh, int face_num,
  const double(* integrals)[face_num][3], const DENSE_MAT* lp_mat);

#endif  
