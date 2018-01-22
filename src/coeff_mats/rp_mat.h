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

#ifndef _RP_MAT
#define _RP_MAT

#include "../constants.h"
#include "../mesh.h"
#include "../utils/lin_algebra.h"

/* Function prototypes */

/**
 * Fills the partial resistance sparse matrix.
 *
 * \param[in] cs Constants containing surface thickness and conductivity.
 * \param[in] mesh Mesh and mesh properties.
 * \param[out] rp_mat Pointer to a \c SPARSE_MAT instance to be filled with
 *             partial resistances.
 * \return \c 0 on success, \c 1 on failure.
 */
int fill_rp_mat(CONSTS cs, const MESH* mesh, SPARSE_MAT* rp_mat);

#endif
