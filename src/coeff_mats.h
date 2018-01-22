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

#ifndef _COEFF_MATS
#define _COEFF_MATS

#include "utils/lin_algebra.h"
#include "constants.h"
#include "mesh.h"

/* Data structures */

/**
 * \brief Contains electrical coefficient matrices.
 */
typedef struct {
  SPARSE_MAT* rp; /**< Resistive coefficient sparse matrix. */
  DENSE_MAT* pp; /**< Potential coefficient dense matrix. */
  DENSE_MAT* lp; /**< Inductive coefficient dense matrix. */
} COEFF_MATS;

/* Functions */

/**
 * \brief Builds the electrical coefficient matrices.
 *
 * \param[in] cs Physical and geometrical constants.
 * \param[in] mesh Mesh, geometrical and topological properties.
 * \return Pointer to a \c COEFF_MATS structure on success, \c NULL otherwise.
 */
COEFF_MATS* build_coeff_mats(CONSTS cs, const MESH* mesh);

/**
 * \brief Destroys a \c COEFF_MATS instance.
 *
 * Frees memory.
 *
 * \param[in] coeff_mats Electrical coefficient matrices.
 */
void destroy_coeff_mats(COEFF_MATS* coeff_mats);

#endif
