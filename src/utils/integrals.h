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

#ifndef _INTEGRALS
#define _INTEGRALS

#include "../mesh.h"

/* Functions */

/**
 * \brief Computes integrals for the evaluation of Lp and Pp terms.
 *
 * The choice of which integrals are evaluated differ upon compilation.
 * Two options are possible: RWG(default) or MACKENZIE.
 * In both cases, three integrals for face are evaluated, for RWG these are
 * associated with nodes while for MACKENZIE the correspondence is with edges.
 *
 * \param[in] mesh Mesh and mesh properties.
 * \param[in] face_num Number of faces.
 * \param[out] integrals Integrals collected in a 3D matrix [face,
 *  observation point, face node/edge].
 */
int compute_integrals(const MESH* mesh, int face_num,
  double(* integrals)[face_num][3]);

#endif
