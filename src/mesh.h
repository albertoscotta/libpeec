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

#ifndef _MESH
#define _MESH

/* Data structures */

/**
 * \brief Node coordinates.
 */
typedef struct {
  double x;
  double y;
  double z;
} NODE;

/**
 * \brief Edge topological and geometrical properties.
 */
typedef struct {
  int node[2]; /**< Sorted node indices of the extrema. */
  int face[2]; /**< Sorted face indices of the neighboring triangles. */
  double midpoint[3]; /**< Coordinates of the midpoint. */
  double length; /**< Length. */
} EDGE;

/**
 * \brief Face topological and geometrical properties.
 */
typedef struct {
  int node[3]; /**< Sorted node indices defining the triangular face. */
  int iedge[3]; /**< Neighboring interior edge indices, \c -1 if the ith
                 * edge is a boundary one.
\verbatim
Local edge-to-node connectivity matrix:

edge index -> node indices
0 -> 1, 2
1 -> 0, 2
2 -> 0, 1
\endverbatim
                  */
  double center[3]; /**< Coordinates of the barycenter. */
  double area; /**< Area. */
} FACE;

/**
 * \brief Mesh topological and geometrical properties.
 */
typedef struct {
  int node_num; /**< Number of nodes. */
  NODE* node; /**< Node coordinates. */
  int iedge_num; /**< Number of interior edges. */
  EDGE* iedge; /**< Interior edges topological and geometrical properties. */
  int face_num; /**< Number of faces. */
  FACE* face; /**< Face topological and geometrical properties. */
} MESH;

/* Functions */

/**
 * \brief Builds the structure containing mesh and mesh properties.
 *
 * Computes, additional to those supplied, topological and geometrical
 * properties.
 *
 * \param[in] node_num Number of nodes.
 * \param[in] face_num Number of faces.
 * \param[in] node Vector of node coordinates.
 * \param[in] face2node Face-to-node connectivity matrix.
 * \return Pointer to a \c MESH structure on success, \c NULL otherwise.
 */
MESH* build_mesh(int node_num, int face_num, const double(* node)[3],
  const int(* face2node)[3]);

/**
 * \brief Destroys a \c MESH instance.
 *
 * Frees memory.
 *
 * \param[in] mesh \c MESH instance to be freed.
 */
void destroy_mesh(MESH* mesh);

#endif
