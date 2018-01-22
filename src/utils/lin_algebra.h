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

#ifndef _LIN_ALGEBRA
#define _LIN_ALGEBRA

/* Data structures */

/** \brief Contains a sparse matrix. */
typedef struct {
  int row_num; /**< Number of rows. */
  int col_num; /**< Number of columns. */
  int nz; /**< Number of non-zero elements. */
  int nz_max; /**< Maximum number of storable non-zero elements, before
               * resizing. */
  int* i; /**< Vector of element row indices. */
  int* j; /**< vector of element column indices. */
  double* value; /**< Vector of element values. */
} SPARSE_MAT;

/** \brief Contains a dense matrix. */
typedef struct {
  int row_num; /**< Number of rows. */
  int col_num; /**< Number of columns. */
  double* data; /**< Flat matrix. */
} DENSE_MAT;

/* Functions */

/**
 * \brief Creates a sparse matrix.
 *
 * \param[in] row_num Number of rows.
 * \param[in] col_cum Number of columns.
 * \return Pointer to a \c SPARSE_MAT structure on success, \c NULL otherwise.
 */
SPARSE_MAT* create_sparse_mat(int row_num, int col_num);

/**
 * \brief Destroys a \c SPARSE_MAT instance.
 *
 * Frees memory.
 *
 * \param[in] mat \c SPARSE_MAT instance to be freed.
 */
void destroy_sparse_mat(SPARSE_MAT* mat);

/**
 * \brief Adds element to a sparse matrix.
 *
 * Adds \c value to element [\c row,\c col].
 *
 * \param[in] mat Sparse matrix.
 * \param[in] row Element row index.
 * \param[in] col Element column index.
 * \param[in] value Element value.
 * \return \c 0 on success, \c 1 on failure.
 */
int sparse_mat_add(SPARSE_MAT* mat, int row, int col, double value);

/**
 * \brief Creates a dense matrix.
 *
 * \param[in] row_num Number of rows.
 * \param[in] col_cum Number of columns.
 * \return Pointer to a \c DENSE_MAT structure on success, \c NULL otherwise.
 */
DENSE_MAT* create_dense_mat(int row_num, int col_num);

/**
 * \brief Destroys a \c DENSE_MAT instance.
 *
 * Frees memory.
 *
 * \param[in] mat \c DENSE_MAT instance to be freed.
 */
void destroy_dense_mat(DENSE_MAT* mat);

/* Auxiliary functions */

/**
 * \brief Initializes a 3D vector.
 * 
 * Makes the vector [0, 0, 0].
 *
 * \param[in, out] a Vector.
 */
void vector_init(double* a);

/**
 * \brief Copies 3D vector.
 *
 * \param[out] a Copy-to vector.
 * \param[in] b Copy-from vector.
 */
void vector_copy(double* a, const double* b); 

/**
 * \brief Adds two 3D vectors.
 *
 * \c res = \c a + \c b.
 *
 * \param[out] res Result.
 * \param[in] a First addend.
 * \param[in] b Second addend.
 */
void vector_add(double* res, const double* a, const double* b);

/**
 * \brief Subtracts a 3D vector from another.
 *
 * \c res = \c a - \c b.
 *
 * \param[out] res Result.
 * \param[in] a Minuend.
 * \param[in] b Subtrahend.
 */
void vector_sub(double* res, const double* a, const double* b);

/**
 * \brief Scales 3D vector.
 *
 * \c a = \c a * \c x
 *
 * \param[in, out] a Vector.
 * \param[in] x Scale factor.
 */
void vector_scale(double* a, const double x);

/**
 * \brief Computes the dot product of two 3D vectors.
 *
 * \param[in] a First vector.
 * \param[in] b Second vector.
 * \result \c a \f$\cdot\f$ \c b.
 */
double dot_product(const double* a, const double* b);

/**
 * \brief Computes the cross product of two 3D vectors.
 *
 * \c res = \c a x \c b.
 *
 * \param[out] res Result vector.
 * \param[in] a First vector.
 * \param[in] b Second vector.
 */
void cross_product(double* res, const double* a, const double* b);

/**
 *
 * \brief Computes the product of two matrices.
 *
 * \c res = \c a * \c b;
 *
 * \param[in] row_num Number of rows of \c res.
 * \param[in] col_num Number of columns of \c res.
 * \param[in] prod_dim Dimension along which the product is performed, i.e.
 *            Number of columns of \c a or Number of rows of \c b.
 * \param[out] res Resultant matrix.
 * \param[in] a First operand.
 * \param[in] b Second operand.
 */
void mat_product(int row_num, int col_num, int prod_dim, double(* res)[col_num],
  const double(* a)[prod_dim], const double(* b)[col_num]);

/**
 * \brief Transposes a matrix.
 *
 * \c trans_mat = \c mat '
 *
 * \param[in] row_trans Number of transpose matrix rows 
 * \param[in] col_trans Number of transpose matrix columns.
 * \param[out] trans_mat Resultant transposed matrix.
 * \param[in] mat Matrix to be transposed.
 */
void mat_transpose(int row_trans, int col_trans, double(* trans_mat)[col_trans],
  const double(* mat)[row_trans]);

/**
 * \brief Computes the rotation matrix between two reference frames.
 *
 * Computes the global-to-local rotation matrix.
 *
 * Local point = \c rot_mat * global point.
 *
 * \param[in] p0 Origin of the local frame.
 * \param[in] pz Point on the local z-axis.
 * \param[in] px Point on the local x-axis.
 * \param[out] rot_mat Rotation matrix [3x3].
 */
void get_rotation_matrix(double p0[3], double pz[3], double px[3],
  double(* rot_mat)[3]);

#endif
