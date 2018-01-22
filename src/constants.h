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

#ifndef _CONSTANTS
#define _CONSTANTS

#include <math.h>

/** \f$\pi\f$. */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
/** \f$\varepsilon_0\f$, vacuum permittivity [\f$\mathrm{\frac{F}{m}}\f$]. */
#define EPSILON_0 (1e-9/(36.0*M_PI))
/** \f$\mu_0\f$, vacuum permeability [\f$\mathrm{\frac{H}{m}}\f$]. */
#define MU_0 (4.0*M_PI*1e-7)

/* Data structures */

/**
 * \brief Contains physical and geometrical constants.
 */
typedef struct {
  double epsilon_r; /**< Relative dielectric constant. */
  double sigma; /**< Electrical conductivity [S]. */
  double delta; /**< Conductive surface thickness [m]. */
} CONSTS;

/* Functions */

/**
 * \brief Sets constants.
 * 
 * \param[in] epsilon_r Relative dielectric constant.
 * \param[in] sigma Electrical conductivity [S].
 * \param[in] delta Conductive surface thickness [m].
 * \param[in,out] cs Structure to be filled with constants.
 */
void set_constants(double epsilon_r, double sigma, double delta, CONSTS* cs);

#endif
