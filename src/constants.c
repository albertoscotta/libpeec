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

#include <stdio.h>
#include "constants.h"

void set_constants(double epsilon_r, double sigma, double delta, CONSTS* cs) {

  if(!cs) {
    #ifdef DEBUG
    fprintf(stderr, "Error: set_constants(), NULL pointer detected in ");
    fprintf(stderr, "input arguments.\n");
    #endif
    return;
  }

  cs->epsilon_r = epsilon_r;
  cs->sigma = sigma;
  cs->delta = delta;

  return;
}
