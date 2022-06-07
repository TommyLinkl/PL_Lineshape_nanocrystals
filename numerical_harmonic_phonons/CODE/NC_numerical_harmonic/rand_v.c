#include "md.h"

/* return a random number from a standard normal */
double randv() {
  double rnd1, rnd2;
  void gauss_rand(double *, double *);

  rnd1 = ran();
  rnd2 = ran();
  gauss_rand(&rnd1, &rnd2);

  return rnd1;
}

/*****************************************************************************/

/* sample a standard normal distribution via Box-Muller transform */
void gauss_rand(double *rnd1, double *rnd2) {
  double tmp1 = *rnd1, tmp2 = *rnd2;

  *rnd1 = sqrt(-2.0 * log(tmp1)) * cos(TWOPI * tmp2);
  *rnd2 = sqrt(-2.0 * log(tmp1)) * sin(TWOPI * tmp2);
  return;
}
