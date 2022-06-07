#include "md.h"

/* computes Delta function at each time stamp
 * sum( g * x(t) ) */
double delta(sym_st *r, double *r0, par_st *par, int NB) {
    int i;
    double gx_sum = 0.0;

    for (i = 0; i < NB; i++) {
      // Delta = \sum_\alpha { w^2*g*Q + 1/2*w^2*g^2}
      gx_sum += par[i].gama * r[i].old * par[i].omega * par[i].omega + 0.5 * par[i].omega * par[i].omega * par[i].gama * par[i].gama;
    }
    return (gx_sum);
}
