#include "md.h"

/* computes forces on bath particles */
void force(double *F,sym_st *r,sym_st *p,par_st *par)
{
  int    i;

  for (i = 0; i < NB; i++){
    // F_b
    F[i]  = -par[i].k * r[i].new - 0.5 * (sqr(par[i].omega) * par[i].gama);
  }
  return;
}





