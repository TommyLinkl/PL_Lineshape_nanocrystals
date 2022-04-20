#include "md.h"

/* computes total energy of bath */
double energy(sym_st *r,sym_st *p,par_st *par)
{
  int    i;
  double tmp, ene = 0.0;

  for (i = 0; i < NB; i++){
    // H_b
    tmp = (0.5 * (sqr(p[i].old) / par[i].mass +
		   par[i].mass * sqr(r[i].old * par[i].omega)));
    ene += tmp;
  }
  return (ene);
}
