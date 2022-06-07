#include "md.h"

/* set initial equilibrium distributions of solute and bath particles */
void init(sym_st *r, sym_st *p, par_st *par, int NB)
{
  int k;
  double c;

  for (k = 0; k < NB; k++) {
    // printf("k = %d \n", k);
    // initialize position
    r[k].old = randr(par[k]);
    // printf("DONE r_init \n");

    // initialize momentum
    c = sqrt(KB * T);
    p[k].old = c*randv();
    // printf("DONE p_init \n");
  }
  return;
}
