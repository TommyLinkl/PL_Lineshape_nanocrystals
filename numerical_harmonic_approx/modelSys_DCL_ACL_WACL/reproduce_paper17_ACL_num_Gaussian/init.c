#include "md.h"

double gaussian(double om) {
  return ( 2/sqrt(TWOPI*0.01) * exp(-sqr(om-1)/(2*0.01)) );
}

/******************************************************************************/

/* initialize parameters for solute and bath particles */
void init_par(par_st *par)
{
  int i;
 
  // bath
  double domega = OM_MAX / NB; // 
  for (i = 0; i < NB; i++){
    par[i].omega = (i+1)*domega;
    par[i].mass  = 1.0;
    par[i].k     = par[i].mass * sqr(par[i].omega);
    /* system-bath coupling strength
     * g = sqrt(J(om) * dom) */
    par[i].gama  = sqrt(gaussian(par[i].omega)*domega);
  }
  return;
}

/******************************************************************************/

/* set initial equilibrium distributions of solute and bath particles */
void init(sym_st *r,sym_st *p,par_st *par)
{
  int i;
  double c;

  for (i = 0; i < NB; i++) {

    // initialize position
    r[i].old = randr(par[i]);
    // initialize momentum
    c = sqrt(KB * T);
    p[i].old = c*randv();
  }
  return;
}
