#include "md.h"

/* Ohmic spectral density
 * J(om) = om*exp(-om/om_c) */
double ohmic(double om) {
  //return ( 0.25*om * exp(-om / OM_C) );
  return 0.05 * pow(om, 3.0) * exp(-om);
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
    par[i].gama  = sqrt(ohmic(par[i].omega)*domega);
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
    // initialize momentum for a Boltzmann dist
    // c = sqrt(KB * T);
    // p[i].old = c*randv();
    
    // initialize momentum for a Wigner dist
    c = sqrt(HBAR * par[i].omega / (2 * tanh(HBAR * par[i].omega / (2 * KB * T))));
    p[i].old = c*randv();
  }
  return;
}
