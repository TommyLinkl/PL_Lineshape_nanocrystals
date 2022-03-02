#include "md.h"

/* initialize position of harmonic oscillator 
 * according to equilibrium statistics */
double randr(par_st par) {
  int i;
  double V, Vexp, q, dq, qmin, qmax, emax;

  dq = 10.0;
  emax = 20.0*KB*T;

  // find width of potential to 20KbT
  qmin = 0.0;
  qmax = 0.0;
  // lower bound
  for (V = 0.5*par.mass*sqr(par.omega * qmin); V < emax; qmin -=dq) {
    V = 0.5*par.mass*sqr(par.omega * qmin);
  }
  // upper bound
  for (V = 0.5*par.mass*sqr(par.omega * qmax); V < emax; qmax +=dq) {
    V = 0.5*par.mass*sqr(par.omega * qmax);
  }

  // sample initial position
  // accepts if Boltzmann factor is less than random number (i.e. Metropolis)
  for (V = 0.1, Vexp = 0.0, i = 0; Vexp < V; i++) {
    q = qmin + (qmax-qmin) * ran(); // pick random position in harmonic oscillator
    Vexp = exp(-0.5*par.mass*sqr(par.omega * q) / (KB * T)); // Boltzmann factor
    V = ran();
  }
  return (q);
}
