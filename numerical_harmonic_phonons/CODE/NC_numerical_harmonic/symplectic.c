#include "md.h"

/* parameters for symplectic integrator */
void symplectic(double *a,double *b)
{
  a[0] = 0.205177661542290;
  a[1] = 0.403021281604210;
  a[2] = -0.120920876338910;
  a[3] = 0.512721933192410;
  a[4] = 0.0;

  b[0] = 0.061758858135626;
  b[1] = 0.338978026553640;
  b[2] = 0.614791307175580;
  b[3] = -0.140548014659370;
  b[4] = 0.125019822794530;
  return;
}

/******************************************************************************/

/* symplectic integrator to update p and r */
void symp_prop(sym_st *r,sym_st *p,double *a,double *b,par_st *par, int NB)
{
  int i, j;
  double F[NB];

  for (i = 0; i < M; i++) {
    for (j = 0; j < NB; j++)
      r[j].new = r[j].old + b[i] * dt * p[j].old / par[j].mass;
    force(F,r,p,par,NB);
    for (j = 0; j < NB; j++) {
      p[j].new = p[j].old + a[i] * dt * F[j];
      r[j].old = r[j].new;
      p[j].old = p[j].new;
    }
  }
  return;
}
