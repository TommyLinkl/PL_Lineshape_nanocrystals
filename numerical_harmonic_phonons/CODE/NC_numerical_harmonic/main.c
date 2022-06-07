#include "md.h"

void main() {
  FILE *pf, *sd, *testing;
  int i, j, k, Nph, NB; 
  sym_st *r, *p;
  par_st *par;
  double ene, *a, *b, *r0, *cor, I, domega, v, cutoff, slope, omega, omega_0, *sum_re, *sum_im; 
  double **integral;
  clock_t runtime;

  runtime = clock();

  Nph = read_ph_number("w.dat");
  NB = Nph - 6;
  printf("Number of bath modes: %d \n", NB);

  /*** allocating memory ***/
  r        = (sym_st*)calloc(NB,sizeof(sym_st));
  r0       = (double*)calloc(NB,sizeof(double));
  p        = (sym_st*)calloc(NB,sizeof(sym_st));
  par      = (par_st*)calloc(NB,sizeof(par_st));
  a        = (double*)calloc(M,sizeof(double));
  b        = (double*)calloc(M,sizeof(double));
  cor      = (double*)calloc(nts,sizeof(double));
  integral = (double**)calloc(NT,sizeof(double*));
  for (i = 0; i<NT; i++) {
    integral[i] = (double*) calloc(nts, sizeof(double));
  }
  sum_re   = (double*)calloc(nts,sizeof(double));
  sum_im   = (double*)calloc(nts,sizeof(double));
  for (i=0; i<NT; i++) {
    integral[i][0] = 0;
  }

  read_ph_freq("w.dat", NB, par);
  read_coupling("Vklq-diabatic.dat", NB, par);

  /*** Initialize other quantities ***/
  for (i=0; i<NB; i++) {
    par[i].mass = 1.0;
    par[i].k = par[i].mass * sqr(par[i].omega);
  }

  symplectic(a,b);

  printf("DONE initialization \n");

  /*** main propagation loop ***/
  for (j = 0; j < NT; j++) {
    // printf(" j = %d (max to NT) \n", j);
    init(r,p,par,NB); // configuration
    
    for (i = 1; i < nts; i++) { 
      integral[j][i] = integral[j][i-1] + delta(r,r0,par,NB) * dt;
      symp_prop(r,p,a,b,par,NB);
    }
    if (j%100==0) {
      printf("DONE with traj #j = %d (max to %d) \n", j, NT);
    }
  }

  printf("###################\n##  DONE calculating integrals\n###################\n");

  for (i = 0; i < nts; i++) {
    sum_re[i] = 0;
    sum_im[i] = 0;
  }

  /*
  printf("Integral at nts = 0; for the 0th traj = %g \n", integral[0][0]);
  printf("Integral at nts = 0; for the 1th traj = %g \n", integral[1][0]);
  printf("Integral at nts = 0; for the 2th traj = %g \n", integral[2][0]);
  printf("Integral at nts = 1; for the 0th traj = %g, cos of which = %g \n", integral[0][1], cos(integral[0][1]/HBAR));
  printf("Integral at nts = 1; for the 1th traj = %g, cos of which = %g \n", integral[1][1], cos(integral[1][1]/HBAR));
  printf("Integral at nts = 1; for the 2th traj = %g, cos of which = %g \n", integral[2][1], cos(integral[2][1]/HBAR));
  */
  
  /*** Calculate the real and imaginary parts of the sum of the exponentials ***/
  for (i = 0; i < nts; i++) {
    for (j = 0; j < NT; j++) {
      sum_re[i] += cos(integral[j][i] / HBAR);
      sum_im[i] -= sin(integral[j][i] / HBAR);
    }
  }
    
  printf("DONE calculating the dephasing function \n");
  printf("Start with output to corr.dat \n");

  runtime = clock() - runtime; 
  /*** write to file ***/
  pf = fopen("corr.dat", "w");
  fprintf(pf, "#Number of particles / phonon modes : %d\n", NB);
  fprintf(pf, "#Number of trajectories             : %d\n", NT);
  fprintf(pf, "#Temperature                        : %g\n", T);
  fprintf(pf, "#Number of timestamps               : %d\n", nts);
  fprintf(pf, "#Runtime in seconds                 : %f\n", ((double)runtime)/CLOCKS_PER_SEC);
  fprintf(pf, "# \n#The Dephasing function F(t) \n");
  fprintf(pf, "#time t           real_part F_Re(t)          imaginary_part F_Im(t)\n");
  fprintf(pf, "\n");
  for (i = 0; i < nts; i++) {
    fprintf(pf, "%10g            %10g            %10g\n", (double)(i)*dt, sum_re[i]/NT, sum_im[i]/NT);
    fflush(pf);
  }
  fclose(pf);
  
  /*** free allocated memory ***/
  free(r);
  free(r0);
  free(p);
  free(par);
  free(a);
  free(b);
  free(cor);
    for (i = 0; i < NT; i++) { 
    free(integral[i]);
  }
  free(integral);
  free(sum_re);
  free(sum_im);

  exit(0);
}
