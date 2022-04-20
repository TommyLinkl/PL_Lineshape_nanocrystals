#include "md.h"

void main()
{
  FILE *pf, *sd, *testing;
  int i, j, k;
  sym_st *r, *p;
  par_st *par;
  double ene, *a, *b, dt, *r0, *cor, I, domega, v, cutoff, slope, omega, omega_0, *sum_re, *sum_im;
  double **integral;
  clock_t runtime;

  runtime = clock();
  /*** allocating memory ***/
  r        = (sym_st*)calloc(NB,sizeof(sym_st));
  r0       = (double*)calloc(NB,sizeof(double));
  p        = (sym_st*)calloc(NB,sizeof(sym_st));
  par      = (par_st*)calloc(NB,sizeof(par_st));
  a        = (double*)calloc(M,sizeof(double));
  b        = (double*)calloc(M,sizeof(double));
  cor      = (double*)calloc(nts,sizeof(double));
  integral = (double**)calloc(NT, sizeof(double*));
  for (i = 0; i < NT; i++) { 
    integral[i] = (double*)calloc(nts, sizeof(double));
  }
  sum_re   = (double*)calloc(nts,sizeof(double));
  sum_im   = (double*)calloc(nts,sizeof(double));
  
  for (i = 0; i < NT; i++) {
    integral[i][0] = 0;
  }

  /*** initial conditions ***/
  init_par(par); // parameters
  
  dt = 0.005;
  symplectic(a,b);

  /*** main propagation loop ***/
  for (j = 0; j < NT; j++) { 
    
    init(r,p,par);  // configuration
        
    for (i = 1; i < nts; i++) { 
      integral[j][i] = integral[j][i-1] + delta(r,r0,par) * dt;
      symp_prop(r,p,a,b,dt,par);
    }
  }

  for (i = 0; i < nts; i++) {
    sum_re[i] = 0;
    sum_im[i] = 0;
  }

  /*** Calculate the real and imaginary parts of the sum of the exponentials ***/
  for (i = 0; i < nts; i++) {
    for (j = 0; j < NT; j++) {
      sum_re[i] += cos(integral[j][i]);
      sum_im[i] -= sin(integral[j][i]);
    }
  }
    
  runtime = clock() - runtime; 
  /*** write to file ***/
  pf = fopen("corr.dat", "w");
  fprintf(pf, "#Number of particles   : %d\n", NB);
  fprintf(pf, "#Number of trajectories: %d\n", NT);
  fprintf(pf, "#Temperature           : %g\n", T);
  fprintf(pf, "#Number of timestamps  : %d\n", nts);
  fprintf(pf, "#Runtime in seconds    : %f\n", ((double)runtime)/CLOCKS_PER_SEC);
  fprintf(pf, "#time       real_part      imaginary_part\n");
  fprintf(pf, "\n");
  for (i = 0; i < nts; i++) { // timestep
    fprintf(pf, "%10g   %10g   %10g\n", (double)(i)*dt, sum_re[i]/NT, sum_im[i]/NT);
    fflush(pf);
  }
  fclose(pf);
  
  /***
  //fourier transform to obtain the spectrum I(omega)
  omega_0 = 1;
  for (omega = -30.0; omega < 30.0; omega+=0.01){ //iterate through some range of omega
    I = 0;
    for (i = 0; i < nts; i++) {
      I += 1/PI * dt * cort_erf[i] * cos(omega*(double)(i)*dt);
    }
    fprintf(sd, "%g %g\n", omega, I);
    fflush(pf);
  }
  fclose(sd);
  ***/
  
  /*** free allocated memory ***/
  free(r);
  free(r0);
  free(p);
  free(par);
  free(a);
  free(b);
  free(cor);
  free(integral);
  for (i = 0; i < NT; i++) { 
    free(integral[i]);
  }
  free(sum_re);
  free(sum_im);

  exit(0);
}
