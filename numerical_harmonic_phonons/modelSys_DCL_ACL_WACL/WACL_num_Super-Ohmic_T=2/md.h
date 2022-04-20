#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <assert.h>
#include <sys/time.h>
#include <time.h>

/* structure with
 * omega : frequency
 * mass : mass
 * k : mass*omega^2
 * gama : system-bath coupling
 */
struct st1 {
  double omega, mass, k, gama;
};

/* structure with
 * old : variable with last time
 * new : variable at current time
 */
struct st2 {
  double old, new;
};

typedef struct st1 par_st;
typedef struct st2 sym_st;

#define NB         500         // number of bath particles
#define NT         300000         // number of trajectories
#define M          5            // parameter for integrator
#define OM_MAX     20.0         // maximum frequency in bath
#define OM_C       1.0          // characteristic bath frequency
#define KB         1.0          // Boltzmann constant
#define T          2.0         // temperature
#define nts        10000       // number of timestamps
#define PI         3.14159265358979323846
#define TWOPI      6.28318530717958647692
#define HBAR       1    // instead of 1.054571817e-34
 
#define sqr(x)     ((x) * (x))
#define tanh(x)    (exp(2*x) - 1) / (exp(2*x) + 1)

double ran();
double randv();
double randr(par_st);
void gauss_rand(double *, double *);

double energy(sym_st *, sym_st *, par_st *);

void init_par(par_st *);
void init(sym_st *, sym_st *, par_st *);
void nerror(char *);
void symplectic(double *, double *);
void symp_prop(sym_st *, sym_st *, double *, double *, double, par_st *);

double delta(sym_st *, double *, par_st *);
