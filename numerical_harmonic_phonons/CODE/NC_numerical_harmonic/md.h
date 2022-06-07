#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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

#define NT         100000         // number of trajectories
#define nts        10000
#define dt         1e-16   // in s
#define M          5            // parameter for integrator
#define KB         1.380649e-23    // J/K
#define HBAR       1.054571817e-34   //SI unit: Js
#define T          300.0         // temperature in K
#define PI         3.14159265358979323846
#define TWOPI      6.28318530717958647692
 
#define sqr(x)     ((x) * (x))

double ran();
double randv();
double randr(par_st);

int read_ph_number(char *);
void read_ph_freq(char *, int, par_st *); 
void read_coupling(char *, int, par_st *); 

void gauss_rand(double *, double *);

void force(double *,sym_st *,sym_st *, par_st *,int);

double energy(sym_st *, sym_st *, par_st *, int);

void init(sym_st *, sym_st *, par_st *, int);
void nerror(char *);
void symplectic(double *, double *);
void symp_prop(sym_st *, sym_st *, double *, double *, par_st *, int);

double delta(sym_st *, double *, par_st *, int);
