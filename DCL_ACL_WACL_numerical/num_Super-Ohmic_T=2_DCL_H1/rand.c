#include <sys/time.h>
#include <math.h>
#include <stdio.h>


/********************************************************/
/*							*/
/*	THIS PROGRAM RETURNS A NUMBER ("rnd") 		*/
/*	BETWEEN 0.0 and 1.0 .				*/
/*							*/
/********************************************************/

double ran()
{
  double     rnd;
  int static once = 0;
  
  if(!once) {
    struct timeval tv;
    struct timezone tpz;

    gettimeofday(&tv,&tpz);
    srandom(tv.tv_sec);
    once = 1;
  }
  rnd = (double)random() / 2147483647.0;

  return rnd;
}
