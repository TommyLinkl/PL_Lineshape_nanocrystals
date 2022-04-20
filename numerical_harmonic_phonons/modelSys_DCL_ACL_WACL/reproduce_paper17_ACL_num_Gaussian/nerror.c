#include "md.h"

/* prints error to file */
void nerror(char *str)
{
  FILE *pf;

  pf = fopen("error","w");
  fprintf(pf,"%s\n",str);
  exit(0);
}
