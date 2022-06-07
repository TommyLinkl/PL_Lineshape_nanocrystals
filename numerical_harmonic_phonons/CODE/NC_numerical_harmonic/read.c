#include "md.h"

/*** Read number of phonons  ***/
int read_ph_number(char phonon_filename[32]) {
  FILE *w;
  int Nph = 0;
  char ch;

  w = fopen(phonon_filename,"r");
  if (w == NULL){
    printf("Error! Can't open file w.dat \n");
    exit(1);
  }
  for (ch = getc(w); ch != EOF; ch = getc(w)) {
    if (ch == '\n') {
      Nph ++;
    }
  }
  printf("Number of phonons: %d \n", Nph);
  fclose(w);

  return Nph;
}


/*** Read phonon frequencies and convert to Hz (omega) ***/
void read_ph_freq(char phonon_filename[32], int NB, par_st *par) {
  FILE *w;
  char read[128], *line; 
  int ph;

  /*** Read the phonon frequencies from w.dat  ***/
  w = fopen(phonon_filename,"r");
  if (w == NULL){
    printf("Error! Can't open file w.dat \n");
    exit(1);
  }
  while (fgets(read, 128, w)!=NULL) {
    line = strtok(read, " ");
    ph = atoi(line);
    if (ph >= 6) {
      line = strtok(NULL, " ");
      par[ph - 6].omega = atof(line) * TWOPI * 1e12;  //frequency converted to Hz (omega) 
    }
    //printf("Phonon frequency: %f \n", atof(line));
  }
  fclose(w);
  printf("DONE READING phonon frequencies \n");

  return;
}


/*** Read the coupling from Vklq.dat, and return the E17 coupling constants in the SI units ***/
void read_coupling(char Vklq_filename[32], int NB, par_st *par) {
  FILE *Vklq;
  char read[128], *line; 
  int ph;

  Vklq = fopen(Vklq_filename,"r");
  if (Vklq == NULL){
    printf("Error! Can't open file Vklq.dat \n");
    exit(1);
  }
  while (fgets(read, 128, Vklq)!=NULL) {
    line = strtok(read, " ");
    ph = atoi(line);
    if (ph >= 6) {
      char *el_1 = strtok(NULL, " ");
      char *el_2 = strtok(NULL, " ");
      char *coupling = strtok(NULL, " ");
      if ((atoi(el_1)==0) && (atoi(el_2)==0)) {
        // printf("Phonon number: %d. Is it 0,0: %d %d. Coupling: %f \n", ph, atoi(el_1), atoi(el_2), atof(coupling));
        par[ph - 6].gama = atof(coupling) / sqr(par[ph - 6].omega);   //SI units now
        // printf("%f /n", par[ph - 6].gama);
      }
    }
  }
  fclose(Vklq);
  printf("DONE READING coupling Vklq \n");
  return;
}