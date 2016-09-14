#ifndef DEF_PRINT
#define DEF_PRINT


#include "variance_random.h"


void affiche_var(int novar);
void affiche_var_vide();
void affiche_CL(u64* comb);
void affiche_attaque(attaque A);


void affiche_mat_feistel();
void affiche_attaque_feistel(attaque A);
void affiche_binaire(unsigned long int a, int nb_bit);
void affiche_partition(unsigned long int* part);
void affiche_attaque2(attaque A);

#endif
