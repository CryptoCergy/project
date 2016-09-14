#ifndef DEF_FEISTEL
#define DEF_FEISTEL


#include "variance_random.h"

int representant(u32 lp);
void xor_var(u64* comb, int novar);
void Str2CL_feistel(char *s, u64* comb, int* egalite);
int genere_non_egalites2_feistel(attaque A, int point1, int point2, int no_sys, attaque** res);
int genere_non_egalites_feistel(attaque A, attaque** res);
int NoVar_feistel(int point);

// Tool functions to compute Q(N)
void init_mat_feistel();
int genere_partitions2(int no_part,int no_groupe, unsigned long int* part);
int genere_partitions();
int genere_systeme_feistel2(int no_sys, attaque A, int no_tour, attaque** res3);
int count_and_genere_systeme_feistel2(int no_sys, attaque A, int no_tour);
int genere_systeme_feistel(int no_sys, attaque A, attaque** res3);
int count_and_genere_systeme_feistel(int nb_sys, attaque A);
void manage_feistel_schema();

#endif
