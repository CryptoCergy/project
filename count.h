#ifndef DEF_COUNT
#define DEF_COUNT

#include "common.h"
#include "print.h"
#include "variance_random.h"


void count3(attaque* A,  int var_left, polynome *p);
void count(attaque *A, int var_left, polynome *p);
void count_feistel(attaque* A, int var_left, polynome *p);

int var_en_commun(u64* liste1, u64* liste2);
void fusionne(u64* liste1, u64* liste2);
void copy_liste(u64* liste1, u64* liste2);
int nb_var_presentes(u64* liste);
void creer_attaque(u64* liste, attaque* A, attaque* B);

#endif
