#ifndef DEF_RAND
#define DEF_RAND

#include "common.h"
#include "count.h"
#include "polynome.h"
#include "polynome2var.h"
#include "print.h"


// Tool functions for an attack
attaque* attaque_init();
void attaque_free(attaque* A);

// Tool functions to compute P(N,m) and P(N)
unsigned int index_base(unsigned int a, unsigned int* base, int dim);
int base(unsigned int* liste, int maxliste ,unsigned int* base);
void polynome2polynome2var(polynome p, int** p2);
int min(int a, int b);
int present(int novar, u64* comb);
void enleve_var(u64* comb, int novar);
void ajoute_toutes_var(u64* comb, u64* comb2);
void ajoute_var(u64* comb, int novar);
void init_comb(u64* comb);
void Str2CL(char *s, u64* comb, int* egalite);
int NoVar(int point, int IouS);
int NoVarp(int point, int IouS);
void copie_attaque(attaque *A, attaque *B);
void ajoute_egalite(u64* comb, unsigned int val, attaque *A);
void ajoute_inegalite(u64* comb, unsigned int val, attaque *A);
void delete_egalite(int no_eg, attaque *A);
void delete_non_egalite(int no_ineg, attaque *A);
void echange(int i, int j, attaque *A);
void additionne_egalite(int i, int j, attaque *A);
void additionne_egalite_non_egalite(int i, int j, attaque *A);
void additionne_non_egalite_non_egalite(int i, int j, attaque *A);
int no_variable(u64* liste);
int solo_variable(u64* liste, int no_var);
int identique_liste_var(u64* liste1, u64* liste2);
int simplifie_systeme( attaque *A);
int genere_non_egalites2(attaque* A, int point1, int point2, int IouS, int no_sys, attaque **res);
int genere_non_egalites(attaque* A, int IouS, attaque **res);
int nb_inegalites(int novar, attaque *A);
int nombre_distincts(unsigned int* liste, int max_liste);
void init_ph(int*** ph);
int Point_present(int noPoint, unsigned int masque);
unsigned int ajoute_point(int noPoint, unsigned int masque);
int hamming(unsigned int u);
int trouve_Point(int i, unsigned int Ipris);
void genere_systemes2(int couple, int h, unsigned int Ipris, unsigned int Ippris, attaque** mix_sys);
void genere_systemes(int noPoint, unsigned int Ipris, unsigned int  Ippris, attaque** mix_sys, attaque** mix_sys2);


#endif
