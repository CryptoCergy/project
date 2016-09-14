#ifndef DEF_POLY2
#define DEF_POLY2



#include "common.h"

// Tool functions for an u32 array and matrix
u32* init_array32(int size);
void free_array32(u32* array);
u32** init_matrix32(int size1, int size2);
void free_matrix32(u32** matrix, int size1);

// Tool functions for an u64 array and matrix
u64* init_array64(int size);
void free_array64(u64* array);
u64** init_matrix64(int size1, int size2);
void free_matrix64(u64** matrix, int size1);

// Tool functions for an int array and matrix
int* init_array(int size);
void free_array(int* array);
int** init_matrix(int size1, int size2);
void free_matrix(int** matrix, int size1);

// Tool functions for a 2variables polynome
void moins2var(int** p) ;
void somme_polynome2var(int** p, int** q, int** r);
void produit_polynome2var (int** p, int** q, int** r) ;
void affiche_pui(int i, char c) ;
void affiche_polynome2var (int** p);
void init_polynome2var(int** p);
void copie_polynome2var(int** p,int** q) ;



#endif
