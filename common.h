#ifndef DEF_COMMON
#define DEF_COMMON


#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "sms4.h"

//#define DEBUG // You can uncomment this to see more details
//#define BIJECTION // You can uncomment this if round functions are bijections

#define SIZE_INPUT 100
//#define LATEX2
//#define LATEX 1
#define SCILAB


#define TRUE 1
#define FALSE 0

#ifdef LATEX
#define SIGNESOMME "\\oplus "
#else
#define SIGNESOMME "+"
#endif

// Default values for an attack
#define SIZE_SCHEMA 3
#define K 3 
#define PHI 2 
#define D 11
#define FEISTEL_TYPE 1

#define FEISTEL_TYPE1 1 // type 1
#define FEISTEL_TYPE2 2 // type 2
#define FEISTEL_TYPE3 3 // type 3
#define FEISTEL_TYPEEXP 4 // expansif
#define FEISTEL_TYPECON 5 // contractant
#define FEISTEL_TYPEEGFN 6

// Default values in configuration file
#define MAX_REL 2 // A regler pour que MAX_REL*PROC >= 4*PHI*K
#define PROC 64 /* 32 ou 64 bits */
#define MAX_EGALITES 1000
#define MAX_NON_EGALITES 1000
#define MAX_SYS 100
#define MAXSYSAFF 100
#define MAX_PART 203
#define MAX_D 30
#define MAX_PHI 8
#define MAX_DEG 300

#define EXACT TRUE
#define NBTERMES 2

typedef uint8_t u8;
typedef uint32_t u32;
typedef uint64_t u64;






typedef struct {
	int coef[MAX_DEG]; /* coef[0] est la constante du polynome */
	int degre;
} polynome;



typedef struct {
	int feistel_type;
	int d;
	int k;
	int size_schema;
	int phi;
	char** schema;
	
	int max_rel;
	int proc;
	int max_egalites;
	int max_non_egalites;
	int max_sys;
	int max_sysaff;
	int max_part;
	int max_d;
	int max_phi;
	int max_deg;
}Input;


typedef struct {
    int nb_egalites;
    u64** egalites;
    u32* valeur_egalites; // 0 par defaut
    int nb_non_egalites;
    u64** non_egalites;
    u32* valeur_non_egalites; //0 par défaut

	int nb_variables;
    u32** chemin;
} attaque;





//Schema par default si on n'en donne pas.
extern char *schema[SIZE_SCHEMA];

extern int debug;
extern int* indice;
extern int compteur_general;
extern char nomvar[9];
extern unsigned int* val_ineg;

//Management variables for feistel schema
extern u64**** mat_feistel;
extern unsigned long int** liste_partitions; // Liste toutes les partitions
extern int nb_partitions;
extern polynome p_somme,r_int;

extern Input* input;

extern polynome* variance;


extern unsigned int* basei;
extern int* compteur;
extern int NbVar; // Nombres de variables dans le systeme 
extern int ISInd; // TRUE si dans le schema d attaque I et S sont independants (pas de combinaisons lineaire entre I et S 

extern int barre_egalites; // barre a ne pas depasser pour avoir une approximation dans l esperance
extern int barre_egalites_var;





#endif
