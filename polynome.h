#ifndef DEF_POLY
#define DEF_POLY

#include "common.h"


polynome UN,X,X_1;

extern void moins(polynome *p);

extern void surX(polynome *p);

extern void somme_polynome(polynome p, polynome q,polynome *r);

extern void produit_polynome (polynome p, polynome q,polynome *r);

extern void affiche_polynome (polynome *p);

extern void copie_polynome(polynome p,polynome *q);

extern int egal_polynome(polynome p, polynome q);

extern void init_polynome();

extern void XpuissanceK( int pui, polynome *p);

#endif
