#include "include/polynome.h"

#define LATEX 1

extern polynome UN,X,X_1;

void init_polynome() {
    UN.degre=0;
    UN.coef[0]=1;
    X.degre = 1;
    X.coef[0]=0;
    X.coef[1]=1;
    X_1.degre = 1;
    X_1.coef[0]=-1;
    X_1.coef[1]=1;
}



//donne l oppose d un polynome
void moins(polynome *p) {
	int i;
	for(i=0; i<=p->degre;i++) {
		p->coef[i] = - p->coef[i];
	}
	return;
}



//divise un polynome par X
void surX(polynome *p) {
	int i;
	if (p->coef[0] !=0) {
		perror("Division d un polynome non nul en 0 par X\n");
		exit(-1);
	}
	if (p->degre > 0) {
		p->degre --;
		for(i=0; i <= p->degre; i++) {
			p->coef[i] = p->coef[i+1];
		}
		return;
	}
	return;
}



void somme_polynome(polynome p, polynome q, polynome *r) {
	int i;
	if (p.degre < q.degre) r->degre = q.degre; else r->degre = p.degre;
	for(i=0;i<=r->degre;i++) (r->coef)[i] = 0;
	for(i=0; i<=p.degre; i++) (r->coef)[i] += (p.coef)[i];
	for(i=0;i<=q.degre;i++) (r->coef)[i] += (q.coef)[i];
	i=r->degre;
	while( i>0 && (r->coef)[i] ==0) i--;
	r->degre = i;
	return;
}



void produit_polynome (polynome p, polynome q,polynome *r) {
	int i,j,k,dp,dq;
	dp=p.degre; dq=q.degre;

	if (dp+dq >= input->max_deg) {
		perror("Augmenter MAX_DEG\n");
		exit(-1);
	}
	r->degre=dp+dq;
	for (i=0 ; i<=r->degre; i++) r->coef[i]=0;
	for(j=0; j<=dp; j++) 
		for(k=0; k<=dq; k++)
			(r->coef)[j+k] += (p.coef)[j] * (q.coef)[k];
	return;
}



void affiche_polynome (polynome *p) {
	int i,j=0;
#ifdef LATEX2
    printf("$");
#endif
	if  (p -> degre == 0) {
		printf("%d",p -> coef[0]);
    }
    else {
		
        if (p->coef[0]!=0) {
            printf("%d",p->coef[0]);
            j=1;
        }
        for (i=1; i<= p->degre; i++) {
            if (p->coef[i]>0 && j==1) {
                printf("+");
            }
            if (p->coef[i] != 1 && p->coef[i] !=0) {
                if (p->coef[i] == -1) printf("-");
#ifdef SCILAB
                else printf("%d*",p->coef[i]);
#else
                else printf("%d",p->coef[i]);
#endif
            }
            if (p->coef[i] !=0) {
                j=1;
                printf("N");
#ifdef LATEX2
                if (i > 1) printf("^{%d}",i);
#elif defined(SCILAB)
				if(i > 1) printf("^(%d)",i);
#else
                if (i > 1) printf("^%d",i);
#endif
            }
        }
    }

#ifdef LATEX2
    printf("$");
#endif
    printf(" ");
}



void copie_polynome(polynome p, polynome *q) {
	int i;
	q->degre=p.degre;

	for(i=0;i<=p.degre;i++)
	{
		q->coef[i]=p.coef[i];
	}

	return;
}



int egal_polynome(polynome p, polynome q) {
	int i=0;
	if (p.degre != q.degre) return(0);
	for(i=0; (i<=p.degre) && (p.coef[i] == q.coef[i]); i++) {}
	if (i == p.degre +1) return(1); else return(0);
}



void XpuissanceK(int pui, polynome *p) {
    int i=0;
    p -> degre = pui;
    if (pui >= input->max_deg) {
		perror("Augmenter MAX_DEG\n");
		exit(-1);
	}

	for (i=0; i<pui; i++) {
        (p -> coef)[i]=0;
    }

    (p -> coef)[pui] = 1;
}


