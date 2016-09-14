#include "include/common.h"
#include "include/print.h"





void affiche_var(int novar)
{
    int indice=0;
    int no_point=0;

    if (novar >= (input->k+input->d)*input->phi) {
        indice = input->k + (novar - (input->k+input->d)*input->phi)/(input->phi-1)+1;
        no_point = (novar - (input->k+input->d)*input->phi) % (input->phi-1) + input->phi + 1;
    }
    else {
        no_point=(novar/(input->k+input->d))+1;
        indice=novar % (input->k+input->d);
    }
    if (indice<input->k) {
        printf("I");
    }
    else printf("K");
#ifdef LATEX
    printf("_{");
#endif
    if (indice < input->k) {
        printf("%d",indice+1);
    }
    else printf("%d",indice-input->k+1);
#ifdef LATEX
    printf("}");
#endif
    printf("(%d)",no_point);
}





void affiche_var_vide()
{
#ifndef LATEX
    //printf("      ");
#endif
}





void affiche_CL(u64* comb)
{
	int novar;
    int flag=FALSE;

    for (novar=0; novar<NbVar; novar++) {
#ifdef LATEX
        if (novar>0) {
            printf("&");
        }
#endif
        if (present(novar,comb)) {
            if (flag) {
                printf(SIGNESOMME);
            }

            flag = TRUE;
            affiche_var(novar);
        }
        else {
            affiche_var_vide();
        }
    }
}





void affiche_attaque(attaque A)
{
	int i;
#ifdef LATEX
    printf("$\\tiny\\left\\{\\begin{array}{*{%d}{@{}c}}\n",input->k*input->phi*4+2);
#endif
    for (i=0; i<A.nb_egalites; i++) {
        affiche_CL(A.egalites[i]);
#ifdef LATEX
        printf("&=&");
#else
        printf(" = ");
#endif
        printf("%d",A.valeur_egalites[i]);
#ifdef LATEX
        printf("\\\\");
#endif
        printf("\n");
    }
    for (i=0; i<A.nb_non_egalites; i++) {
        affiche_CL(A.non_egalites[i]);
#ifdef LATEX
        printf("&\\ne&");
#else
        printf(" <> ");
#endif
        printf("%d",A.valeur_non_egalites[i]);
#ifdef LATEX
        printf("\\\\");
#endif
        printf("\n");
    }
#ifdef LATEX
    printf("\\end{array}\\right.$\n\n");
#endif
}





void affiche_mat_feistel()
{
//#ifdef DEBUG
    int i=0,j=0,t=0;
    
    for (t=0; t<input->phi ; t++) {
        printf("Point %d\n",t+1);
    
        for (j=0; j<=input->d; j++) {
            printf("Tour %d\n",j);
            for (i=0; i<input->k ; i++) {
                affiche_CL(mat_feistel[i][j][t]);
                printf("\n");
            }
        }
    }
//#endif
}





void affiche_attaque_feistel(attaque A)
{
	int i=0;
	
#ifdef LATEX
    if ((NbVar>9)&&(NbVar<16)) {
        printf("$\\tiny\\left\\{\\begin{array}{*{%d}{@{}c}}\n",NbVar+2);
    }
    else if ((NbVar<=9)) {
        printf("$\\left\\{\\begin{array}{*{%d}{@{}c}}\n",NbVar+2);
    }
    else {
        printf("\\begin{eqnarray*}\n");
    }
#endif

    for (i=0; i<A.nb_egalites; i++)
    {
        affiche_CL(A.egalites[i]);
        
#ifdef LATEX
        printf("&=&");
#else
        printf(" = ");
#endif

        printf("%u",A.valeur_egalites[i]);
        
#ifdef LATEX
        printf("\\\\");
#endif

        printf("\n");
    }
    
    for (i=0; i<A.nb_non_egalites; i++)
    {
        affiche_CL(A.non_egalites[i]);
#ifdef LATEX
        printf("&\\ne&");
#else
        printf(" <> ");
#endif

        printf("%u",A.valeur_non_egalites[i]);

#ifdef LATEX
        printf("\\\\");
#endif

        printf("\n");
    }

#ifdef LATEX
    if (NbVar<16) {
        printf("\\end{array}\\right.$\n\n");
    }
    else {
        printf("\\end{eqnarray*}\n\n");
    }
#endif

}





void affiche_binaire(unsigned long int a, int nb_bit)
{
#ifdef DEBUG
    unsigned long int u=1UL << (nb_bit-1);
    while (u != 0UL) {
        if ((u & a) == 0UL) {
            printf("0");
        }
        else printf("1");
        u = u >> 1;
    }
#endif
}





void affiche_partition(unsigned long int* part)
{
#ifdef DEBUG
    int i=0;
    for (i=0; i<input->phi; i++) {
        affiche_binaire(part[i], input->phi);
        printf("\n");
    }
#endif
}





void affiche_attaque2(attaque A)
{
#ifdef DEBUG
    int i=0,j=0;
    printf("%d variables\n", A.nb_variables);
    for (i=0; i<input->phi ; i++) {
        for (j=0; j<input->d ; j++) {
            affiche_binaire(A.chemin[j][i], input->phi);
            printf(" ");
        }
        printf("\n");
    }
#endif
}
