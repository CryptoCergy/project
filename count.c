#include "include/common.h"
#include "include/count.h"




void count3(attaque* A,  int var_left, polynome *p)
{
    //Les seconds membres deviennent inutiles avec cette methode.
    // mais elle marche quand meme avec des seconds membres non nuls.
    int i=0, j=0, c=0;
    polynome q, *r=NULL;
    attaque *B=NULL;

    //---------- PHASE I -----------------------------
    // On repere et supprime les inegalites sans variable.
    // Si une inegalite equivaut a 0<>0 on renvoie 0
    //------------------------------------------------
    for (i= A->nb_non_egalites-1; i>=0; i--) {
        if (no_variable(A->non_egalites[i])) {
            if (A->valeur_non_egalites[i] == 0UL) {
                p->coef[0]=0;p->degre=0;
                return ;
            }
            else delete_non_egalite(i, A);
        }
    }

     
    //---------- PHASE II -----------------------------
    // Cas ou on n a pas d inegalite. On renvoie X^nb_var_restantes
    //------------------------------------------------
     
    if (A->nb_non_egalites == 0) {
        if (var_left >= input->max_deg) {
            perror("Augmenter le degre max des polynomes\n");
            exit(-1);
        }
        XpuissanceK(var_left, p);
        return ;
    }
    //---------- PHASE III -----------------------------
    // On repere une variable utilisee une seule fois
    // alors on supprime l inegalite et
    //    on renvoie (X-1) x nb_sol(nouveau_sys,nb_var_restante-1)
    //------------------------------------------------
     
    for (i=0; i< NbVar; i++) {
        c = -1;
        for (j=0; j<A->nb_non_egalites; j++) {
            if (present(i,A->non_egalites[j])) {
                if (c==-1) c=j; else {c=-1; break;}
            }
        }
        if (c!=-1) break;
    } //on se place sur la premiere variable utilisee  dans une seule inegalite
    if (c!=-1)
	{
        delete_non_egalite(c, A);
        count3(A, var_left-1, &q);
        produit_polynome(X_1, q, p);
        return;
    }

 
    //---------- PHASE IV (nouvelle) -----------------------------
    // Fusion / supression d'arretes de coloriage
    // On remplace une non egalite par une egalite ou par rien du tout
    // On fait la difference entre les deux resultats
    //---------- PHASE V -----------------------------
     
    B = attaque_init(); 
    copie_attaque(A, B);
    //On supprime la derniere non egalite
    c=B->nb_non_egalites-1;
    delete_non_egalite(c,B);
    
    //q=(polynome *)malloc(sizeof(polynome));
    count3(B,var_left,&q);
     
    copie_attaque(A,B);
    //on repere une variable dans la derniere non-egalite que l on va supprimer dans les autres.
    for (i=0; i< NbVar; i++) {
        if (present(i, B -> non_egalites[c])) break;
    }
    if (i==NbVar) {
        perror("Cas impossible, inegalite sans variable\n");
        exit(-1);
    }
    
    for (j=0; j<B-> nb_non_egalites-1; j++) {
        if (present(i,B-> non_egalites[j])) {
            additionne_non_egalite_non_egalite(c,j,B);
        }
    }
    
    delete_non_egalite(c,B);
    r  = (polynome *)malloc(sizeof(polynome));

    count3(B,var_left-1,r);
    moins(r);
    somme_polynome(q,*r,p);

	free(r);
    attaque_free(B);
     
    return;
}





void count(attaque* A, int var_left, polynome *p)
{
	attaque *B=NULL,*C=NULL;
    int i=0,j=0,k=0,ii=0;
    u64** liste = calloc(NbVar, sizeof(u64*));

	for(i=0 ; i < NbVar ; i++)
	{
		liste[i] = calloc(input->max_rel, sizeof(u64));
	}

    B = attaque_init();
    copie_attaque(A, B);

    var_left = var_left - (B->nb_egalites);
    B->nb_egalites = 0;

    for (i=0; i< NbVar; i++) {
        for (k=0; k<input->max_rel; k++) {
            liste[i][k] = ~0U; //liste des variables qui apparaissent en meme temps que la variable i (au depart toutes)
        }
        for (j=0; j<B->nb_non_egalites; j++) {
            if (present(i,B->non_egalites[j])) {
                for (k=0; k<input->max_rel; k++) {
                    liste[i][k] = liste[i][k] & B->non_egalites[j][k];
                }
            }
        }
    }

    for (i=0; i<NbVar-1; i++) {
        for (ii=i+1; ii<NbVar; ii++) {
            if (present(i,liste[ii])&& present(ii,liste[i])) {// i et ii sont toujours presentes ensemble -> on supprime i
                ii = NbVar; // on passe a la variable suivante.
                for (j=0; j<B->nb_non_egalites; j++) {
                    if (present(i,B->non_egalites[j])) {
                        enleve_var(B->non_egalites[j],i);
                    }
                }
            }
        }
    }
    //On va separer le systeme en systemes independants
    //On fabrique donc les groupes de variables connexes
    int nb_groupe=0;
    int dernier_groupe;
    for (j=0; j<B->nb_non_egalites; j++) {
        dernier_groupe=-1;
        for (i=0; i<nb_groupe; i++) {
            if (var_en_commun(B->non_egalites[j],liste[i])) {
                if (dernier_groupe==-1) {
                    fusionne(B->non_egalites[j],liste[i]);
                    dernier_groupe=i;
                }
                else {
                    fusionne(liste[i],liste[dernier_groupe]);
                    for (ii=i; ii<nb_groupe-1; ii++) {
                        copy_liste(liste[ii+1],liste[ii]);
                    }
                    nb_groupe--;
                    i--;
                }
            }
        }
        if (dernier_groupe==-1) {
            copy_liste(B->non_egalites[j],liste[nb_groupe]);
            nb_groupe++;
        }
    }

    if (nb_groupe<=1)
	{
        count3(B,var_left,p);
    }

    else
	{
        polynome *q=NULL;
        q  = (polynome *)malloc(sizeof(polynome));

        p ->degre = 0;
        p ->coef[0]=1;
        C = attaque_init();

        for (i=0; i<nb_groupe; i++)
		{
            j = nb_var_presentes(liste[i]);
            var_left -= j;
            creer_attaque(liste[i],B,C); //copie les inegalites de B contenant

            //au moins une variable de liste[i] dans l attaque C
            count3(C,j,q);
            produit_polynome(*q,*p,p);
        }
        XpuissanceK(var_left, q);
        produit_polynome(*q  , *p , p);

        free(q);
        attaque_free(C);
    }


	for(i=0 ; i < NbVar ; i++)
	{
		free(liste[i]);
	}

    free(liste);
	attaque_free(B);

}





void count_feistel(attaque *A, int var_left, polynome *p)
{ 
    attaque *B=NULL;
    B = attaque_init();
    copie_attaque(A, B);
    var_left = var_left - (B -> nb_egalites);
    B->nb_egalites = 0;
    count3(B, var_left, p);
    attaque_free(B);
	return;
}





int var_en_commun(u64* liste1, u64* liste2)
{
    //Renvoie true si les listes ont une variables en commun
    int k=0;
    for(k=0;k<input->max_rel;k++)
	{
        if ((liste1[k]&liste2[k])!=0) return(TRUE);
    }
    return(FALSE);
}




void fusionne(u64* liste1, u64* liste2)
{
    //fusionne la liste1 avec la liste2. Resultat dans la liste2
    int k=0;
    for (k=0; k<input->max_rel; k++)
	{
        liste2[k] = liste2[k] | liste1[k];
    }
}




 
void copy_liste(u64* liste1, u64* liste2)
{
    // Copie la liste1 sur la liste2 qui est effacee
    int k=0;
    for (k=0; k<input->max_rel; k++)
	{
        liste2[k] = liste1[k];
    }
}





int nb_var_presentes(u64* liste)
{
    int poids=0,k=0;
    unsigned long int a=0;

    poids=0;
    for (k=0; k<input->max_rel; k++) {
        a = liste[k];
        while (a!=0) {
            if ((a & 1UL)!=0){
                poids++;
            }
            a = a >>1;
        }
    }
    return(poids);
}





void creer_attaque(u64* liste, attaque* A, attaque* B)
{
    //copie l'attaque A dans l attaque B a condition que au moins
    //une variable soit dans liste
    int i=0,j=0, ni=0;
    B->nb_egalites = A->nb_egalites;

    for (i=0; i< B->nb_egalites; i++) {
        for (j=0; j<input->max_rel; j++) {
            B->egalites[i][j] = A->egalites[i][j];
        }
        B -> valeur_egalites[i] = A -> valeur_egalites[i];
    }    ni = 0;
 
    for (i=0; i< A->nb_non_egalites; i++) {
        if (var_en_commun(liste, A->non_egalites[i])) {
            for (j=0; j<input->max_rel; j++) {
                B->non_egalites[ni][j] = A->non_egalites[i][j];
            }
            B -> valeur_non_egalites[ni] = A -> valeur_non_egalites[i];
            ni++;
        }
    }
    B->nb_non_egalites = ni;
    B -> nb_variables = A -> nb_variables;
    for (i=0; i<input->d  ; i++) {
        for (j=0; j<input->phi; j++) {
            B->chemin[i][j]=A->chemin[i][j];
        }
    }
}
 



