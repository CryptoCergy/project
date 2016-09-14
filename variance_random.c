#include "include/common.h"
#include "include/variance_random.h"



attaque* attaque_init()
{
	// Allocation
	attaque* A = malloc(sizeof(attaque));
	A->egalites = init_matrix64(input->max_egalites, input->max_rel);
	A->non_egalites = init_matrix64(input->max_non_egalites, input->max_rel);
	A->chemin = init_matrix32(input->max_d, input->max_phi);

	A->valeur_egalites = init_array32(input->max_egalites);
	A->valeur_non_egalites = init_array32(input->max_non_egalites);

	// Initialization
	A->nb_egalites = 0;
	A->nb_non_egalites = 0;
	A->nb_variables = 0;

	return A;
}





void attaque_free(attaque* A)
{
	if(A != NULL)
	{
		free_matrix64(A->egalites, input->max_egalites);
		free_matrix64(A->non_egalites, input->max_non_egalites);
		free_matrix32(A->chemin, input->max_d);
		free_array32(A->valeur_egalites);
		free_array32(A->valeur_non_egalites);
		free(A);
	}

	return;
}




unsigned int index_base(unsigned int a, unsigned int* base, int dim)
{
    // Retourne les coordonnees de a dans la base sous forme d un entier unsigned.
    int j=0;
    unsigned int c=0UL;

    while (a != 0UL && j<dim ) {
        if ((a ^ base[j]) < a) {
            a ^= base[j];
            c ^= (1UL << j);
        }
        j++;
    }
    if (a == 0UL) {
        return(c);
    }
    else {
        perror("Pas dans le sev engendre par la base");
        exit(-1);
    }
}





int base(unsigned int* liste,int maxliste ,unsigned int* base)
{
	// transforme une liste d elements en base triangulaire (du plus gros au plus petit) et retourne la dimension
    int i,j,k, dim;
    unsigned int a;
    dim = 0;

    for (i=0; i<maxliste; i++)
	{
        a = liste[i];

        for (j=0; j<dim; j++)
		{
            if ((a ^ base[j]) < a) a ^= base[j];
        }

        if (a != 0UL) 
		{
            j = 0;
            while ((a < base[j]) && j<dim) {
                j++;
            }

    		for (k=dim; k>j; k--) {
                base[k] = base[k-1];
            }
            dim++;
            base[j] = a;
        }
    }
    return(dim);
}





void polynome2polynome2var ( polynome p, int** p2)
{
	// transforme un polynome en polynome a deux variables
    int i,j;
    for (i=0; i<input->max_deg; i++) {

        for (j=0; j<input->max_deg; j++) {
            p2[i][j] = 0;
        }
    }


    for (i=0; i<= p.degre; i++) {
        p2[i][0] = p.coef[i];
    }
}





int min(int a, int b)
{
    if (a<b) return(a); else return(b);
}





int present(int novar, u64* comb)
{
	if (comb[novar /input->proc] & (u64)(1ULL << (u64)(novar % input->proc))) return(TRUE);
    else return(FALSE);
}





void enleve_var(u64* comb, int novar)
{
    // enleve la variable de la combinaison
    if (present(novar,comb)) {
        comb[novar /input->proc] ^= (u64)(1ULL << (u64)(novar % input->proc));
    }

    else {
        perror("On veut enlever une variable non presente");
        exit(-1);
    }

}





void ajoute_toutes_var(u64* comb, u64* comb2)
{
	int novar = 0;
	for(novar=0 ; novar < NbVar ; novar++)
	{
		if(present(novar, comb2))
		{
			//Si absent on ajoute
			if(!present(novar, comb))
			{
				ajoute_var(comb, novar);
			}

			//Si elle est deja presente, on la retire
			else
			{
				enleve_var(comb, novar);
			}
		}
	}
	return;
}





void ajoute_var(u64* comb, int novar)
{
    comb[novar /input->proc] |= (u64)(1ULL << (u64)(novar % input->proc));
}





void init_comb(u64* comb)
{
	//initialise une combinaison a 0UL
    int i=0;

    for (i=0; i<input->max_rel; i++) {
        comb[i] = 0ULL;
    }
}





void Str2CL(char *s, u64* comb, int* egalite)
{
	int novar=0, i=0;
    int flagI = FALSE; int flagS = FALSE;

    for (i=0; i<input->max_rel; i++) {
        comb[i] = 0UL;
    }

    i=0;
    while (s[i] != '=' && s[i] != '!' && s[i] != '\0') {
        switch (s[i]) {
            case 'I':
                novar = 0;
                flagI = TRUE;
                break;
            case 'S':
                flagS = TRUE;
                novar = input->k * input->phi;
                break;
            default:
                perror("chaine non conforme");
                exit(-1);
                break;
        }
        i++;
        switch (s[i]) {
            case '1':
                break;
            case '2':
                novar +=1*input->phi;
                break;
            case '3':
                novar +=2*input->phi;
                break;
            case '4':
                novar +=3*input->phi;
                break;
            case '5':
                novar +=4*input->phi;
                break;
            case '6':
                novar +=5*input->phi;
                break;
            case '7':
                novar +=6*input->phi;
                break;
            case '8':
                novar +=7*input->phi;
                break;
            case '9':
                novar +=8*input->phi;
                break;
            default:
                perror("chaine non conforme");
                exit(-1);
                break;
        }
        i++;
        switch (s[i]) {
            case '(':
                break;
            default:
                perror("chaine non conforme");
                exit(-1);
                break;
        }
        i++;
        switch (s[i]) {
            case '1':
                break;
            case '2':
                novar +=1;
                break;
            case '3':
                novar +=2;
                break;
            case '4':
                novar +=3;
                break;
            case '5':
                novar +=4;
                break;
            case '6':
                novar +=5;
                break;
            case '7':
                novar +=6;
                break;
            case '8':
                novar +=7;
                break;
            case '9':
                novar +=8;
                break;
            default:
                perror("chaine non conforme");
                exit(-1);
                break;
        }
        i++;
        switch (s[i]) {
            case ')':
                break;
            default:
                perror("chaine non conforme");
                exit(-1);
                break;
        }
        i++;
        switch (s[i]) {
            case '+':
                i++;
                break;
            case '=':
				*egalite = TRUE;
                break;

			case '!':
				*egalite = FALSE;
				break;

            default:
                perror("chaine non conforme");
                exit(-1);
                break;
        }
        ajoute_var(comb, novar);
    }
    if ((flagI == TRUE) && (flagS == TRUE)) ISInd = FALSE;
    return;
}





int NoVar(int point, int IouS)
{
    // IouS = 0 -> I, IouS = 1 -> S
    return( (point - 1) + IouS * input->k * input->phi);
}





int NoVarp(int point, int IouS)
{
    // IouS = 0 -> I, IouS = 1 -> S
    return( (point - 1) + (IouS + 2) * input->k * input->phi);
}





void copie_attaque(attaque *A, attaque *B)
{
	int i=0, j=0;

    B->nb_egalites = A->nb_egalites;
    for (i=0; i< B->nb_egalites; i++) {

        for (j=0; j<input->max_rel; j++) {
            B->egalites[i][j] = A->egalites[i][j];
        }
        B -> valeur_egalites[i] = A -> valeur_egalites[i];
    }
    B->nb_non_egalites = A->nb_non_egalites;
    for (i=0; i< B->nb_non_egalites; i++) {

        for (j=0; j<input->max_rel; j++) {
            B->non_egalites[i][j] = A->non_egalites[i][j];
        }
        B->valeur_non_egalites[i] = A->valeur_non_egalites[i];
    }

 	B->nb_variables = A->nb_variables;
    for (i=0; i<input->d ; i++) {
        for (j=0; j<input->phi; j++) {
            B->chemin[i][j]=A->chemin[i][j];
        }
    }

}





void ajoute_egalite(u64* comb, unsigned int val, attaque *A)
{
	// ajoute comb = val
    int k=0;
    if (A->nb_egalites >= input->max_egalites) {
        perror("Trop d egalites");
        exit(-1);
    }

    for (k=0 ; k<input->max_rel ; k++) {
        A->egalites[A->nb_egalites][k] = comb[k];
    }
    A->valeur_egalites[A->nb_egalites] = val;
    (A->nb_egalites)++;
}





void ajoute_inegalite(u64* comb, unsigned int val, attaque *A)
{
	// ajoute comb = val
    int k;
    if (A->nb_non_egalites >= input->max_non_egalites) {
        perror("Trop d inegalites");
        exit(-1);
    }

    for (k=0 ; k<input->max_rel ; k++) {
        A->non_egalites[A->nb_non_egalites][k] = comb[k];
    }

    A->valeur_non_egalites[A->nb_non_egalites] = val;
    (A->nb_non_egalites)++;
}





void delete_egalite(int no_eg, attaque *A)
{
	//supprime la non egalite numero no_ineg de l attaque A
    int j=0,k=0;
    for (j=no_eg+1 ; j<A->nb_egalites ; j++) {
        for (k=0 ; k<input->max_rel ; k++) {
            A->egalites[j-1][k] = A->egalites[j][k];
        }
        A -> valeur_egalites[j-1] = A -> valeur_egalites[j];
    }
    A->nb_egalites --;
}





void delete_non_egalite(int no_ineg, attaque *A)
{
	//supprime la non egalite numero no_ineg de l attaque A
    int j=0,k=0;
    for (j=no_ineg+1; j<A->nb_non_egalites; j++) {
        for (k=0; k<input->max_rel; k++) {
            A->non_egalites[j-1][k] = A->non_egalites[j][k];
        }
        A -> valeur_non_egalites[j-1] = A -> valeur_non_egalites[j];
    }
    A->nb_non_egalites --;
}





void echange(int i, int j, attaque *A)
{
	//echange les egalites i et j de l attaque A
    int k=0;
    u64 temp=0;

    for (k=0; k<input->max_rel; k++) {
        temp = A->egalites[j][k];
        A->egalites[j][k] = A->egalites[i][k];
        A->egalites[i][k] = temp;
    }

    temp = (A -> valeur_egalites)[j];
    (A -> valeur_egalites)[j] = (A -> valeur_egalites)[i];
    (A -> valeur_egalites)[i] = temp;
}





void additionne_egalite(int i, int j, attaque *A)
{
	//additionne l egalite i a l egalite j de l attaque A
    int k=0;

    for (k=0; k<input->max_rel; k++) {
        A->egalites[j][k] ^= A->egalites[i][k];
    }
    (A -> valeur_egalites)[j] ^= (A -> valeur_egalites)[i];
}





void additionne_egalite_non_egalite(int i, int j, attaque *A)
{
	//additionne l egalite i a la non-egalite j de l attaque A
    int k=0;
    for (k=0 ; k<input->max_rel ; k++) {
        A->non_egalites[j][k] ^= A->egalites[i][k];
    }
    (A -> valeur_non_egalites)[j] ^= (A -> valeur_egalites)[i];
}





void additionne_non_egalite_non_egalite(int i, int j, attaque *A)
{
    //additionne l egalite i a la non-egalite j de l attaque A
    int k=0;
    for (k=0; k<input->max_rel; k++)
    {
        A->non_egalites[j][k] ^= A->non_egalites[i][k];
    }
    (A -> valeur_non_egalites)[j] ^= (A -> valeur_non_egalites)[i];
}





int no_variable(u64* liste)
{
	//renvoie TRUE ou FALSE suivant que la liste ne contient que 0UL ou non
    int k=0;
    for (k=0 ; k<input->max_rel ; k++) {
        if (liste[k] != 0ULL) return(FALSE);
    }
    return(TRUE);
}





int solo_variable(u64* liste, int no_var)
{
	// renvoie true si l equation contient uniquement la variable no_var
    // cet procedure n est utilisee que si on sait que l equation contient
    // deja la variable
    int res=0;
    enleve_var(liste, no_var);
    res = no_variable(liste);
    ajoute_var(liste,no_var);
    return(res);
}





int identique_liste_var(u64* liste1, u64* liste2)
{
	// renvoie TRUE si les deux listes sont identiques (meme combinaison lineaire de variable)
    int k=0;
    for (k=0; k<input->max_rel; k++) {
        if (liste1[k] != liste2[k]) return(FALSE);
    }
    return(TRUE);
}





int simplifie_systeme( attaque *A)
{
	int novar=0, i=0, j=0, m=0;
    j=0; // nombre d equations deja triangularisees (numero pivot)

    for (novar=0; novar < NbVar; novar++)
	{
        i=j;
        while (i < A->nb_egalites && !present(novar, A->egalites[i]) )
		{
            // on cherche la premiere equqtion contenant la variable novar
            i++;
        }

        if (i < A->nb_egalites) { // on a trouve une equation
            if (i>j) echange(i, j, A);

            // On supprime la variable novar presente dans l equation j dans toutes les autres equations et inequations
            for (i=j+1; i<A->nb_egalites; i++) {
                if (present(novar, A->egalites[i])) additionne_egalite(j, i, A);
            }
            for (i=0; i<A->nb_non_egalites; i++) {
                if (present(novar, A->non_egalites[i])) additionne_egalite_non_egalite(j,i,A);
            }
            j++;
        }
    }

    // On va supprimer les egalites et inegalites evidentes 0 = 0
    for (i=A->nb_egalites-1; i>=0; i--) {
        if (no_variable((A -> egalites)[i])) {
            if ((A-> valeur_egalites)[i] == 0UL) delete_egalite(i, A);
            else return(0);
        }
    }

    for (i=A->nb_non_egalites-1; i>=0; i--) {
        if (no_variable((A -> non_egalites)[i]) ) {
            if ((A -> valeur_non_egalites)[i] == 0UL) return(0);
            else delete_non_egalite(i, A);
        }
    }

    // On va supprimer les inegalites en doublon
    for (i=A->nb_non_egalites-1 ; i>=0 ; i--) {
        for (m=0 ; m<i ; m++) {
            if (identique_liste_var(A->non_egalites[m], A->non_egalites[i])
                 && (A->valeur_non_egalites[m] == A->valeur_non_egalites[i])) {
                delete_non_egalite(i, A);
                break;
            }
        }
    }
    return(1);
}





int genere_non_egalites2(attaque* A, int point1, int point2, int IouS, int no_sys, attaque **res)
{
    u64* comb = NULL;
	attaque* B = NULL;
    int j=0, cas=0, no_var1=0, no_var2=0;

    if (point2 == (input->phi + 1))
	{
        point1 ++;

        if (point1 == input->phi)
		{
            IouS ++;

            if (IouS == 2 || ISInd == TRUE)
			{
                if (no_sys == input->max_sys)
				{
                    printf("depassement %d cas\n",no_sys);
                    perror("Pas assez de places pour tous les cas");
                    exit(-1);
                }

                copie_attaque(A, *(res+no_sys));

                no_sys ++;
				if(B != NULL)
				{
					attaque_free(B);
					B = NULL;
				}

                return(no_sys);
            }
            point1 = 1;
        }
        point2 = point1 + 1;
    }

	B = attaque_init();
	comb = calloc(input->max_rel, sizeof(u64));

    for (cas = 0; cas < input->k; cas++)//Parti pris de choisir la premiere differente et les autres qqconques puis
	{

        // la premiere egale et la deuxieme differente etc
        copie_attaque(A, B);

        no_var1 = NoVar( point1, IouS); no_var2 = NoVar( point2, IouS);
        for (j = 0; j <= cas; j++)
		{

			init_comb(comb);
			ajoute_var(comb, no_var1);ajoute_var(comb, no_var2);

			if  (j < cas)// cas d egalite
			{
				ajoute_egalite(comb, 0ULL, B);
			}

			else// cas de l inegalite
			{
				ajoute_inegalite(comb, 0ULL, B);
			}

			no_var1 += input->phi; no_var2 += input->phi;
        }

        if (simplifie_systeme(B) != 0)
		{
            //DL
            if (!EXACT) {
                if (B->nb_egalites > barre_egalites) {
                    break; // On sort de la boucle de cas
                }
            }
            //DL
			free(comb);
			comb = NULL;
	        no_sys = genere_non_egalites2(B , point1, point2 + 1, IouS, no_sys, res);
			comb = calloc(input->max_rel, sizeof(u64));
        }
    }

	if(comb != NULL)
	{
		free(comb);
	}

	if(B != NULL)
	{
		attaque_free(B);
		B = NULL;
	}

    return( no_sys );
}





int genere_non_egalites(attaque* A, int IouS, attaque **res)
{

    // Pour tenir compte du fait que les entrees sont distinctes et les sorties aussi
    return(genere_non_egalites2(A, 1, 2, IouS, 0, res));
}





int nb_inegalites(int novar, attaque *A)
{
	// compte le nombre de non egalites contenant la variable novar
    int i=0,c=0;
    for (i=0; i< A->nb_non_egalites; i++) {
        if (present(novar, (A -> non_egalites)[i])) c++;
    }
    return(c);
}





int nombre_distincts(unsigned int* liste, int max_liste)
{
	//Compte le nombre d elements distincts dans la liste
    int res=max_liste, i=0, j=0;

    for (i=0; i<max_liste-1; i++)
	{
        for (j=i+1; j<max_liste; j++)
		{
            if (liste[i] == liste[j]) {
                res--;
                break;
            }
        }
    }
    return(res);
}





void init_ph(int*** ph)
{
	//initialise les facteurs polynomiaux pour le calcul de la variance
    int i=0,h=0;
	int** p1 = init_matrix(input->max_deg, input->max_deg);
	int** p2 = init_matrix(input->max_deg, input->max_deg);

    for (h=0; h<=input->phi ; h++) {
#ifdef DEBUG
        printf("h = %d\n",h);
#endif
        init_polynome2var(p1);

        p1[0][0] = 1;
        init_polynome2var(p2);
        p2[input->k][0] = 1;
        for (i=1; i<=h; i++) {
            p2[0][0] = -2*input->phi + i;
            produit_polynome2var(p1,p2,p1);
            produit_polynome2var(p1,p2,p1);
        }
        p2[input->k][0] = 0;
        p2[0][1] = 1;
        for (i=h+1; i<=input->phi ; i++) {
            p2[0][0] = -2*input->phi + i;
            produit_polynome2var(p1,p2,p1);
        }

        copie_polynome2var(p1,ph[h]);
#ifdef DEBUG
        affiche_polynome2var(ph[h]);
        printf("\n");
#endif
    }
	free_matrix(p1, input->max_deg);
	free_matrix(p2, input->max_deg);

	return;
}





int Point_present(int noPoint, unsigned int masque)
{
    if ((masque & (1UL <<(noPoint -1))) == 0UL) return(FALSE); else return(TRUE);
}





unsigned int ajoute_point(int noPoint, unsigned int masque)
{
    return( masque ^ (1UL <<(noPoint-1)));
}





int hamming(unsigned int u)
{
	//calcule le poids de hamming de u seulement sur les PHI premiers bits
    unsigned int a=1UL;
    int h=0, i=0;
    for (i=0 ; i<input->phi ; i++)
	{
        if ((u & a ) != 0UL) {
            h++;
        }
        a = a << 1;
    }
    return(h);
}





int trouve_Point(int i, unsigned int Ipris)
{
	//trouve le i eme point qui n est pas dans la liste indiquee par Ipris (bit =1 => Pris)
    int noPoint = 1;
    while (Point_present(noPoint, Ipris) || (i>0))
	{
        if (!Point_present(noPoint, Ipris)) i--;
        noPoint ++;
        if (noPoint > input->phi) {
            perror("Numero de point trop eleve");
            exit(-1);
        }
    }
    return(noPoint);
}





void genere_systemes2(int couple, int h, unsigned int Ipris, unsigned int Ippris, attaque** mix_sys)
{
	polynome r;
    int noPoint=0, noPointp=0, k=0, k1=0, l=0, l1=0;
    u64* comb = NULL;

    if (couple ==0)
	{
        compteur_general++;
		count(*mix_sys, 4 * input->k * input->phi, &r);
        somme_polynome(variance[h], r , variance+h);
        compteur[h]++;

#ifdef DEBUG
        if (h ==2)
		{ 
            printf("cas numero %d, h = %d \n\n", compteur[h],h);
            affiche_attaque(**mix_sys);
            printf("\n");
            affiche_polynome(&r);printf(".\n\n");
        }
#endif
    }

    else
	{
        // On va prendre le i+1 eme point libre et le j+1 eme point prime libre ou on a transforme couple = (i,j)
        noPoint = trouve_Point((couple - 1) / (input->phi - h), Ipris);
        noPointp =  trouve_Point((couple - 1) % (input->phi - h), Ippris);
		comb = calloc(input->max_rel, sizeof(u64));

        // On fabrique tous les cas possibles pour I(noPoint) <> I'(noPointp) et S(noPoint) <> S'(noPointp)
        for (k=0 ; k<input->k ; k++) {
            for (l=0 ; l<input->k ; l++) {
                copie_attaque(*mix_sys, *(mix_sys+1));

                for (k1=0; k1<k ; k1++) {
                    init_comb(comb);
                    ajoute_var(comb, NoVar(noPoint,0) + k1 * input->phi);
                    ajoute_var(comb, NoVarp(noPointp,0) + k1 * input->phi);
                    ajoute_egalite(comb, 0UL, *(mix_sys+1));// On ajoute les egalite I_k1(noPoint)=I'_k1(noPointp)
                }

                init_comb(comb);
                ajoute_var(comb,NoVar(noPoint, 0) + k * input->phi);
                ajoute_var(comb,NoVarp(noPointp, 0) + k * input->phi);
                ajoute_inegalite(comb,0UL, *(mix_sys+1)); //On rajoute l inegalite I_k(noPoint)<> I'_k(

				for (l1=0; l1<l  ; l1++) {
                    init_comb(comb);
                    ajoute_var(comb,NoVar(noPoint,1) + l1 * input->phi);
                    ajoute_var(comb,NoVarp(noPointp,1) + l1 * input->phi);
                    ajoute_egalite(comb,0UL, *(mix_sys+1)); // On ajoute les egalite S_l1(noPoint)=S'_l1(noPointp)
                }

                init_comb(comb);
                ajoute_var(comb,NoVar(noPoint, 1) + l * input->phi);
                ajoute_var(comb,NoVarp(noPointp, 1) + l * input->phi);
                ajoute_inegalite(comb,0UL, *(mix_sys+1)); //On rajoute l inegalite S_l(noPoint)<> S'_l(noPointp)

                if (simplifie_systeme(*(mix_sys+1))==1) {
					free(comb);
					comb = NULL;
                    genere_systemes2(couple-1, h, Ipris, Ippris, mix_sys+1);
					comb = calloc(input->max_rel, sizeof(u64));
                }
            }
        }
    }

	if(comb != NULL)
	{
		free(comb);
		comb = NULL;
	}

	return;
}





void genere_systemes(int noPoint, unsigned int Ipris, unsigned int  Ippris, attaque** mix_sys, attaque** mix_sys2)
{
	/* noPoint indique le numero du point considere, de 1 a phi inclus, par exemple noPoint=1 concerne I(1) et S(1) donc I_1(1), I_2(1) etc
      noPointp: entre 0 et phi inclus,  phi+1 possibilites au maximum: 0 -> aucune egalite, sinon egalite entre I(noPoint) et I'(noPointp)
     Ipris donne les numero des variables deja prises pour une egalite du cote I
     Ippris fait la meme chose pour I'
     */
    int noPointp=0, i=0, IouS=0, h=0;
    u64* comb = NULL;

    if (noPoint == input->phi + 1)
	{
        h = hamming(Ipris);
        copie_attaque(*(mix_sys + input->phi), *mix_sys2);
        genere_systemes2( (input->phi-h) * (input->phi-h), h, Ipris, Ippris, mix_sys2);

        return;
    }

    for (noPointp = 0; noPointp <= input->phi; noPointp++)
	{
        copie_attaque(*(mix_sys+noPoint-1), *(mix_sys+noPoint));
        if (noPointp == 0) genere_systemes(noPoint+1, Ipris, Ippris, mix_sys,mix_sys2);

        else
		{
            if (!Point_present(noPointp, Ippris))
			{
				comb = calloc(input->max_rel, sizeof(u64));
                for (i=0; i<input->k; i++)
				{
                    for (IouS=0; IouS <2; IouS++) {
                        init_comb(comb);
                        ajoute_var(comb,NoVar(noPoint,IouS) + i * input->phi);
                        ajoute_var(comb,NoVarp(noPointp,IouS) + i * input->phi);
                        ajoute_egalite(comb,0UL, *(mix_sys+noPoint));
                    }
                }

				free(comb);
				comb = NULL;

                if (simplifie_systeme(*(mix_sys+noPoint)) == 1) {
                    genere_systemes(noPoint+1, ajoute_point(noPoint, Ipris), ajoute_point(noPointp, Ippris), mix_sys, mix_sys2);
                }
            }
        }
    }
}





