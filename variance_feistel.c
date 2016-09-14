#include "include/variance_feistel.h"


int nb_points;



int representant(u32 lp)
{
    //renvoie un representant du groupe lp
    u32 u;
    u32 no_point;
    u = 1UL << (nb_points-1);
    no_point = 0;

    while ((u != 0) && ((u & lp)==0))
	{
        u = u >> 1;
        no_point++;
    }

    if (u == 0)
	{
        printf("PROBLEME Impossible de trouver un representant\n");
        perror("BUG 2");
        exit(-1);
    }

    return(no_point);
}





void xor_var(u64* comb, int novar)
{
    comb[novar /input->proc] ^= (u64)(1ULL << (u64)(novar % input->proc));
}





void Str2CL_feistel(char *s, u64* comb, int* egalite)
{
    int novar, i, j, k ;
    for (i=0; i<input->max_rel; i++) {
        comb[i] = 0UL;
    }
    i=0;

    while (s[i] != '=' && s[i] != '!' && s[i] != '\0') {
        switch (s[i]) {
            case 'I':
                novar = 0;
                break;
            case 'S':
                novar = (input->k+input->d)*input->phi;
                break;
            default:
                perror("chaine non conforme");
                exit(-1);
                break;
        }
        i++;
        j = s[i]-'1';
        novar += j;
        if (j>9 || j<0) {
                perror("chaine non conforme");
                exit(-1);
        }
        i++;
        if (s[i] !='(') {
                perror("chaine non conforme");
                exit(-1);
        }
        i++;
        j = s[i]-'1';

        if (j>9 || j<0) {
            perror("chaine non conforme");
            exit(-1);
        }

        novar += (j)*(input->k+input->d);
        i++;
        if (s[i] !=')') {
            perror("chaine non conforme");
            exit(-1);
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

        if (novar < (input->k+input->d)*input->phi)
		{
			xor_var(comb, novar);
		}

        else
		{
            novar -= (input->k+input->d)*input->phi;
            j = novar / (input->k+input->d) ; //numero du point
            for (k=0; k<input->max_rel; k++) {
                comb[k] ^= mat_feistel[novar % (input->k+input->d)][input->d][j][k];
            }
        }
    }
    return;
}





int genere_non_egalites2_feistel(attaque A, int point1, int point2, int no_sys, attaque** res)
{
	attaque* B=NULL;
    int j=0, cas=0, no_var1=0, no_var2=0;
    u64* comb = calloc(input->max_rel, sizeof(u64));
    if (point2 == input->phi +1)
	{
        point1 ++;
        if (point1 == input->phi)
		{
            if (no_sys == input->max_sys)
			{
                printf("depassement %d cas\n",no_sys);
                perror("Pas assez de places pour tous les cas");
                exit(-1);
            }

            copie_attaque(&A, *(res+no_sys));
            no_sys++;
            free(comb);

            return(no_sys);
        }
        point2 = point1 + 1;
    }

	B = attaque_init();
    for (cas = 0; cas < input->k; cas++)
	{ 
        copie_attaque(&A, B);
        no_var1 = NoVar_feistel(point1); no_var2 = NoVar_feistel(point2);

         for (j = 0; j <= cas; j++) {
             init_comb(comb);
             ajoute_var(comb, no_var1); ajoute_var(comb, no_var2);

            if  (j < cas) {// cas d egalite
                ajoute_egalite(comb, 0UL, B);
            }
            else {// cas de l inegalite
                ajoute_inegalite(comb, 0UL, B);
            }
            no_var1 += 1; no_var2 += 1;
        }

        if (simplifie_systeme(B) != 0) {
//DL
            if (!EXACT) {
                if (B->nb_egalites > barre_egalites) {
                    break; // On sort de la boucle de cas
                }
            }
//DL
            no_sys = genere_non_egalites2_feistel(*B, point1, point2 + 1, no_sys, res);
        }
    }
    
    free(comb);
	attaque_free(B);

    return( no_sys );
}





int genere_non_egalites_feistel(attaque A, attaque** res)
{
	// Pour tenir compte du fait que les entrees sont distinctes et les sorties aussi
    return(genere_non_egalites2_feistel(A, 1, 2, 0, res));
}





int NoVar_feistel(int point)
{
    return( (point - 1)*(input->k+input->d));
}





void init_mat_feistel()
{
	int i=0,j=0,t=0,m=0,k=0;

    // Initialization with 0
    for (i=0; i<input->k ; i++) {
        for (j=0; j<input->d+1 ; j++) {
            for (t=0; t<input->phi; t++) {
                init_comb(mat_feistel[i][j][t]);
            }
        }
    }

    // Premier tour, I_(i,0)=I_i
    for (i=0; i<input->k ; i++) {
        for (t=0; t<input->phi; t++) {
			if(input->feistel_type != FEISTEL_TYPE3)
			{
	            ajoute_var(mat_feistel[i][0][t],i+(input->k+input->d)*t );
			}

			else
			{
	            ajoute_var(mat_feistel[i][0][t],((i+(input->k-2))%input->k)+(input->k+input->d)*t );
			}
        }
    }

    // tours suivants
    if (input->feistel_type==FEISTEL_TYPE1)
	{
        for (j=1; j<=input->d ; j++) {
            for (t=0; t<input->phi; t++) {
                // on fait glisser
                for (i=0; i<input->k-1; i++) {
                    for (m=0; m<input->max_rel; m++) {
                        mat_feistel[i][j][t][m] = mat_feistel[i+1][j-1][t][m];
                    }
                }
                for (m=0; m<input->max_rel; m++) {
                    mat_feistel[input->k-1][j][t][m] = mat_feistel[0][j-1][t][m];
                }

                // puis on rajoute la constante au premier
                ajoute_var(mat_feistel[0][j][t],input->k+t*(input->k+input->d)+j-1);
            }
        }
    }


    else if (input->feistel_type==FEISTEL_TYPEEXP) //FEISTEL EXPANSIF : On double le nombre de tour
	{
        for (t=0; t<input->phi ; t++)
		{
            for (j=0; j<input->d/(input->k -1); j++) // On groupe les tours K-1 par K-1
			{
                for (i=1; i<input->k; i++)
				{
                    for (k=0; k<input->k ; k++) //On recopie sur tous les sous-tours
					{
                        for (m=0; m<input->max_rel; m++)
						{
                            mat_feistel[k][i+j*(input->k-1)][t][m] = mat_feistel[k][i-1+j*(input->k-1)][t][m];
                        }
                    }
                    ajoute_var(mat_feistel[i][i+j*(input->k-1)][t],input->k+t*(input->k+input->d)+j*(input->k-1)+i-1);
                }
                
                // Puis on glisse en changeant la derniere ligne
                for (i=1; i<input->k ; i++) {
                    for (m=0; m<input->max_rel; m++) {
                        mat_feistel[i-1][(j+1)*(input->k-1)][t][m] = mat_feistel[i][(j+1)*(input->k-1)][t][m];
                    }
                }
                for (m=0; m<input->max_rel; m++) {
                    mat_feistel[input->k-1][(j+1)*(input->k-1)][t][m] = mat_feistel[0][j*(input->k-1)][t][m];
                }
            }
        }
    }


    else if (input->feistel_type==FEISTEL_TYPE2) { //type-2 en fait
        for (t=0; t<input->phi ; t++) {
            for (j=0; j<input->d/(input->k/2); j++) { // On groupe les tours K/2 par K/2
                for (i=1; i<input->k/2; i++) {
                    for (k=0; k<input->k ; k++) { //On recopie sur tous les sous-tours
                        // en decalant de 2 vers la gauche
                        for (m=0; m<input->max_rel; m++) {
                            mat_feistel[k][i+j*(input->k/2)][t][m] = mat_feistel[(k+2)%input->k ][i-1+j*(input->k/2)][t][m];
                        }
                    }
                    ajoute_var(mat_feistel[input->k-1][i+j*(input->k/2)][t],input->k+t*(input->k+input->d)+j*(input->k/2)+i-1);
                }
                // Mouvement different au dernier pas :
                // On glisse 3 fois vers la gauche et on rajoute une constante
                // au bon endroit
                for (k=0; k<input->k ; k++) {
                    for (m=0; m<input->max_rel; m++) {
                        mat_feistel[k][(j+1)*(input->k/2)][t][m] = mat_feistel[(k+3)% input->k ][(j+1)*(input->k/2)-1][t][m];
                    }
                }
                ajoute_var(mat_feistel[input->k-2][(j+1)*(input->k/2)][t],input->k+t*(input->k+input->d)+(j+1)*(input->k/2)-1);
            }
        }
    }

    if(input->feistel_type == FEISTEL_TYPE3)
    {
        for (j=1; j<=input->d ; j++) {
            for (t=0; t<input->phi; t++) {
                // on fait glisser
                for (i=1; i<=input->k-1; i++) {
                    for (m=0; m<input->max_rel; m++) {
                        mat_feistel[i][j][t][m] = mat_feistel[i-1][j-1][t][m];
                    }
                }
                for (m=0; m<input->max_rel; m++) {
                    mat_feistel[0][j][t][m] = mat_feistel[input->k-1][j-1][t][m];
                }
 
                // puis on rajoute la constante au premier
                ajoute_var(mat_feistel[2][j][t],input->k+t*(input->k+input->d)+j-1);
            }
        }
 
//      if((input->d/(input->k-1))%2==1) // si r est impair
//      {puts("########");
            for(j=0; j < input->phi ; j++)
            {
                u64 aa1 = mat_feistel[0][input->d][j][0];
                u64 aa2 = mat_feistel[0][input->d][j][1];
                u64 bb1 = mat_feistel[1][input->d][j][0];
                u64 bb2 = mat_feistel[1][input->d][j][1];
 
                for(t=0 ; t < input->max_rel ; t++){
                mat_feistel[0][input->d][j][t] = mat_feistel[2][input->d][j][t];
                mat_feistel[1][input->d][j][t] = mat_feistel[3][input->d][j][t];}
 
                mat_feistel[2][input->d][j][0] = aa1;
                mat_feistel[3][input->d][j][0] = bb1;
                mat_feistel[2][input->d][j][1] = aa2;
                mat_feistel[3][input->d][j][1] = bb2;
            }
//      }
    }

	else if (input->feistel_type==FEISTEL_TYPECON)
	{ //Contractant avec inversion droite gauche
        // I1    I2    I3    I4
        // I4+K1 I1    I2    I3 avec K1=f1(I1I2I3)
        // I3+K2 I4+K1 I1    I2 avec K2=f2((I4+K1)I1I2) etc
        // I2+K3 I3+K2 I4+K1 I1
        for (j=1; j<=input->d ; j++) {
            for (t=0; t<input->phi ; t++) {
                // on fait glisser
                for (i=1; i<input->k ; i++) {
                    for (m=0; m<input->max_rel ; m++) {
                        mat_feistel[i][j][t][m] = mat_feistel[i-1][j-1][t][m];
                    }
                }
                for (m=0; m<input->max_rel ; m++) {
                    mat_feistel[0][j][t][m] = mat_feistel[input->k-1][j-1][t][m];
                }
                // puis on rajoute la constante au premier
                ajoute_var(mat_feistel[0][j][t],input->k+t*(input->k+input->d)+j-1);
            }
        }
    }

	else if(input->feistel_type == FEISTEL_TYPEEGFN)
	{ //Feistel de Berger a 8 branches
		int round=0;


		for (t=0; t<input->phi ; t++)
		{
            for (j=0; j<input->d/(input->k/2) ; j++) // On groupe les tours K/2 par K/2
			{
				round = 0;
                for (i=1; i<=input->k/2; i++)
				{
                    for (k=0; k<input->k ; k++) //On recopie sur tous les sous-tours
					{
                        for (m=0; m<input->max_rel; m++)
						{
							mat_feistel[k][i+j*(input->k/2)][t][m] = mat_feistel[k][i-1+j*(input->k/2)][t][m];
//                            mat_feistel[k][i][t][m] = mat_feistel[k][i-1][t][m];
                        }
                    }

                	ajoute_var(mat_feistel[(i-1+input->k/2)%input->k][i+j*(input->k/2)][t],input->k+t*(input->k+input->d)+j*(input->k/2)+i-1);
					if(round > 0)
					{
						ajoute_toutes_var(mat_feistel[(i-1+input->k/2)%input->k][i+j*(input->k/2)][t], mat_feistel[input->k/2-1][i+j*input->k/2][t]);
//						ajoute_var(mat_feistel[(i-1+input->k/2)%input->k][i+j*(input->k/2)][t], input->k/2-1);
					}

					if(round == input->k/2-1)
					{
						for(k=input->k/2-2 ; k >0 ; k--){
//							ajoute_var(mat_feistel[(i-1+input->k/2)%input->k][i+j*(input->k/2)][t],k);
							ajoute_toutes_var(mat_feistel[(i-1+input->k/2)%input->k][i+j*(input->k/2)][t], mat_feistel[k][i+j*input->k/2][t]);
						}
					}

					round++;
                }



				//On inverse les coordonnees
/*
				if(j != 0)
				{
*/
					for(i=0 ; i < input->k/2 ; i+=2)
					{
						u64* value_first = calloc(input->max_rel, sizeof(u64));
						u64* value_second = calloc(input->max_rel, sizeof(u64));

						for(m=0 ; m < input->max_rel ; m++)
						{
/*							value_first[m] = mat_feistel[i][input->d][t][m];
							value_second[m] = mat_feistel[i+1][input->d][t][m];

							mat_feistel[i][input->d][t][m] = mat_feistel[i+input->k/2][input->d][t][m];
							mat_feistel[i+1][input->d][t][m] = mat_feistel[i+input->k/2+1][input->d][t][m];

							mat_feistel[i+input->k/2][input->d][t][m] = value_first[m];
							mat_feistel[i+input->k/2+1][input->d][t][m] = value_second[m];
*/
							value_first[m] = mat_feistel[i][(j+1)*input->k/2][t][m];
							value_second[m] = mat_feistel[i+1][(j+1)*input->k/2][t][m];

							mat_feistel[i][(j+1)*input->k/2][t][m] = mat_feistel[i+input->k/2][(j+1)*input->k/2][t][m];
							mat_feistel[i+1][(j+1)*input->k/2][t][m] = mat_feistel[i+input->k/2+1][(j+1)*input->k/2][t][m];

							mat_feistel[i+input->k/2][(j+1)*input->k/2][t][m] = value_first[m];
							mat_feistel[i+input->k/2+1][(j+1)*input->k/2][t][m] = value_second[m];
						}
						free(value_first);
						free(value_second);
					}
	//			}
            }

/*				for (i=1; i<input->k/2 ; i++) {
                    for (m=0; m<input->max_rel ; m++) {
                        mat_feistel[i+input->k/2-1][j][t][m] = mat_feistel[i-1][j-1][t][m];
                    }
                }
*/
/*
int iii=0;
while(iii < 25)
{
	affiche_var(iii++);puts("");
}printf("%d\n", input->k+t*(input->k+input->d)+j*(input->k-1)+i-1);exit(1);
*/
                
                // Puis on glisse en changeant la derniere ligne
                /*for (i=1; i<input->k ; i++)
				{
                    for (m=0; m<input->max_rel; m++) {
                        mat_feistel[i-1][(j+1)*(input->k-1)][t][m] = mat_feistel[i][(j+1)*(input->k-1)][t][m];
                    }
                }
                for (m=0; m<input->max_rel; m++) {
                    mat_feistel[input->k-1][(j+1)*(input->k-1)][t][m] = mat_feistel[0][j*(input->k-1)][t][m];
                }*/
        }
	}
}





int genere_partitions2(int no_part, int no_groupe, unsigned long int* part)
{
	int i=0;
    unsigned long int j=0;
    unsigned long int max=0,masque=0;

    if (no_groupe==input->phi) {
        for (i=0; i<input->phi; i++)
        {
            liste_partitions[no_part][i] = part[i];
        }
        
        return(no_part+1);
    }

    if (no_groupe==0)
    {
        max=(1UL <<input->phi)-1; // max atteint
        masque=0UL;
    }

    else
    {
        masque = part[0];
        for (i=1; i<no_groupe; i++) {
            masque |= part[i];
        }
        max = masque ^ ((1UL << input->phi)-1);
    }
    
    j = 1UL <<(input->phi-1);

    while ((j & masque)!=0UL && (j !=0UL))
    {
        j = j >> 1;
    }

    if (j==0UL)
    { // Plus de place
        for (i=no_groupe; i<input->phi; i++) {
            part[i] = 0UL;
        }
        return(genere_partitions2(no_part,input->phi,part));
    }

    while (j<=max) {
        if ((j & masque) == 0UL) {
            part[no_groupe]=j;
            no_part=genere_partitions2(no_part,no_groupe+1,part);
        }
        j++;
    }
    
    return(no_part);
}





int genere_partitions()
{
	unsigned long int* part = NULL;
	part = malloc(input->phi * sizeof(unsigned long int));

    // Va stocker tous les types de partitions
    int result = genere_partitions2(0, 0, part);

	if(part != NULL)
	{
		free(part);
	}

	return result;
}





int genere_systeme_feistel2(int no_sys, attaque A, int no_tour, attaque** res3)
{
	attaque *B = NULL, *C=NULL;
    u64 u=0, *comb=NULL;
    int i=0,j=0,k=0,j2=0,no_tour_suiv=0;
    int no_point1=0,no_point2=0;

    int nb_groupes,nb_couples,nb_cas,no_cas;
    int ii, jj,  cas;

    if (no_tour == input->d ) 
	{
    
#ifdef DEBUG
        if ((no_sys % input->max_sys)==0)
		{
            printf("depassement %d cas\n",no_sys);
        }
#endif

        if (no_sys < input->max_sys) {
            copie_attaque (&A, *(res3+no_sys));
        }

        no_sys ++;
#ifdef DEBUG
        printf("+");
#endif
        return(no_sys);
    }

	B = attaque_init();
	comb = calloc(input->max_rel, sizeof(u64));

    // on determine le tour suivant
    if (2*no_tour < input->d-1) no_tour_suiv = input->d-1-no_tour;
    else if (2*no_tour >input->d) no_tour_suiv = input->d-no_tour;
    else no_tour_suiv = input->d; 

    for (i=0; i<nb_partitions; i++)
	{
        copie_attaque(&A,B);
        for (j=0; j<nb_points; j++)
        {
            if (liste_partitions[i][j] != 0)
            {
                // On va rajouter les egalites entre les X (et les K)
                // du groupe j de la partition numero i
                u = 1UL << (nb_points-1);
                no_point1=0;

                while ((u != 0) && ((u & liste_partitions[i][j])==0)) {
                    u = u >> 1;
                    no_point1++;
                }

                if (u == 0) {
                    printf("Error in genere_feistel2\n");
                    exit(-1);
                }
                no_point2 = no_point1 + 1;
                u = u >> 1;

                while ((u != 0))
                {
                    if ((u & liste_partitions[i][j]) != 0) {

                        // On rajoute l egalite des X(no_point1) X(no_point2)
                        for (k=0; k<input->max_rel ; k++) {
                            comb[k] = mat_feistel[0][no_tour][no_point1][k] ^ mat_feistel[0][no_tour][no_point2][k];
                        }
                        ajoute_egalite(comb,0UL,B);

						if (input->feistel_type==FEISTEL_TYPECON) {
                            for (j2=1; j2<input->k-1; j2++) {
                                for (k=0; k<input->max_rel ; k++) {
                                    comb[k] = mat_feistel[j2][no_tour][no_point1][k] ^ mat_feistel[j2][no_tour][no_point2][k];
                                }
                                ajoute_egalite(comb,0UL,B);
                            }
                        }

                        init_comb(comb);
                        ajoute_var(comb,(input->k+input->d)*no_point1+input->k+no_tour);
                        ajoute_var(comb,(input->k+input->d)*no_point2+input->k+no_tour);
                        ajoute_egalite(comb,0UL,B);
                        B->nb_variables++;
                    }
                    u = u >> 1;
                    no_point2++;
                }

                // On va rajouter des inegalites cas non contractant
				if(input->feistel_type != FEISTEL_TYPECON)
				{
		            for (j2=j+1; j2<nb_points; j2++)
					{
		                if (liste_partitions[i][j2] != 0)
						{
							// Inegalite entre un representant du groupe j avec
                            // un representant du groupe j2
                            no_point2=representant(liste_partitions[i][j2]);

                            for (k=0; k<input->max_rel; k++)
							{
                                comb[k] = mat_feistel[0][no_tour][no_point1][k] ^ mat_feistel[0][no_tour][no_point2][k];
                            }
                            ajoute_inegalite(comb,0UL,B);
#ifdef BIJECTION
                            init_comb(comb);
                            ajoute_var(comb,(input->k+input->d)*no_point1+input->k+no_tour);
                            ajoute_var(comb,(input->k+input->d)*no_point2+input->k+no_tour);
                            ajoute_inegalite(comb,0UL,B);
#endif
						}
#ifdef BIJECTION
                        else
						{
                            init_comb(comb);
                            ajoute_var(comb,(input->k+input->d)*no_point1+input->k+no_tour);
                            ajoute_var(comb,(input->k+input->d)*input->phi+(input->phi-1)*no_tour+j2-1);
                            ajoute_inegalite(comb,0UL,B);
                    	}
#endif
/*
	                    // Inegalite entre un representant du groupe j avec
	                    // un representant du groupe j2
	                    u = 1UL << (input->phi-1);
	                    no_point2 = 0;

	                    while ((u != 0) && ((u & liste_partitions[i][j2])==0))
	                    {
	                        u = u >> 1;
	                        no_point2++;
	                    }

	                    if (u == 0)
	                    {
	                        printf("PROBLEME Numero point dans genere_feistel2\n");
	                        perror("BUG 2");
	                        exit(-1);
	                    }

	                    for (k=0; k<input->max_rel; k++)
	                    {
	                        comb[k] = mat_feistel[0][no_tour][no_point1][k] ^
	                        mat_feistel[0][no_tour][no_point2][k];
	                    }

	                    ajoute_inegalite(comb,0UL,B);
*/
	                }
	            } //FIN FEISTEL NON CONTRACTANT
			} // FIN IF LISTE_PARTITION
#ifdef BIJECTION //Revoir ce passage
	        else
			{
	            for (j2=j+1; j2<input->phi; j2++)
				{
	                init_comb(comb);
	                ajoute_var(comb,(input->k+input->d)*input->phi+(input->phi-1)*no_tour+j-1);
	                ajoute_var(comb,(input->k+input->d)*input->phi+(input->phi-1)*no_tour+j2-1);
	                ajoute_inegalite(comb,0UL,B);
	            }
	        }
#endif

        } // FIN for nb_points
	    /* A partir d ici on a tenu compte de la partition au tour no_tour dans le cas non contractant
	    Dans le cas contractant on n a tenu compte que des egalites dans les groupes
	     pas des inegalites entre les groupes.
	     Pour chaque couple de groupes il y a K-1 possibilites distinctes
	    */
	    for (k=0; k<nb_points; k++) {
	        B->chemin[no_tour][k] = liste_partitions[i][k];
	    }

		if (input->feistel_type == FEISTEL_TYPECON)
		{
	        nb_groupes=0;
	        for (j=0; j<nb_points; j++)
			{
	           if (liste_partitions[i][j] != 0) nb_groupes++;
	        }

	        if (nb_groupes==0) {
	            perror("CAS IMPOSSIBLE");
	            exit(-1);
	        }

	        if (nb_groupes==1)
			{
	            if (simplifie_systeme(B) !=0)
				{
	                no_sys = genere_systeme_feistel2(no_sys,*B,no_tour_suiv,res3);
	            }
	        }

	        else
			{ //calcul du nombre de cas : (K-1)^nb_couples
	            nb_couples=nb_groupes*(nb_groupes-1)/2;
	            nb_cas=1;
	            for (j=0; j<nb_couples; j++)
				{
	                nb_cas = nb_cas * (input->k-1); //RISQUE DE DEPASSEMENT ???
	            }

				C = attaque_init();

	            for (no_cas=0; no_cas<nb_cas; no_cas++)
				{//no_cas doit etre vu comme un nombre en base K-1
	                copie_attaque(B,C);
	                cas=no_cas;

	                for (ii=0; ii<nb_groupes-1; ii++)
					{
	                    no_point1=representant(liste_partitions[i][ii]);

	                    for (jj=ii+1; jj<nb_groupes; jj++)
						{
	                        no_point2=representant(liste_partitions[i][jj]);

	                        for (j2=0; j2<(cas%(input->k-1)); j2++)
							{
	                            for (k=0; k<input->max_rel; k++)
								{
	                                comb[k] = mat_feistel[j2][no_tour][no_point1][k] ^ mat_feistel[j2][no_tour][no_point2][k];
	                            }
	                            ajoute_egalite(comb,0UL,C);
	                        }

	                        for (k=0; k<input->max_rel; k++)
							{
	                            comb[k] = mat_feistel[cas%(input->k-1)][no_tour][no_point1][k] ^ mat_feistel[cas%(input->k-1)][no_tour][no_point2][k];
	                        }
	                        ajoute_inegalite(comb,0UL,C);
	                        cas = cas / (input->k-1);
	                    }
	                }
	                if (simplifie_systeme(C))
					{
	                    no_sys=genere_systeme_feistel2(no_sys,*C,no_tour_suiv,res3);
	                }
	            }

				attaque_free(C);
				C = NULL;
	        }
    	} // FIN FEISTEL CONTRACTANT

		else
		{
	        if (simplifie_systeme(B) !=0)
			{
	            no_sys = genere_systeme_feistel2(no_sys,*B,no_tour_suiv,res3);
	            //printf("Retour\n");
	        }
        }


	} //FIN NB_PARTITION

     /*   // A partir d ici on a tenu compte de la partition au tour no_tour
        if (simplifie_systeme(B) !=0)
        {
            for (k=0; k<input->phi; k++)
            {
                B->chemin[no_tour][k] = liste_partitions[i][k];
            }

            no_sys = genere_systeme_feistel2(no_sys,*B,no_tour_suiv,res3);
			attaque_free(B);
			B = NULL;
			B = attaque_init();
        }
    */

	if(comb != NULL)
	{
		free(comb);
	}

	attaque_free(B);

    return(no_sys);
}





int count_and_genere_systeme_feistel2(int no_sys, attaque A, int no_tour)
{
	attaque* B = attaque_init();
    u64* comb = NULL;
	u64 u=0;
    int i=0,j=0,k=0, j2=0, no_tour_suiv=0;
    int no_point1=0,no_point2=0;

    if (no_tour == input->d ) {
        copie_attaque(&A,B);
        count_feistel(B,B->nb_variables, &r_int);
        somme_polynome(p_somme, r_int ,&p_somme);
#ifdef DEBUG
        if (r_int.degre== p_somme.degre) { // on affiche que les cas les + significatifs
            printf("Path %d\n",no_sys+1);
            affiche_attaque_feistel(A);
            affiche_attaque2(A);
            affiche_polynome(&r_int);printf("\n");
            affiche_polynome(&p_somme);
            printf("\n");
        }
#endif
        no_sys++;
		attaque_free(B);
		
        return(no_sys);
    }

    // on determine le tour suivant
    if (2*no_tour < input->d-1) no_tour_suiv = input->d-1-no_tour;
    else if (2*no_tour >input->d) no_tour_suiv = input->d-no_tour;
    else no_tour_suiv = input->d; 

	comb = calloc(input->max_rel, sizeof(u64));

    for (i=0; i<nb_partitions; i++)
	{

        copie_attaque(&A,B);
        for (j=0; j<input->phi; j++) {
            if (liste_partitions[i][j] != 0) {
                // On va rajouter les egalites entre les X (et les K)
                // du groupe j de la partition numero i
                u = (u64)(1ULL << (input->phi-1));
                no_point1=0;
                while ((u != 0) && ((u & liste_partitions[i][j])==0)) {
                    u = (u64)(u >> 1);
                    no_point1++;
                }

                if (u == 0) 
				{
                    exit(-1);
                }

                no_point2 = no_point1 + 1;
                u = (u64)(u >> 1);

                while ((u != 0)) {
                    if ((u & liste_partitions[i][j]) != 0) {
                        // On rajoute l egalite des X(no_point1) X(no_point2)
                        for (k=0; k<input->max_rel; k++) {
                            comb[k] = mat_feistel[0][no_tour][no_point1][k] ^
                            mat_feistel[0][no_tour][no_point2][k];
                        }
                        ajoute_egalite(comb,0UL,B);

                        init_comb(comb);
                        ajoute_var(comb,(input->k+input->d)*no_point1+input->k+no_tour);
                        ajoute_var(comb,(input->k+input->d)*no_point2+input->k+no_tour);
                        ajoute_egalite(comb,0UL,B);
                        B->nb_variables++;
                    }
                    u = u >> 1;
                    no_point2++;
                }

                // On va rajouter des inegalites
                for (j2=j+1; j2<input->phi; j2++) {
                    if (liste_partitions[i][j2] != 0) {
                        // Inegalite entre un representant du groupe j avec
                        // un representant du groupe j2
                        u = 1UL << (input->phi-1);
                        no_point2 = 0;

                        while ((u != 0) && ((u & liste_partitions[i][j2])==0))
						{
                            u = u >> 1;
                            no_point2++;
                        }

                        if (u == 0)
						{
                            exit(-1);
                        }

                        for (k=0; k<input->max_rel; k++)
                        {
                            comb[k] = mat_feistel[0][no_tour][no_point1][k] ^ mat_feistel[0][no_tour][no_point2][k];
                        }
                        ajoute_inegalite(comb, 0UL, B);
#ifdef BIJECTION
                        init_comb(comb);
                        ajoute_var(comb,(input->k+input->d)*no_point1+input->k+no_tour);
                        ajoute_var(comb,(input->k+input->d)*no_point2+input->k+no_tour);
                        ajoute_inegalite(comb,0UL,B);
#endif
                    }
#ifdef BIJECTION
                    else {
                        init_comb(comb);
                        ajoute_var(comb,(input->k+input->d)*no_point1+input->k+no_tour);
                        ajoute_var(comb,(input->k+input->d)*input->phi+(input->phi-1)*no_tour+j2-1);
                        ajoute_inegalite(comb,0ULL,B);
                    }
#endif
                }
            }
#ifdef BIJECTION
            else {
                for (j2=j+1; j2<input->phi; j2++) {
                    init_comb(comb);
                    ajoute_var(comb,(input->k+input->d)*input->phi+(input->phi-1)*no_tour+j-1);
                    ajoute_var(comb,(input->k+input->d)*input->phi+(input->phi-1)*no_tour+j2-1);
                    ajoute_inegalite(comb,0ULL,B);
                }
            }
#endif
        }

        // A partir d ici on a tenu compte de la partition au tour no_tour
        if (simplifie_systeme(B) !=0) {
            for (k=0; k<input->phi; k++) {
                B->chemin[no_tour][k] = liste_partitions[i][k];
            }
            no_sys = count_and_genere_systeme_feistel2(no_sys, *B, no_tour_suiv);
        }
    }

	free(comb);
	attaque_free(B);
    
    return(no_sys);
}





int genere_systeme_feistel(int no_sys, attaque A, attaque** res3)
{
	// On va rajouter a l attaque A tous les cas possibles
    // suivant les egalites des parties gauches ou non
    // pour chacun des tours.
    // C est a dire au plus (B_PHI)^D possibilites.
    return(genere_systeme_feistel2(no_sys, A, 0, res3));
}





int count_and_genere_systeme_feistel(int nb_sys, attaque A)
{
	// On va rajouter a l attaque A tous les cas possibles
    // suivant les egalites des parties gauches ou non
    // pour chacun des tours.
    // C est a dire au plus (B_PHI)^D possibilites.
    return(count_and_genere_systeme_feistel2(nb_sys,A,0));
}





void manage_feistel_schema()
{
	int i=0,j=0,k=0;
	nb_points = input->phi;
	mat_feistel = malloc(input->k*sizeof(u64***));

	for(i=0 ; i < input->k ; i++)
	{
		mat_feistel[i] = malloc((input->d+1) * sizeof(u64**));
		for(j=0 ; j < input->d+1 ; j++)
		{
			mat_feistel[i][j] = malloc(input->phi * sizeof(u64*));
			for(k=0 ; k < input->phi ; k++)
			{
				mat_feistel[i][j][k] = calloc(input->max_rel, sizeof(u64));
			}
		}
	}

	liste_partitions = malloc(input->max_part*sizeof(unsigned long int*)); 
	for(i=0 ; i < input->max_part ; i++)
	{
		liste_partitions[i] = calloc(input->phi, sizeof(unsigned long int));
	}

    int nb_sys=0, nb_sys3=0;
    attaque* A = attaque_init(); 
	attaque** res = malloc(input->max_sys*sizeof(attaque*));
	attaque** res1 = malloc(input->max_sys*sizeof(attaque*));
	u64* comb = calloc(input->max_rel, sizeof(u64));
 
	for(i=0 ; i < input->max_sys ; i++)
	{
		res[i] = attaque_init();
		res1[i] = attaque_init();
	}

    // Nbvar initialization
#ifdef BIJECTION
	NbVar = (input->k+input->d) * input->phi + input->d * (input->phi -1);
#else
    NbVar = (input->k+input->d) * input->phi; 
#endif
    if (((NbVar-1) / input->proc) >= input->max_rel) {
        perror("PAS ASSEZ DE PLACE POUR LES VARIABLES");
        exit(1);
    }
    init_polynome();
    init_mat_feistel();
 
 
//    affiche_mat_feistel();
//exit(1);

    A->nb_egalites = 0;
    A->nb_non_egalites = 0;

    // Transformation du schema de l attaque en liste d egalites
	int equality = TRUE;
    for (i=0; i<input->size_schema; i++)
    {
        Str2CL_feistel(input->schema[i], comb, &equality);
		if(equality == TRUE)
		{
	        ajoute_egalite(comb, 0ULL, A);
		}

		else
		{
			ajoute_inegalite(comb, 0ULL, A);
		}
    }  
   
    for (i=0; i< A->nb_egalites; i++)
	{
        A->valeur_egalites[i] = 0;
    }

#ifdef DEBUG
    affiche_attaque_feistel(*A);
#endif
    // On genere les systemes miroir en tenant compte que les entrees  doivent etre distinctes
    nb_sys = genere_non_egalites_feistel(*A, res1); 

	attaque_free(A);

#ifdef DEBUG
    for (i=0; i<nb_sys; i++)
	{
        printf("Mirror system %d :\n\n",i+1);
        affiche_attaque_feistel(*res1[i]);
    }
#endif
  
    nb_partitions=genere_partitions();

#ifdef DEBUG
    printf("Il y a %d partitions\n",nb_partitions);
#endif

    p_somme.degre = 0;
    p_somme.coef[0]=0; 
    nb_sys3=0;
	for(i=0;i<nb_sys; i++)
	{
#ifdef DEBUG
        printf("Mirror root-system %d :\n\n",i+1);
        affiche_attaque(*res1[i]);
		puts("");
#endif
puts("test");
        res1[i]->nb_variables = (input->k+input->d)*input->phi;
        // pour Psi5 : (2+5)x2=14 (3 conditions -> 14-3=11)
        //Pour F3,6 (3+12)x4=60 var. 10 c. -> 60-10=50
// /*       nb_sys3=count_and_genere_systeme_feistel(nb_sys3,*res1[i]);   */
        nb_sys3=genere_systeme_feistel(nb_sys3,*res1[i],res);
    }
    printf("TOTAL (%d paths):\n",nb_sys3);

// /*
    for(i=0;i<nb_sys3; i++) {
//printf("i======%d\n", i);
#ifdef DEBUG
        if (i<input->max_sysaff) {
            printf("Mirror system %d :\n\n",i);
            affiche_attaque(*res[i]);
            affiche_attaque2(*res[i]);
        }
#endif
        count(res[i],res[i]->nb_variables, &r_int);
        somme_polynome(p_somme, r_int ,&p_somme);
#ifdef DEBUG
        if (i<input->max_sysaff){
            affiche_polynome(&r_int);printf("\n");
            affiche_polynome(&p_somme);
            printf("\n");
        }
#endif
    }
// */
    printf("Q(N)=");
    affiche_polynome(&p_somme);
    printf("\n");
/*
    for(i=0;i<nb_sys; i++) {
#ifdef DEBUG
        printf("Mirror system %d :\n\n",i+1);
        affiche_attaque_feistel(*res[i]);
#endif
        res[i]->nb_variables = (input->k+input->d)*input->phi;
        nb_sys3=count_and_genere_systeme_feistel(nb_sys3,*res[i]);
    }

    printf("Q(N)=\n");
    affiche_polynome(&p_somme);
    printf("\n");
*/
	nb_points = 2*input->phi;

	for(i=0 ; i < input->k ; i++)
	{
		for(j=0 ; j < input->d+1 ; j++)
		{
			for(k=0 ; k < input->phi ; k++)
			{
				if(mat_feistel[i][j][k] != NULL)
				{
					free(mat_feistel[i][j][k]);
				}
			}
			if(mat_feistel[i][j] != NULL)
			{
				free(mat_feistel[i][j]);
			}
		}
		if(mat_feistel[i] != NULL)
		{
			free(mat_feistel[i]);
		}
	}

	if(mat_feistel != NULL)
	{
		free(mat_feistel);
	}
 
	for(i=0 ; i < input->max_part ; i++)
	{
		if(liste_partitions[i] != NULL)
		{
			free(liste_partitions[i]);
		}
	}

	if(liste_partitions != NULL)
	{
		free(liste_partitions);
	}

	free(comb);

	for(i=0 ; i < input->max_sys ; i++)
	{
		attaque_free(res[i]);
		attaque_free(res1[i]);
	}

	free(res);
	free(res1);

	return;
}





