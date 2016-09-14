#include "include/polynome.h"
#include "include/input.h"
#include "include/variance_feistel.h"
#include "include/variance_random.h"
#include "include/print.h"



//Global variables
//Default scheme
char *schema[SIZE_SCHEMA]={"I1(1)+I1(2)=0","I2(1)+I2(2)+S3(1)+S3(2)=0","S2(1)+S2(2)=0"};

int debug=40;
int* indice = NULL;
int compteur_general=0;
char nomvar[9]=".ijklmnop";
unsigned int* val_ineg = NULL;

// Input management
Input* input = NULL;


//Schema par default si on n'en donne pas.
char *schema[SIZE_SCHEMA];

//Management variables for feistel schema
u64**** mat_feistel;
unsigned long int** liste_partitions; 
int nb_partitions;
polynome p_somme,r_int;

// Result
polynome* variance;

// Others
unsigned int* basei;
int* compteur;
int NbVar; // Number of variables into the systeme 
int ISInd; // TRUE si dans le schema d attaque I et S sont independants (pas de combinaisons lineaire entre I et S 

int barre_egalites=30000; 
int barre_egalites_var=30000;




int main(int argc, char* argv[])
{
	// Input management
	input = calloc(1, sizeof(Input));
	int compute_variance = FALSE;
	int compute_feistel = FALSE;
	look_for_option(argc, argv, &compute_variance, &compute_feistel);


	// Initialization
	indice = init_array(input->max_non_egalites);
	u64* comb = calloc(input->max_rel, sizeof(u64));
    attaque* A = attaque_init(); 
    attaque* A2 = attaque_init();
	attaque** res = malloc(input->max_sys*sizeof(attaque*));
	attaque** res2 = malloc(input->max_sys*sizeof(attaque*));
	basei = calloc(input->max_non_egalites, sizeof(unsigned int));
	val_ineg = calloc(input->max_non_egalites, sizeof(unsigned int));

    int nb_sys=0, i=0, j=0, novar=0;
    ISInd = TRUE;
    polynome p,r ;

	for(i=0 ; i < input->max_sys ; i++)
	{
		res[i] = attaque_init();
		res2[i] = attaque_init();
	}

    //Polynome in formula
	int*** ph = malloc((input->phi+1) * sizeof(int**));
    int*** p_res = malloc((input->phi+1) * sizeof(int**));
	for(i=0 ; i < (input->phi+1) ; i++)
	{
		ph[i] = init_matrix(input->max_deg, input->max_deg);
		p_res[i] = init_matrix(input->max_deg, input->max_deg);
	}
 
	int** s = init_matrix(input->max_deg, input->max_deg);
	int k_phi = input->k * input->phi;
    init_polynome();

    // From scheme into list of equalities and inequalities
	int is_equality = TRUE;
    for (i=0 ; i < input->size_schema ; i++)
	{
        Str2CL(input->schema[i], comb, &is_equality);
		if(is_equality == TRUE)
		{
	        ajoute_egalite(comb, 0UL, A);
		}
		else
		{
			ajoute_inegalite(comb, 0UL, A);
		}
    }

	free(comb);
	comb = NULL;
    //DL
    if (ISInd) {
        if (!EXACT) {
            perror("not yet implemented");
            exit(-1);
        }
    }

	ISInd = FALSE; 
         
    NbVar = 4 * k_phi; // K*PHI for I, K*PHI for S, K*PHI for I', K*PHI for S'
    if (((NbVar-1) / input->proc) > input->max_rel) {
        perror("Not enough space for variables");
        exit(1);
    }

    //DL
    if (!EXACT)
	{
        barre_egalites = A->nb_egalites + NBTERMES;
        barre_egalites_var = 2*A->nb_egalites + NBTERMES;
    }

    // Mirror systems generation
    nb_sys = genere_non_egalites(A, 0, res);
	attaque_free(A);

#ifdef DEBUG
	printf("%d systemes a compter\n",nb_sys);
#endif

	attaque_free(A2);
    p.degre = 0;
    p.coef[0]=0;

	// P(N) computation
	for (i=0; i<nb_sys; i++)
    {
#ifdef DEBUG
        if (i<input->max_sysaff)
		{
            printf("Mirror system %d :\n\n",i+1);
            affiche_attaque(*res[i]);
        }
#endif
 
       count(*(res+i),2*input->k*input->phi, &r);

#ifdef DEBUG
        if (i<input->max_sysaff)
		{
              affiche_polynome(&r);printf("\n");
        }
#endif
        somme_polynome(p, r ,&p);        
    }
    
    
	printf("P(N) : \n");
    affiche_polynome(&p);
	puts("\n");


	// Computation of Q(N)
	if(compute_feistel == TRUE)
	{
		manage_feistel_schema();
		puts("\n");
	}

	// Computation of P(N,m)
	if(compute_variance == TRUE)
	{
		// calcul des coefficients multiplicateurs (polynomes en N et m) une fois pour toute.
		init_ph(ph);
		
		// On va fabriquer les systemes pour I' et S', il s'agit de decaler les numeros des variables des systemes trouves de 2*K*PHI
		// systeme primes pour le calcul de covariance
		attaque** resp = malloc(input->max_sys * sizeof(attaque*));
		for(i=0 ; i < input->max_sys ; i++)
		{
			resp[i] = attaque_init();
		}

		for (i=0; i<nb_sys; i++) {
		    copie_attaque(*(res+i) , *(resp+i) );

		    for (j=0; j<resp[i]->nb_egalites; j++) {
		        for (novar=0; novar < 2 * k_phi ; novar++) {

		            if (present(novar, resp[i]->egalites[j])) {
		                enleve_var(resp[i]->egalites[j], novar);
		                ajoute_var(resp[i]->egalites[j], novar + 2 * input->k * input->phi);
		            }
		        }
		    }

		    for (j=0; j<resp[i]->nb_non_egalites; j++) {
		        for (novar=0; novar < 2 * k_phi; novar++) {
		            if (present(novar, resp[i]->non_egalites[j])) {
		                enleve_var(resp[i]->non_egalites[j], novar);
		                ajoute_var(resp[i]->non_egalites[j], novar + 2 * input->k * input->phi);
		            }
		        }
		    }
		}

		// On fabrique maintenant les systemes combines des precedents systemes
		int ii=0;
		attaque** mix_sys = malloc( (input->phi+1) * sizeof(attaque*) );
		attaque** mix_sys2 = malloc( (input->phi * input->phi +1) * sizeof(attaque*) );
		compteur = init_array(input->phi + 1);
		variance = calloc(input->phi + 1, sizeof(polynome));

		for(i=0 ; i < input->phi+1 ; i++)
		{
			mix_sys[i] = attaque_init();
		}
		for(i=0 ; i < input->phi*input->phi +1 ; i++)
		{
			mix_sys2[i] = attaque_init();
		}
	
		for (i=0 ; i<=input->phi ; i++) { 
			//pour toutes les valeurs de h possibles
		    // on initialise les polynomes servant au calcul de variance a 0
		    // on multipliera a la fin par les polynomes constants ph
		    //On a un compteur pour le nombre de systemes analyses dans chaque cas
		    compteur[i] = 0;
		    variance[i].degre = 0;
		    variance[i].coef[0] = 0;
		}

		for (i=0; i<nb_sys; i++) {// i pour le premier systeme
		    for (j=0; j<nb_sys; j++) { // j pour le second
#ifdef DEBUG
		        printf("Systeme numero %d a analyser\n\n",i*nb_sys+j+1); //INFO
#endif

		        copie_attaque(*(res+i),*mix_sys);

		        for (ii=0; ii<resp[j]->nb_egalites; ii++) {
		            ajoute_egalite(resp[j]->egalites[ii], resp[j]->valeur_egalites[ii], *mix_sys);
		        }

		        for (ii=0; ii<resp[j]->nb_non_egalites; ii++) {
		            ajoute_inegalite(resp[j]->non_egalites[ii],resp[j]->valeur_non_egalites[ii], *mix_sys);
		        }
#ifdef DEBUG
		        affiche_attaque(**mix_sys);
#endif
		        genere_systemes(1,0UL,0UL, mix_sys, mix_sys2);
		    }
		}
//#ifdef DEBUG
		printf("(h, nb_sys) :\n");
//#endif
		for (i=0; i<=input->phi; i++)
		{
//#ifdef DEBUG
			printf("(%d,%d)\n",i,compteur[i]);
//#endif
		    polynome2polynome2var( variance[i], s);
		    produit_polynome2var(s, ph[i], p_res[i]);
//#ifdef DEBUG
		    affiche_polynome2var(s );
			puts("\n");
			printf(" x \n");affiche_polynome2var(ph[i]);printf("\n");
		    printf("=\n");
		    affiche_polynome2var(p_res[i]);printf("\n");
//#endif
		}

		printf("P(N,m) =\n");

		copie_polynome2var(p_res[0],s);

		for (i=1; i<=input->phi; i++)
		{
		    somme_polynome2var(p_res[i],s,s);
		}

		affiche_polynome2var(s);printf("\n");


		// We free all ressources needed in P(N,m)
		free_array(compteur);
		free(variance);
		for(i=0 ; i < input->max_sys ; i++)
		{
			attaque_free(resp[i]);
		}
		free(resp);

		for(i=0 ; i < input->phi+1 ; i++)
		{
			attaque_free(mix_sys[i]);
		}
		free(mix_sys);

		for(i=0 ; i < input->phi*input->phi +1 ; i++)
		{
			attaque_free(mix_sys2[i]);
		}
		free(mix_sys2);
	}


	// Now we free ressources
	for(i=0 ; i < input->phi+1 ; i++)
	{
		if(ph[i] != NULL)
		{
			free_matrix(ph[i], input->max_deg);
		}

		if(p_res[i] != NULL)
		{
			free_matrix(p_res[i], input->max_deg);
		}
	}

	if(ph != NULL)
	{
		free(ph);
	}
	
	if(p_res != NULL)
	{
		free(p_res);
	}

	free_matrix(s, input->max_deg);

	free_array(indice);

	free(basei);
	free(val_ineg);

	for(i=0 ; i < input->max_sys ; i++)
	{
		attaque_free(res[i]);
		attaque_free(res2[i]);
	}

	if(res != NULL)
		free(res);

	if(res2 != NULL)
		free(res2);

	// We have to free the input at the end only
	free_inputs(&input);

	return(0);
}



