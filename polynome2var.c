#include "include/polynome2var.h"


char NOM_VAR1='N';
char NOM_VAR2='m';


u32* init_array32(int size)
{
	u32* result = NULL;
	result = malloc(size * sizeof(u32));
	if(result == NULL)
	{
		fprintf(stderr, "Error when allocating array\n");
		perror("");
		exit(EXIT_FAILURE);
	}

	return result;
}





void free_array32(u32* array)
{
	if(array != NULL)
	{
		free(array);
	}

	return;
}



u32** init_matrix32(int size1, int size2)
{
	int i=0;
	u32** result = NULL;

	result = (u32**)calloc(size1, sizeof(u32*));
	if(result == NULL)
	{
		fprintf(stderr, "Error when allocating matrix\n");
		perror("");
		exit(EXIT_FAILURE);
	}

	for(i=0 ; i < size1 ; i++)
	{
		result[i] = NULL;
		result[i] = (u32*)calloc(size2, sizeof(u32));
		if(result[i] == NULL)
		{
			fprintf(stderr, "Error when allocating matrix\n");
			perror("");
			exit(EXIT_FAILURE);
		}
	}

	return result;
}



void free_matrix32(u32** matrix, int size1)
{
	int i=0;
	for(i=0 ; i < size1 ; i++)
	{
		if(matrix[i] != NULL)
		{
			free(matrix[i]);
		}
	}

	if(matrix != NULL)
	{
		free(matrix);
	}

	return;
}


u64* init_array64(int size)
{
	u64* result = NULL;
	result = malloc(size * sizeof(u64));
	if(result == NULL)
	{
		fprintf(stderr, "Error when allocating array\n");
		perror("");
		exit(EXIT_FAILURE);
	}

	return result;
}





void free_array64(u64* array)
{
	if(array != NULL)
	{
		free(array);
	}

	return;
}



u64** init_matrix64(int size1, int size2)
{
	int i=0;
	u64** result = NULL;

	result = (u64**)calloc(size1, sizeof(u64*));
	if(result == NULL)
	{
		fprintf(stderr, "Error when allocating matrix\n");
		perror("");
		exit(EXIT_FAILURE);
	}

	for(i=0 ; i < size1 ; i++)
	{
		result[i] = NULL;
		result[i] = (u64*)calloc(size2, sizeof(u64));
		if(result[i] == NULL)
		{
			fprintf(stderr, "Error when allocating matrix\n");
			perror("");
			exit(EXIT_FAILURE);
		}
	}

	return result;
}



void free_matrix64(u64** matrix, int size1)
{
	int i=0;
	for(i=0 ; i < size1 ; i++)
	{
		if(matrix[i] != NULL)
		{
			free(matrix[i]);
		}
	}

	if(matrix != NULL)
	{
		free(matrix);
	}

	return;
}



int* init_array(int size)
{
	int* result = NULL;
	result = malloc(size * sizeof(int));
	if(result == NULL)
	{
		fprintf(stderr, "Error when allocating array\n");
		perror("");
		exit(EXIT_FAILURE);
	}

	return result;
}





void free_array(int* array)
{
	if(array != NULL)
	{
		free(array);
	}

	return;
}





int** init_matrix(int size1, int size2)
{
	int i=0;
	int** result = NULL;

	result = (int**)calloc(size1, sizeof(int*));
	if(result == NULL)
	{
		fprintf(stderr, "Error when allocating matrix\n");
		perror("");
		exit(EXIT_FAILURE);
	}

	for(i=0 ; i < size1 ; i++)
	{
		result[i] = NULL;
		result[i] = (int*)calloc(size2, sizeof(int));
		if(result[i] == NULL)
		{
			fprintf(stderr, "Error when allocating matrix\n");
			perror("");
			exit(EXIT_FAILURE);
		}
	}

	return result;
}





void free_matrix(int** matrix, int size1)
{
	int i;
	for(i=0 ; i < size1 ; i++)
	{
		if(matrix[i] != NULL)
		{
			free(matrix[i]);
		}
	}

	if(matrix != NULL)
	{
		free(matrix);
	}

	return;
}




//donne l oppose d un polynome
void moins2var(int** p) {
	int i=0,j=0;
	for(i=0; i<input->max_deg;i++) {
        for (j=0; j<input->max_deg; j++) {
            p[i][j] = -p[i][j];
        }
	}
	return;
}



void somme_polynome2var(int** p, int** q, int** r) {
	int i=0,j=0;
	for(i=0; i<input->max_deg;i++) {
        for (j=0; j<input->max_deg; j++) {
            r[i][j] = p[i][j] + q[i][j];
        }
	}
	return;
}



void init_polynome2var(int** p) {
    int i=0, j=0;

    for(i=0 ; i < input->max_deg ; i++)
	{
        for (j=0 ; j < input->max_deg ; j++)
		{
            p[i][j] = 0;
        }
	}    
}



void copie_polynome2var(int** p, int** q) {
	int i,j;

    for (i=0; i<input->max_deg; i++) {
        for (j=0; j<input->max_deg; j++) {
            q[i][j] = p[i][j];
        }
    }
	return;
}



void produit_polynome2var (int** p, int** q, int** r) {
	int i=0,j=0,k=0,l=0;
    int prod=0;

	int** s = init_matrix(input->max_deg, input->max_deg);
	init_polynome2var(s);

	for(i=0; i<input->max_deg;i++) {
        for (j=0; j<input->max_deg; j++) {
            for (k=0; k<input->max_deg; k++) {

				if((i+k) >= input->max_deg)
				{
					continue;
				}

                for (l=0; l<input->max_deg; l++) {

					if((j+l) >= input->max_deg)
					{
						continue;
					}

                    prod = p[i][j] * q[k][l];
                    if ( ((i+k >= input->max_deg) || (j+l >= input->max_deg)) && prod != 0 )
					{
                        perror("Augmenter MAX_DEG\n");
                        exit(-1);
                    }
                    s[i+k][j+l] += prod;
                }
            }
        }
	}

    copie_polynome2var(s,r);

	free_matrix(s, input->max_deg);

	return;
}



void affiche_pui(int i, char c) {
    // affiche X^i avec le charactere c pour X
    switch (i) {
        case 1:
/*#ifdef SCILAB
			if(c=='m')
			{
				printf("*%c", c);
			}
			else
			{
				printf("%c", c);
			}
#else
*/
            printf("%c",c);
//#endif
            break;
        case 0:
            break;
        default:
#ifdef LATEX2
            printf("%c^{%d}",c,i);
#elif defined(SCILAB)
/*			if(c=='m')
			{
				printf("*%c^(%d)",c,i);
			}
			else
			{*/
				printf("%c^(%d)",c,i);
//			}

#else
            printf("%c^%d",c,i);
#endif
            break;
    }   
}



void affiche_polynome2var(int** p) {
	int i,j;
    int flag=0; // 0 si on a rien affiche encore.
#ifdef LATEX
    printf("$");
#endif
    for (i=0; i<input->max_deg; i++) {
        for (j=0; j<input->max_deg; j++) {

            if (p[i][j] != 0) {
                if (flag == 1) {
					if (p[i][j]>0){
						printf("+");
					}
				}
				else flag = 1;

                if (i>0 || j>0) {
                    if (p[i][j] !=1) {
                        if (p[i][j]==-1) {
                            printf("-");
                        }
                        else{
#ifdef SCILAB
							if(i>0)
							{
								 printf("%d*",p[i][j]);
							}
							else
							{
								printf("%d", p[i][j]);
							}
#else
							 printf("%d",p[i][j]);
#endif
						}
                    }
                }
                else{
					 printf("%d",p[i][j]);
				}

                affiche_pui(i,NOM_VAR1);
#ifdef SCILAB
				if( (j>0) && ((i>0) || ((i==0) && (p[i][j]!=1) && (p[i][j] !=-1))) )
				{
					printf("*");
				}
#endif
                affiche_pui(j,NOM_VAR2);
            }
        }
    }
    if (flag ==0){
		 printf("0");
	}
#ifdef LATEX
    printf("$");
#endif
}


