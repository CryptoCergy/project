#ifndef DEF_SMS4
#define DEF_SMS4

#include <stdlib.h>

/** UTILISATION **************************************
	//Nombre de tour
	int r = 32;
	int* key = calloc(4, sizeof(int));

	key[0] = 19088743;
	key[1] = 2309737967;
	key[2] = 4275878552;
	key[3] = 1985229328;

	unsigned int* rk = key_schedule_sms4(r, key);

	//La fonction alloue rk il faut liberer l'espace
	free(rk);
	free(key);
*****************************************************/

#define LR(a,b) ((a << b) | (a >> (32-b)))
#define ROT_L(a,b) ((unsigned int) (a << b) | (unsigned int) (a >> (8*sizeof(unsigned int)-b)))


unsigned int L(unsigned int value);
unsigned int S(int value);
void computeCK(unsigned char* CK, int i);
unsigned int ck2Int(unsigned char* CK);
unsigned int* key_schedule_sms4(int rounds, int* key);

#endif
