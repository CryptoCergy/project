#ifndef DEF_INPUT
#define DEF_INPUT

#include "common.h"


// Tool functions used in input management
char* delete_spaces(char* string);
int parse_string(char* string, int default_value);
void manage_string(char* buffer);
void read_inputs(char* name_file);
void free_inputs(Input** input);

// This is the function user has to use to get all inputs (attack and configuration)
void look_for_option(int argc, char* argv[], int* compute_variance, int* compute_feistel);

#endif
