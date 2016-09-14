#include "include/input.h"





char* delete_spaces(char* string)
{
	while(string[0] == ' ')
	{
		string++;
	}

	return string;
}





int parse_string(char* string, int default_value)
{
	int result = 0;
	char* tmp = delete_spaces(string);
	result = atoi(tmp);

	if(result == -1)
	{
		result = default_value;
	}

	return result;
}





void manage_string(char* buffer)
{
	static int current_index = 0;
	int i=0, size_string = strlen(buffer);
	char* current_buffer = strdup(buffer);
	char* tmp = NULL;

	while( (i < size_string) && (current_buffer[i++] == ' ') );

	if(i < size_string)
	{
		tmp = current_buffer;
		strtok(tmp, "=");

		if(strncmp(tmp, "SIZE_SCHEMA", 11) == 0)
		{
			tmp = strtok(NULL, "=");
			input->size_schema = parse_string(tmp, SIZE_SCHEMA);
			input->schema = calloc(input->size_schema, sizeof(char*));
		}

		else if(strncmp(tmp, "K", 1) == 0)
		{
			tmp = strtok(NULL, "=");
			input->k = parse_string(tmp, K);
		}

		else if(strncmp(tmp, "PHI", 3) == 0)
		{
			tmp = strtok(NULL, "=");
			input->phi = parse_string(tmp, PHI);
		}

		else if (strncmp(tmp, "D", 1) == 0)
		{
			tmp = strtok(NULL, "=");
			input->d = parse_string(tmp, D);
		}

		else if (strncmp(tmp, "FEISTEL_TYPE", 12) == 0)
		{
			tmp = strtok(NULL, "=");
			input->feistel_type = parse_string(tmp, FEISTEL_TYPE);
		}

		else if (strncmp(tmp, "MAX_REL", 7) == 0)
		{
			tmp = strtok(NULL, "=");
			input->max_rel = parse_string(tmp, MAX_REL);
		}
		
		else if (strncmp(tmp, "PROC", 4) == 0)
		{
			tmp = strtok(NULL, "=");
			input->proc = parse_string(tmp, PROC);
		}
		
		else if (strncmp(tmp, "MAX_EGALITES", 12) == 0)
		{
			tmp = strtok(NULL, "=");
			input->max_egalites = parse_string(tmp, MAX_EGALITES);
		}
		
		else if (strncmp(tmp, "MAX_NON_EGALITES", 16) == 0)
		{
			tmp = strtok(NULL, "=");
			input->max_non_egalites = parse_string(tmp, MAX_NON_EGALITES);
		}
		
		else if (strncmp(tmp, "MAX_SYS", 7) == 0)
		{
			tmp = strtok(NULL, "=");
			input->max_sys = parse_string(tmp, MAX_SYS);
		}
		
		else if (strncmp(tmp, "MAXSYSAFF", 9) == 0)
		{
			tmp = strtok(NULL, "=");
			input->max_sysaff = parse_string(tmp, MAXSYSAFF);
		}
		
		else if (strncmp(tmp, "MAX_PART", 8) == 0)
		{
			tmp = strtok(NULL, "=");
			input->max_part = parse_string(tmp, MAX_PART);
		}
	
		else if(strncmp(tmp, "MAX_d", 5) == 0)
		{
			tmp = strtok(NULL, "=");
			input->max_d = parse_string(tmp, MAX_D);
		}
		
		else if (strncmp(tmp, "MAX_PHI", 7) == 0)
		{
			tmp = strtok(NULL, "=");
			input->max_phi = parse_string(tmp, MAX_PHI);
		}
		
		else if (strncmp(tmp, "MAX_DEG", 7) == 0)
		{
			tmp = strtok(NULL, "=");
			input->max_deg = parse_string(tmp, MAX_DEG);
		}

		else
		{
			input->schema[current_index++] = strndup(buffer, size_string - 1);
		}
	}

	free(current_buffer);
	return;
}





void read_inputs(char* name_file)
{
	int i=0;
	FILE* file = NULL;
	file = fopen(name_file, "r");
	char* buffer = malloc( (SIZE_INPUT+1) * sizeof(char) );
	char* back = calloc( (SIZE_INPUT+1), sizeof(char));

	if(file != NULL)
	{
		while(fgets(buffer, SIZE_INPUT, file) != NULL)
		{
			if(buffer[0] != ';')
			{		
				manage_string(buffer);
			}
			
			memcpy(buffer, back, SIZE_INPUT);
		}

		fclose(file);
	}
	
	else
	{
		char* default_schema[SIZE_SCHEMA]={"I1(1)+I1(2)=0","I2(1)+I2(2)+S3(1)+S3(2)=0","S2(1)+S2(2)=0"};
		input->feistel_type = FEISTEL_TYPE;
		input->d = D;
		input->k = K;
		input->phi = PHI;
		input->size_schema = SIZE_SCHEMA;
		input->schema = calloc(SIZE_SCHEMA, sizeof(char*));
		for(i=0 ; i < SIZE_SCHEMA ; i++)
		{
			input->schema[i] = strdup(default_schema[i]);
		}
	}

	memcpy(buffer, back, SIZE_INPUT);
	file = fopen("configuration.txt", "r");
	if(file!=NULL)
	{
		while(fgets(buffer, SIZE_INPUT, file) != NULL)
		{
			if(buffer[0] != ';')
			{		
				manage_string(buffer);
			}
			
			memcpy(buffer, back, SIZE_INPUT);
		}
				
		fclose(file);
	}

	else
	{
		input->max_rel = MAX_REL;		
		input->proc = PROC;
		input->max_egalites = MAX_EGALITES;
		input->max_non_egalites = MAX_NON_EGALITES;
		input->max_sys = MAX_SYS;
		input->max_sysaff = MAXSYSAFF;
		input->max_part = MAX_PART;
		input->max_d = MAX_D;
		input->max_phi = MAX_PHI;
		input->max_deg = MAX_DEG;
	}

	free(buffer);
	free(back);

	return;
}





void free_inputs(Input** input)
{
	int i;

	for(i=0 ; i < (*input)->size_schema ; i++)
	{
		free((*input)->schema[i]);
	}

	free((*input)->schema);
	free(*input);
	*input = NULL;

	return;
}




void look_for_option(int argc, char* argv[], int* compute_variance, int* compute_feistel)
{
	int result=0;
	
	if(argc > 4)
	{
		fprintf(stderr, "Usage : ./var [-v] [-f] input_file\n");
		exit(EXIT_FAILURE);
	}

	//If there is no option (set all values to default)
	if(argc == 1)
	{
		read_inputs(argv[2]);
	}

	//If there is one option (file or variance or feistel)
	else if(argc == 2)
	{
		if(strcmp(argv[1], "-v") == 0)
		{
			*compute_variance = TRUE;
			read_inputs(argv[2]);
		}

		else if(strcmp(argv[1], "-f") == 0)
		{
			*compute_feistel = TRUE;
			read_inputs(argv[2]);
		}

		else
		{
			read_inputs(argv[1]);
		}
	}

	//There are two options in (file or variance or feistel)
	else if(argc == 3)
	{
		if(strcmp(argv[1], "-v") == 0)
		{
			result=1;
			*compute_variance = TRUE;
		}

		else if(strcmp(argv[2], "-v") == 0)
		{
			result=2;
			*compute_variance = TRUE;
		}

		if(strcmp(argv[1], "-f") == 0)
		{
			result=1;
			*compute_feistel = TRUE;
		}

		else if(strcmp(argv[2], "-f") == 0)
		{
			result=2;
			*compute_feistel = TRUE;
		}

		if(result == 0)
		{
			fprintf(stderr, "Usage : ./var [-v] [-f] input_file\n");
			exit(EXIT_FAILURE);
		}

		if( (*compute_feistel == FALSE) || (*compute_variance == FALSE) )
		{
			int tmp = (result==1) ? 2 : 1;
			read_inputs(argv[tmp]);
		}

		else
		{
			read_inputs(argv[2]);
		}
	}

	//There are 3 options (file and variance and feistel)
	else
	{
		if(strcmp(argv[1], "-v") == 0)
		{
			*compute_variance = TRUE;
			if(strcmp(argv[2], "-f") == 0)
			{
				*compute_feistel = TRUE;
				read_inputs(argv[3]);
			}

			else if(strcmp(argv[3], "-f") == 0)
			{
				*compute_feistel = TRUE;
				read_inputs(argv[2]);
			}

			else
			{
				fprintf(stderr, "Usage : ./var [-v] [-f] input_file\n");
				exit(EXIT_FAILURE);
			}
		}

		else if(strcmp(argv[2], "-v") == 0)
		{
			*compute_variance = TRUE;
			if(strcmp(argv[1], "-f") == 0)
			{
				*compute_feistel = TRUE;
				read_inputs(argv[3]);
			}

			else if(strcmp(argv[3], "-f") == 0)
			{
				*compute_feistel = TRUE;
				read_inputs(argv[1]);
			}

			else
			{
				fprintf(stderr, "Usage : ./var [-v] [-f] input_file\n");
				exit(EXIT_FAILURE);
			}
			*compute_variance = TRUE;
		}

		else if(strcmp(argv[3], "-v") == 0)
		{
			*compute_variance = TRUE;
			
			if(strcmp(argv[1], "-f") == 0)
			{
				*compute_feistel = TRUE;
				read_inputs(argv[2]);
			}

			else if(strcmp(argv[2], "-f") == 0)
			{
				*compute_feistel = TRUE;
				read_inputs(argv[1]);
			}

			else
			{
				fprintf(stderr, "Usage : ./var [-v] [-f] input_file\n");
				exit(EXIT_FAILURE);
			}
		}
	}

	return;
}




