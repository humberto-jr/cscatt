/******************************************************************************

 About
 -----

 This module is a collection of routines for file manipulation all based upon
 the C library.

******************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdbool.h>

#include "globals.h"
#include "file.h"

/******************************************************************************

 Function file_find_string(): scans a given input file searching for the first
 occurrence of a given pattern. If found it shall return the whole line or an
 empty string otherwise.

 NOTE: Lines starting by '#' are ignored.

******************************************************************************/

char *file_find_string(FILE *input, const char pattern[])
{
	char *line = malloc(sizeof(char)*MAX_LINE_LENGTH);

	ASSERT(line != NULL)
	ASSERT(input != NULL)

	rewind(input);
	while (fgets(line, MAX_LINE_LENGTH, input) != NULL)
	{
		if (line[0] != '#' && strstr(line, pattern) != NULL) return line;
	}

	free(line);
	return "\n";
}

/******************************************************************************

 Function file_get_key(): scans a given input file searching for the first
 occurrence of a keyword with the format "[key] = [value]". If found, and
 if [min] <= [value] <= [max], it shall return [value]. Otherwise, it
 returns a default value.

 NOTE: Lines starting by '#' are ignored.

******************************************************************************/

double file_get_key(FILE *input, const char key[],
                    const double min, const double max, const double default_value)
{
	ASSERT(max > min)

	char *line = file_find_string(input, key);
	char *token = strtok(line, "=");

	while (token != NULL)
	{
		if (strstr(line, key) != NULL)
		{
			token = strtok(NULL, "=");
			const double value = atof(token);

			if (line != NULL) free(line);
			if (value < min) return min;
			if (value > max) return max;

			return value;
		}

		token = strtok(NULL, "=");
	}

	return default_value;
}

/******************************************************************************

 Function file_row_count(): counts the actual number of valid lines in a given
 input data file.

 NOTE: Blank lines or lines starting by '#' are ignored.

******************************************************************************/

int file_row_count(FILE *input)
{
	ASSERT(input != NULL)

	int counter = 0;
	char line[MAX_LINE_LENGTH] = "\n";

	rewind(input);
	while (fgets(line, sizeof(line), input) != NULL)
	{
		if (line[0] != '#' && line[0] != '\n') ++counter;
	}

	return counter;
}

/******************************************************************************

 Function file_col_count(): counts the actual number of columns in a given
 input data file.

 NOTE: Blank lines or lines starting by '#' are ignored.

******************************************************************************/

int file_col_count(FILE *input)
{
	ASSERT(input != NULL)

	int counter = 0;
	char line[MAX_LINE_LENGTH] = "\n";

	rewind(input);
	while (fgets(line, sizeof(line), input) != NULL)
	{
		if (line[0] != '#' && line[0] != '\n') break;
	}

	char *token = strtok(line, " ");

	while (token != NULL)
	{
		++counter;
		token = strtok(NULL, " ");
	}

	return counter;
}

/******************************************************************************

 Function file_open_input(): open a file named filename in write-only mode.

 NOTE: binary format is used if bin_format = true.

******************************************************************************/

FILE *file_open_input(const char filename[], const bool bin_format)
{
	FILE *input
		= (bin_format? fopen(filename, "rb") : fopen(filename, "r"));

	if (input == NULL)
	{
		PRINT_ERROR("unable to open %s\n", filename)
		exit(EXIT_FAILURE);
	}

	return input;
}

/******************************************************************************

 Function file_open_input(): open a file named filename in read-only mode.

 NOTE: binary format is used if bin_format = true.

******************************************************************************/

FILE *file_open_output(const char filename[], const bool bin_format)
{
	FILE *output
		= (bin_format? fopen(filename, "wb") : fopen(filename, "w"));

	if (output == NULL)
	{
		PRINT_ERROR("unable to open %s\n", filename)
		exit(EXIT_FAILURE);
	}

	return output;
}

/******************************************************************************

 Function file_exist(): check if a file named filename exist.

******************************************************************************/

bool file_exist(const char filename[])
{
	return (access(filename, F_OK) == 0);
}

/******************************************************************************

 Function file_end(): check if a file has reached the end.

******************************************************************************/

bool file_end(FILE *input)
{
	return (feof(input) != 0);
}

/******************************************************************************

 Function file_write(): write n data elements with a given data size to a file
 using binary format.

 NOTE: data_size is often the output of sizeof().

******************************************************************************/

void file_write(const char filename[],
                const int n, const int data_size, const void *data[])
{
	ASSERT(n > 0)
	ASSERT(data != NULL)
	ASSERT(data_size > 0)

	FILE *output = file_open_output(filename, true);

	if (fwrite(&data, data_size, n, output) != (size_t) n)
	{
		PRINT_ERROR("unable to write all %d elements\n", n)
	}

	fclose(output);
}

/******************************************************************************

 Function file_read(): read n data elements with a given data size from a file
 using binary format.

 NOTE: data_size is often the output of sizeof().

******************************************************************************/

void *file_read(const char filename[], const int n, const int data_size)
{
	ASSERT(n > 0)
	ASSERT(data_size > 0)

	FILE *input = file_open_input(filename, true);

	void *data = malloc(data_size*n);
	ASSERT(data != NULL)

	if (fread(data, data_size, n, input) != (size_t) n)
	{
		PRINT_ERROR("unable to write all %d elements\n", n)
	}

	fclose(input);
	return data;
}

/******************************************************************************

 Function file_append(): append n data elements with a given data size to an
 already existing binary file.

 NOTE: data_size is often the output of sizeof().

******************************************************************************/

void file_append(const char filename[],
                 const int n, const int data_size, const void *data[])
{
	ASSERT(n > 0)
	ASSERT(data != NULL)
	ASSERT(data_size > 0)

	FILE *output = fopen(filename, "ab");

	if (output == NULL)
	{
		PRINT_ERROR("unable to open %s\n", filename)
		exit(EXIT_FAILURE);
	}

	if(fwrite(&data, data_size, n, output) != (size_t) n)
	{
		PRINT_ERROR("unable to write all %d elements\n", n)
	}

	fclose(output);
}

/******************************************************************************

 Function file_init_stdin(): set a file named filename as the C stdin.

******************************************************************************/

void file_init_stdin(const char filename[])
{
	FILE *input = file_open_input(filename, false);

	stdin = input;
	input = NULL;
}
