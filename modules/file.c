/******************************************************************************

 About
 -----

 This module is a collection of routines for file manipulation all based upon
 the C library.

******************************************************************************/

#include "file.h"

/******************************************************************************

 Function file_exist(): check if a file named filename exist.

******************************************************************************/

bool file_exist(const char filename[])
{
	return (access(filename, F_OK) == 0);
}

/******************************************************************************

 Function file_open_input(): open a file named filename in read-only mode.

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

 Function file_open_input(): open a file named filename in write-only mode.

 NOTE: binary format is used if bin_format = true.

******************************************************************************/

FILE *file_open_output(const char filename[],
                       const bool bin_format, const bool to_append)
{
	FILE *output = NULL;

	if (to_append)
		output = (bin_format? fopen(filename, "ab") : fopen(filename, "a"));
	else
		output = (bin_format? fopen(filename, "wb") : fopen(filename, "w"));

	if (output == NULL)
	{
		PRINT_ERROR("unable to open %s\n", filename)
		exit(EXIT_FAILURE);
	}

	return output;
}

/******************************************************************************

 Wrapper file_close(): an interface to fclose() that also sets the FILE pointer
 to NULL.

******************************************************************************/

void file_close(FILE **stream)
{
	if (*stream != NULL)
	{
		fclose(*stream);
		*stream = NULL;
	}
}

/******************************************************************************

 Function file_init_stdin(): set a file named filename as the C stdin.

******************************************************************************/

void file_init_stdin(const char filename[])
{
	FILE *input = file_open_input(filename, false);

	fclose(stdin);
	stdin = input;
	input = NULL;
}

/******************************************************************************

 Function file_init_stdout(): set a file named filename as the C stdout.

******************************************************************************/

void file_init_stdout(const char filename[], const bool to_append)
{
	FILE *output = file_open_output(filename, false, to_append);

	fclose(stdout);
	stdout = output;
	output = NULL;
}

/******************************************************************************

 Function file_find(): scans a given input file searching for the first
 occurrence of a given pattern. It shall return the whole line, if found, or
 an empty string otherwise.

 NOTE: Lines starting by '#' are ignored.

******************************************************************************/

char *file_find(FILE *input, const char pattern[])
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
 occurrence of a keyword with format "[key] = [value]". If found, and if [min]
 <= [value] <= [max], it shall return [value]. Otherwise, it returns a default
 value.

 NOTE: Lines starting by '#' are ignored.

******************************************************************************/

double file_get_key(FILE *input, const char key[],
                    const double min, const double max, const double default_value)
{
	ASSERT(max > min)

	char *line = file_find(input, key);
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

 Function file_get_key(): scans a given input file searching for the first
 occurrence of a keyword with format "[key] = [value]". If found, and if [min]
 <= [value] <= [max], it shall return [value]. Otherwise, it returns a default
 value.

 NOTE: Lines starting by '#' are ignored.

******************************************************************************/

double file_keyword(FILE *input, const char key[],
                    const double min, const double max, const double default_value)
{
	ASSERT(max > min)

	char *line = file_find(input, key);
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

 Function file_end(): check if a file has reached the end.

******************************************************************************/

bool file_end(FILE *stream)
{
	ASSERT(stream != NULL)
	return (feof(stream) != 0);
}

/******************************************************************************

 Wrapper file_write(): is an interface for fwrite() that check inputs and the
 output for errors, if any.

 NOTE: data_size is often the output of sizeof() for the type of data.

******************************************************************************/

void file_write(const void *data, const int data_size, const int n, FILE *stream)
{
	ASSERT(n > 0)
	ASSERT(data_size > 0)

	ASSERT(data != NULL)
	ASSERT(stream != NULL)

	const int info = fwrite(data, data_size, n, stream);

	if (info < n) PRINT_ERROR("only %d/%d elements written\n", info, n)
}

/******************************************************************************

 Wrapper file_read(): is an interface for fread() which not only check inputs,
 but apply a byte offset before reading, if desired, and check the output for
 errors, if any.

 NOTE: data_size is often the output of sizeof() for the type of data and if
 offset is zero the reading is performed from the beginning of the file.

******************************************************************************/

void file_read(void *data,
               const int data_size, const int n, FILE *stream, const int offset)
{
	ASSERT(n > 0)
	ASSERT(data_size > 0)

	ASSERT(data != NULL)
	ASSERT(stream != NULL)

	if (offset > 0) fseek(stream, offset, SEEK_SET);

	const int info = fread(data, data_size, n, stream);

	if (info < n) PRINT_ERROR("only %d/%d elements read\n", info, n)
}
