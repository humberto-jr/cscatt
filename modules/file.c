/******************************************************************************

 About
 -----

 This module is a collection of routines for file manipulation all based upon
 the C library.

******************************************************************************/

#include "file.h"

/******************************************************************************

 Wrapper file_open(): a safe interface to fopen() that halt the execution in
 the case of error and prints a formatted message in the C stderr.

******************************************************************************/

FILE *file_open(const char filename[], const char mode[])
{
	FILE *stream = fopen(filename, mode);

	if (stream == NULL)
	{
		PRINT_ERROR("unable to open %s with mode %s\n", filename, mode)
		exit(EXIT_FAILURE);
	}

	return stream;
}

/******************************************************************************

 Wrapper file_close(): a safe interface to fclose() that closes a file stream
 and sets the unused pointer to NULL.

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

 Wrapper file_exist(): uses access() to check if a file named filename exists.

******************************************************************************/

bool file_exist(const char filename[])
{
	return (access(filename, F_OK) == 0);
}

/******************************************************************************

 Wrapper file_end(): an interface to feof() that checks if a given file stream
 is at the end.

******************************************************************************/

bool file_end(FILE *stream)
{
	ASSERT(stream != NULL)
	return (feof(stream) != 0);
}

/******************************************************************************

 Wrapper file_delete(): a safe interface to remove() that halt the execution in
 the case of error and prints a formatted message in the C stderr.

******************************************************************************/

void file_delete(const char filename[])
{
	if (remove(filename) == -1)
	{
		PRINT_ERROR("unable to delete %s\n", filename)
		exit(EXIT_FAILURE);
	}
}

/******************************************************************************

 Wrapper file_rename(): a safe interface to rename() that halt the execution in
 the case of error and prints a formatted message in the C stderr.

******************************************************************************/

void file_rename(const char old_filename[], const char new_filename[])
{
	if (rename(old_filename, new_filename) == -1)
	{
		PRINT_ERROR("unable to rename %s to %s\n", old_filename, new_filename)
		exit(EXIT_FAILURE);
	}
}

/******************************************************************************

 Function file_init_stdin(): set a file named filename as the C stdin.

******************************************************************************/

void file_init_stdin(const char filename[])
{
	FILE *input = file_open(filename, "r");

	fclose(stdin);
	stdin = input;
	input = NULL;
}

/******************************************************************************

 Function file_init_stdout(): set a file named filename as the C stdout.

******************************************************************************/

void file_init_stdout(const char filename[])
{
	FILE *output = file_open(filename, "w");

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
	ASSERT(input != NULL)

	char *line = allocate(MAX_LINE_LENGTH, sizeof(char), false);

	rewind(input);
	while (fgets(line, MAX_LINE_LENGTH, input) != NULL)
	{
		if (line[0] != '#' && strstr(line, pattern) != NULL) return line;
	}

	free(line);
	return "\n";
}

/******************************************************************************

 Function file_keyword(): scans a given input file searching for the first
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

 Wrapper file_write(): is an interface to fwrite() that checks both inputs and
 outputs and halt the execution in the case of errors, printing a formatted
 message in the C stderr.

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

 Wrapper file_read(): is an interface to fread() that checks both inputs and
 outputs and halt the execution in the case of errors, printing a formatted
 message in the C stderr.

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
