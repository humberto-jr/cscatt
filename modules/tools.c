#include <omp.h>
#include <time.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>

#include "globals.h"
#include "tools.h"

/******************************************************************************

 Function tools_find_string(): scans a given input file searching for the first
 occurrence of a given pattern. If found, it shall return the whole line, or an
 empty string otherwise.

 NOTE: Lines starting by '#' are ignored.

******************************************************************************/

char *tools_find_string(FILE *input, const char pattern[])
{
	char *line = malloc(sizeof(char)*MAX_LINE_LENGTH);

	ASSERT(line != NULL)
	ASSERT(input != NULL)
	ASSERT(pattern != NULL)

	rewind(input);
	while (fgets(line, MAX_LINE_LENGTH, input) != NULL)
	{
		if (line[0] != '#' && strstr(line, pattern) != NULL) return line;
	}

	free(line);
	return "\n";
}

/******************************************************************************

 Function tools_get_key(): scans a given input file searching for the first
 occurrence of a keyword with the format "[key] = [value]". If found, and
 if [min] <= [value] <= [max], it shall return [value]. Otherwise, it
 returns a default value.

 NOTE: Lines starting by '#' are ignored.

******************************************************************************/

double tools_get_key(FILE *input, const char key[],
                     const double min, const double max, const double default_value)
{
	ASSERT(max > min)

	char *line = tools_find_string(input, key);
	char *token = strtok(line, "=");

	while (token != NULL)
	{
		if (strstr(line, key) != NULL)
		{
			token = strtok(NULL, "=");
			register const double value = atof(token);

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

 Function tools_row_count(): counts the actual number of valid lines in a given
 input file.

 NOTE: Blank lines or lines starting by '#' are ignored.

******************************************************************************/

int tools_row_count(FILE *input)
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

 Function tools_col_count(): counts the actual number of columns in a given
 input file.

 NOTE: Blank lines or lines starting by '#' are ignored.

******************************************************************************/

int tools_col_count(FILE *input)
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

 Function tools_time_stamp(): return the current YMDHMS date as a time stamp.

 Example: 31 May, 2001 09:45:54 AM.

******************************************************************************/

const char *tools_time_stamp()
{
	static char stamp[50];

	time_t now = time(NULL);
	const struct tm *info = localtime(&now);
	strftime(stamp, 50, "%B %d, %Y %I:%M:%S %p", info);

	return stamp;
}

/******************************************************************************

 Function tools_wall_time(): return the wall time in units of seconds.

******************************************************************************/

double tools_wall_time()
{
	return omp_get_wtime();
}
