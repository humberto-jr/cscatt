#include "modules/file.h"
#include "modules/globals.h"

/******************************************************************************

 Function multipole_file(): opens the file for a given arrangement and grid
 index. Where, mode is the file access mode of fopen() from the C library.

 NOTE: mode = "wb" for write + binary format and "rb" for read + binary format,
 extension used is .bin, otherwise .dat is used assuming text mode ("w" or "r").

******************************************************************************/

FILE *multipole_file(const char arrang,
                     const int grid_index, const char mode[])
{
	char filename[MAX_LINE_LENGTH], ext[4];

	if (strlen(mode) > 1 && mode[1] == 'b')
		sprintf(ext, "%s", "bin");
	else
		sprintf(ext, "%s", "dat");

	sprintf(filename, "multipole_arrang=%c_n=%d.%s", arrang, grid_index, ext);

	FILE *stream = fopen(filename, mode);

	if (stream != NULL && mode[0] == 'r') printf("# Reading %s\n", filename);
	if (stream != NULL && mode[0] == 'w') printf("# Writing %s\n", filename);

	return stream;
}

/******************************************************************************

 Function multipole_file(): opens the file for a given arrangement and grid
 index. Where, mode is the file access mode of fopen() from the C library.

 NOTE: mode = "wb" for write + binary format and "rb" for read + binary format.

******************************************************************************/

FILE *open_datafile(const char arrang, const int grid_index, const int lambda)
{
	char filename[MAX_LINE_LENGTH];

	sprintf(filename, "multipole_arrang=%c_n=%d_lambda=%d.dat", arrang, grid_index, lambda);
	printf("# Writing %s\n", filename);

	return fopen(filename, "w");
}

/******************************************************************************

******************************************************************************/

int main(int argc, char *argv[])
{
	ASSERT(argc > 1)
	file_init_stdin(argv[1]);

	const char arrang = 96 + (int) file_keyword(stdin, "arrang", 1.0, 3.0, 1.0);

	int n = 0;
	FILE *input = NULL;

	while ((input = multipole_file(arrang, n, "rb")) != NULL)
	{
		double R = 0.0;
		file_read(&R, sizeof(double), 1, input, 0);

		double r_min = 0.0;
		file_read(&r_min, sizeof(double), 1, input, 0);

		double r_max = 0.0;
		file_read(&r_max, sizeof(double), 1, input, 0);

		double r_step = 0.0;
		file_read(&r_step, sizeof(double), 1, input, 0);

		int lambda_max = 0;
		file_read(&lambda_max, sizeof(int), 1, input, 0);

		int grid_size = 0;
		file_read(&grid_size, sizeof(int), 1, input, 0);

		int lambda = 0;
		file_read(&lambda, sizeof(int), 1, input, 0);

		while (!file_end(input))
		{
			FILE *output = open_datafile(arrang, n, lambda);

			fprintf(output, "# R = %06f\n", R);
			fprintf(output, "# File created on %s\n", time_stamp());

			for (int m = 0; m < grid_size; ++m)
			{
				const double r = r_min + as_double(m)*r_step;

				double multipole = 0.0;
				file_read(&multipole, sizeof(double), 1, input, 0);

/*
 *				NOTE: "% -8e\t" print numbers left-justified with an invisible
 *				plus sign, if any, 8 digits wide in scientific notation + tab.
 */

				fprintf(output, "%06f\t % -8e\t\n", r, multipole);
			}

			fclose(output);
			if (lambda == lambda_max) break;

			file_read(&lambda, sizeof(int), 1, input, 0);
		}

		fclose(input);
		++n;
	}

	return EXIT_SUCCESS;
}
