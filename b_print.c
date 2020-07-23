#include "modules/file.h"
#include "modules/globals.h"

//#include "basis_config.h"

/******************************************************************************

 Function basis_file(): opens the file for the n-th channel, arrangement and
 total angular momentum, J. Where, mode is the file access mode of fopen()
 from the C library.

 NOTE: mode = "wb" for write + binary format and "rb" for read + binary format,
 extension used is .bin, otherwise .dat is used assuming text mode ("w" or "r").

******************************************************************************/

FILE *basis_file(const char arrang,
                 const int n, const int J, const char mode[])
{
	char filename[MAX_LINE_LENGTH], ext[4];

	if (strlen(mode) > 1 && mode[1] == 'b')
		sprintf(ext, "%s", "bin");
	else
		sprintf(ext, "%s", "dat");

	sprintf(filename, "basis_arrang=%c_ch=%d_J=%d.%s", arrang, n, J, ext);

	FILE *stream = fopen(filename, mode);

	if (stream != NULL && mode[0] == 'r') printf("# Reading %s\n", filename);
	if (stream != NULL && mode[0] == 'w') printf("# Writing %s\n", filename);

	return stream;
}

/******************************************************************************

******************************************************************************/

int main(int argc, char *argv[])
{
	ASSERT(argc > 1)
	file_init_stdin(argv[1]);

	const int J = (int) file_keyword(stdin, "J", 0.0, INF, 0.0);

	const char arrang = 96 + (int) file_keyword(stdin, "arrang", 1.0, 3.0, 1.0);

	int n = 0;
	FILE *input = NULL;

	while ((input = basis_file(arrang, n, J, "rb")) != NULL)
	{
		int v = 0;
		file_read(&v, sizeof(int), 1, input, 0);

		int j = 0;
		file_read(&j, sizeof(int), 1, input, 0);

		int l = 0;
		file_read(&l, sizeof(int), 1, input, 0);

		int i = 0;
		file_read(&i, sizeof(int), 1, input, 0);

		double r_min = 0.0;
		file_read(&r_min, sizeof(double), 1, input, 0);

		double r_max = 0.0;
		file_read(&r_max, sizeof(double), 1, input, 0);

		double r_step = 0.0;
		file_read(&r_step, sizeof(double), 1, input, 0);

		double energy = 0.0;
		file_read(&energy, sizeof(double), 1, input, 0);

		int n_max = 0;
		file_read(&n_max, sizeof(int), 1, input, 0);

		FILE *output = basis_file(arrang, n, J, "w");

		ASSERT(output != NULL)

		fprintf(output, "# v = %d\n", v);
		fprintf(output, "# j = %d\n", j);
		fprintf(output, "# l = %d\n", l);

		fprintf(output, "# Component  = %d\n", i);
		fprintf(output, "# Eigenvalue = % -8e\n", energy);

		fprintf(output, "# File created on %s\n", time_stamp());

		for (int n = 0; n < n_max; ++n)
		{
			const double r = r_min + as_double(n)*r_step;

			double wavef = 0.0;
			file_read(&wavef, sizeof(double), 1, input, 0);

/*
 *			NOTE: "% -8e\t" print numbers left-justified with an invisible
 *			plus sign, if any, 8 digits wide in scientific notation + tab.
 */

			fprintf(output, "%06f\t % -8e\t\n", r, wavef);
		}

		fclose(output);
		++n;
	}

	return EXIT_SUCCESS;
}
