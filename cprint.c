#include "modules/file.h"
#include "modules/matrix.h"
#include "modules/globals.h"

#include "coupl_config.h"

int main(int argc, char *argv[])
{
	ASSERT(argc > 1)
	file_init_stdin(argv[1]);

	const int J
		= (int) file_get_key(stdin, "J", 0.0, INF, 0.0);

	const int scatt_grid_size
		= (int) file_get_key(stdin, "scatt_grid_size", 1.0, INF, 500.0);

	const double R_min
		= file_get_key(stdin, "R_min", 0.0, INF, 0.5);

	const double R_max
		= file_get_key(stdin, "R_max", R_min, INF, R_min + 100.0);

	const double R_step
		= (R_max - R_min)/as_double(scatt_grid_size);

	const double shift
		= file_get_key(stdin, "energy_shift", -INF, INF, 0.0);

	const double scale
		= file_get_key(stdin, "energy_scale", -INF, INF, 1.0);

	const char arrang
		= 96 + (int) file_get_key(stdin, "arrangement", 1.0, 3.0, 1.0);

	const bool adiabatic
		= (bool) file_get_key(stdin, "print_adiabatic", 0.0, 1.0, 0.0);

	printf("# J                 = %d\n", J);
	printf("# Arrangement       = %c\n", arrang);
	printf("# Energy shift      = %f\n", shift);
	printf("# Energy scale      = %f\n", scale);
	printf("# Grid points       = %d\n", scatt_grid_size);
	printf("# Representation    = %s\n", (adiabatic? "adiabatic" : "diabatic"));
	printf("#\n");

	int n = 0;

	char filename[MAX_LINE_LENGTH];
	sprintf(filename, CMATRIX_BUFFER_FORMAT, arrang, n, J);

	while (file_exist(filename))
	{
		matrix *c = matrix_load(filename);

		printf("# %5d: reading %s, size = %dx%d\n",
		       n, filename, matrix_row(c), matrix_col(c));

		double *eigenval = NULL;
		if (adiabatic) eigenval = matrix_symm_eigen(c, 'n');

		for (int ch_a = 0; ch_a < matrix_row(c); ++ch_a)
		{
			memset(filename, 0, sizeof(filename));
			sprintf(filename, CMATRIX_DATAFILE_FORMAT, arrang, ch_a, J);

			if (file_exist(filename) && n == 0) remove(filename);
			FILE *output = file_open_output(filename, false, true);

			fprintf(output, "  % 6f\t", R_min + as_double(n)*R_step);

/*
 *			NOTE: "% -8e\t" print numbers left-justified with an invisible
 *			plus sign, if any, 8 digits wide in scientific notation + tab.
 */

			if (eigenval != NULL)
			{
				fprintf(output, "% -8e\n", (eigenval[ch_a] + shift)*scale);
			}
			else
			{
				fprintf(output, "% -8e\t", (matrix_get(c, ch_a, ch_a) + shift)*scale);

				for (int ch_b = (ch_a + 1); ch_b < matrix_col(c); ++ch_b)
				{
					fprintf(output, "% -8e\t", matrix_get(c, ch_a, ch_b)*scale);
				}

				fprintf(output, "\n");
			}

			fclose(output);
		}

		matrix_free(c);
		if (eigenval != NULL) free(eigenval);

		++n;
		memset(filename, 0, sizeof(filename));
		sprintf(filename, CMATRIX_BUFFER_FORMAT, arrang, n, J);
	}

	return EXIT_SUCCESS;
}
