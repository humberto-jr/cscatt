#include "modules/file.h"
#include "modules/matrix.h"
#include "modules/globals.h"

#include "utils.h"

int main(int argc, char *argv[])
{
	ASSERT(argc > 1)
	file_init_stdin(argv[1]);

/*
 *	Total angular momentum, J:
 */

	const int J = (int) file_keyword(stdin, "J", 0.0, INF, 0.0);

/*
 *	Scattering grid:
 */

	const int n_max = (int) file_keyword(stdin, "scatt_grid_size", 1.0, INF, 500.0);

	const double R_min = file_keyword(stdin, "R_min", 0.0, INF, 0.5);

	const double R_max = file_keyword(stdin, "R_max", R_min, INF, R_min + 30.0);

	const double R_step = (R_max - R_min)/as_double(n_max);

/*
 *	Energy scale and shifts:
 */

	const double shift = file_keyword(stdin, "energy_shift", -INF, INF, 0.0);

	const double scale = file_keyword(stdin, "energy_scale", -INF, INF, 1.0);

/*
 *	Arrangement (a = 1, b = 2, c = 3) and representation:
 */

	const char arrang = 96 + (int) file_keyword(stdin, "arrang", 1.0, 3.0, 1.0);

	const bool adiabatic = (bool) file_keyword(stdin, "print_adiabatic", 0.0, 1.0, 0.0);

/*
 *	Print output:
 */

	for (int n = 0; n < n_max; ++n)
	{
		matrix *c = coupling_read(arrang, n, true, J);

		double *eigenval = NULL;
		if (adiabatic) eigenval = matrix_symm_eigen(c, 'n');

		for (int a = 0; a < matrix_row(c); ++a)
		{
			FILE *output = coupling_datafile(arrang, a, J, false, "w");

			if (n == 0)
			{
				fprintf(output, "# Energy shift = %f\n", shift);
				fprintf(output, "# Energy scale = %f\n", scale);
				fprintf(output, "# Representation = %s\n", (adiabatic? "adiabatic" : "diabatic"));
				fprintf(output, "# File created on %s\n", time_stamp());
			}

			fprintf(output, "%06f\t ", R_min + as_double(n)*R_step);

/*
 *			NOTE: "% -8e\t" print numbers left-justified with an invisible
 *			plus sign, if any, 8 digits wide in scientific notation + tab.
 */

			if (eigenval != NULL)
			{
				fprintf(output, "% -8e\n", (eigenval[a] + shift)*scale);
			}
			else
			{
				fprintf(output, "% -8e\t", (matrix_get(c, a, a) + shift)*scale);

				for (int b = (a + 1); b < matrix_col(c); ++b)
				{
					fprintf(output, "% -8e\t", matrix_get(c, a, b)*scale);
				}

				fprintf(output, "\n");
			}

			fclose(output);
		}

		matrix_free(c);
		if (eigenval != NULL) free(eigenval);
	}

	return EXIT_SUCCESS;
}
