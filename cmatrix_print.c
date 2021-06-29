#include "modules/file.h"
#include "modules/matrix.h"
#include "modules/globals.h"

int main(int argc, char *argv[])
{
	ASSERT(argc > 1)

	file_init_stdin(argv[1]);

/*
 *	Total angular momentum, J:
 */

	const size_t J_min = read_int_keyword(stdin, "J_min", 0, 10000, 0);

	const size_t J_max = read_int_keyword(stdin, "J_max", J_min, 10000, J_min);

	const size_t J_step = read_int_keyword(stdin, "J_step", 1, 10000, 1);

/*
 *	Scattering grid:
 */

	const size_t n_max = read_int_keyword(stdin, "scatt_grid_size", 1, 1000000, 500);

	const double R_min = read_dbl_keyword(stdin, "R_min", 0.0, INF, 0.5);

	const double R_max = read_dbl_keyword(stdin, "R_max", R_min, INF, R_min + 30.0);

	const double R_step = (R_max - R_min)/as_double(n_max);

/*
 *	Energy scale and shifts:
 */

	const double shift = read_dbl_keyword(stdin, "energy_shift", -INF, INF, 0.0);

	const double scale = read_dbl_keyword(stdin, "energy_scale", -INF, INF, 1.0);

/*
 *	Arrangement (a = 1, b = 2, c = 3) and representation:
 */

	const char arrang = 96 + read_int_keyword(stdin, "arrang", 1, 3, 1);

	const bool adiabatic = read_int_keyword(stdin, "print_adiabatic", 0, 1, 0);

/*
 *	Print output:
 */

	char filename[MAX_LINE_LENGTH];

	for (size_t J = J_min; J <= J_max; J += J_step)
	{
		for (size_t n = 0; n < n_max; ++n)
		{
			sprintf(filename, "cmatrix_arrang=%c_n=%zu_J=%zu.bin", arrang, n, J);

			matrix *c = matrix_load(filename);

			double *eigenval = NULL;
			if (adiabatic) eigenval = matrix_symm_eigen(c, 'n');

			for (size_t a = 0; a < matrix_rows(c); ++a)
			{
				sprintf(filename, "cmatrix_arrang=%c_n=%zu_J=%zu.dat", arrang, n, J);

				FILE *output = file_open(filename, "w");

				if (n == 0)
				{
					fprintf(output, "# Energy shift = %f\n", shift);
					fprintf(output, "# Energy scale = %f\n", scale);
					fprintf(output, "# Representation = %s\n", (adiabatic? "adiabatic" : "diabatic"));
					fprintf(output, "# File created at %s\n", time_stamp());
				}

				fprintf(output, "%06f\t ", R_min + as_double(n)*R_step);

/*
 *				NOTE: "% -8e\t" print numbers left-justified with an invisible
 *				plus sign, if any, 8 digits wide in scientific notation + tab.
 */

				if (eigenval != NULL)
				{
					fprintf(output, "% -8e\n", (eigenval[a] + shift)*scale);
				}
				else
				{
					fprintf(output, "% -8e\t", (matrix_get(c, a, a) + shift)*scale);

					for (size_t b = (a + 1); b < matrix_cols(c); ++b)
						fprintf(output, "% -8e\t", matrix_get(c, a, b)*scale);

					fprintf(output, "\n");
				}

				file_close(&output);
			}

			matrix_free(c);
			if (eigenval != NULL) free(eigenval);
		}
	}

	return EXIT_SUCCESS;
}
