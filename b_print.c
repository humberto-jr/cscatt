#include "modules/file.h"
#include "modules/matrix.h"
#include "modules/globals.h"

#include "utils.h"

/******************************************************************************
******************************************************************************/

int main(int argc, char *argv[])
{
	ASSERT(argc > 1)
	file_init_stdin(argv[1]);

/*
 *	Total angular momentum, J:
 */

	const int J_min = (int) file_keyword(stdin, "J_min", 0.0, INF, 0.0);

	const int J_max = (int) file_keyword(stdin, "J_max", 0.0, INF, 0.0);

	const int J_step = (int) file_keyword(stdin, "J_step", 1.0, INF, 0.0);

	ASSERT(J_max >= J_min)

/*
 *	Arrangement (a = 1, b = 2, c = 3):
 */

	const char arrang = 96 + (int) file_keyword(stdin, "arrang", 1.0, 3.0, 1.0);

/*
 *	Print the basis functions for each J:
 */

	for (int J = J_min; J <= J_max; J += J_step)
	{
		const int max_ch = basis_count(arrang, J);

		for (int ch = 0; ch < max_ch; ++ch)
		{
			basis b;
			basis_read(arrang, ch, J, &b, true);

			FILE *output = basis_file(arrang, ch, J, "w", true);

			ASSERT(output != NULL)

			fprintf(output, "# v = %d\n", b.v);
			fprintf(output, "# j = %d\n", b.j);
			fprintf(output, "# l = %d\n", b.l);

			fprintf(output, "# Component  = %d\n", b.n);
			fprintf(output, "# Eigenvalue = % -8e\n", b.eigenval);

			fprintf(output, "# File created on %s\n", time_stamp());

			for (int n = 0; n < b.grid_size; ++n)
			{
				const double r = b.r_min + as_double(n)*b.r_step;

/*
 *				NOTE: "% -8e\t" print numbers left-justified with an invisible
 *				plus sign, if any, 8 digits wide in scientific notation + tab.
 */

				fprintf(output, "%06f\t % -8e\t\n", r, b.eigenvec[n]);
			}

			fclose(output);
			free(b.eigenvec);
		}
	}

	return EXIT_SUCCESS;
}
