#include "modules/file.h"
#include "modules/matrix.h"
#include "modules/globals.h"

#include "basis_config.h"

int main(int argc, char *argv[])
{
	ASSERT(argc > 1)
	file_init_stdin(argv[1]);

	const int J
		= (int) file_get_key(stdin, "J", 0.0, INF, 0.0);

	const char arrang
		= 96 + (int) file_get_key(stdin, "arrang", 1.0, 3.0, 1.0);

	int ch = 0;
	while (check_basis(arrang, ch, J))
	{
		scatt_basis basis;
		load_basis(arrang, ch, J, &basis);

		const int grid_size
			 = matrix_row(basis.wavef)/basis.max_state;

		char filename[MAX_LINE_LENGTH];
		sprintf(filename, BASIS_DATAFILE_FORMAT, arrang, ch, J);

		FILE *output
			= file_open_output(filename, false, false);

		printf("# Writing %s\n", filename);

		fprintf(output, "# v = %d\n", basis.v);
		fprintf(output, "# j = %d\n", basis.j);
		fprintf(output, "# l = %d\n", basis.l);
		fprintf(output, "# Spin mult. = %d\n", basis.spin_mult);
		fprintf(output, "# Components = %d\n", basis.max_state);
		fprintf(output, "# Eigenvalue = % -8e\n", basis.energy);
		fprintf(output, "# File created on %s\n", time_stamp());
		fprintf(output, "#\n");

		for (int n = 0; n < grid_size; ++n)
		{
			const double r
				= basis.r_min + as_double(n)*basis.r_step;

/*
 *			NOTE: "% -8e\t" print numbers left-justified with an invisible
 *			plus sign, if any, 8 digits wide in scientific notation + tab.
 */

			if (basis.max_state == 1)
			{
				fprintf(output, "%06f\t % -8e\t\n",
				        r, matrix_get(basis.wavef, n, 0));
			}
			else
			{
				fprintf(output, "%06f\t ", r);

				for (int m = 0; m < basis.max_state; ++m)
				{
					fprintf(output, "% -8e\t ",
					        matrix_get(basis.wavef, m*grid_size + n, 0));
				}

				fprintf(output, "\n");
			}
		}

		matrix_free(basis.wavef);
		fclose(output);
		++ch;
	}

	return EXIT_SUCCESS;
}
