#include "modules/pes.h"
#include "modules/file.h"
#include "modules/globals.h"

int main(int argc, char *argv[])
{
	ASSERT(argc > 1)

	file_init_stdin(argv[1]);

/*
 *	Arrangement (a = 1, b = 2, c = 3) and directory to load multipoles from:
 */

	const char arrang = 96 + read_int_keyword(stdin, "arrang", 1, 3, 1);

	char *dir = read_str_keyword(stdin, "multipole_dir", ".");

/*
 *	Print the multipole coefficients for each R in a datafile with all r (lines)
 *	and lambda (columns) values:
 */

	const size_t n_max = pes_multipole_count(dir, arrang);

	for (size_t n = 0; n < n_max; ++n)
	{
		pes_multipole m;
		pes_multipole_load(&m, dir, arrang, n);

		FILE *output = pes_multipole_file(".", arrang, n, "w", true);

		fprintf(output, "# R = %06f\n", m.R);
		fprintf(output, "# lambda = (%zu, %zu, %zu)\n", m.lambda_min, m.lambda_max, m.lambda_step);
		fprintf(output, "# File created at %s\n", time_stamp());

		for (size_t p = 0; p < m.grid_size; ++p)
		{
			fprintf(output, "%06f\t", m.r_min + as_double(p)*m.r_step);
/*
 *			NOTE: "% -8e\t" print numbers left-justified with an invisible
 *			plus sign, if any, 8 digits wide in scientific notation + tab.
*/
			for (size_t lambda = m.lambda_min; lambda <= m.lambda_max; lambda += m.lambda_step)
				fprintf(output, " % -8e\t", m.value[lambda][p]);

			fprintf(output, "\n");
		}

		file_close(&output);
		pes_multipole_free(&m);
	}

	free(dir);
	return EXIT_SUCCESS;
}
