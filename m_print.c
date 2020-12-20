#include "modules/pes.h"
#include "modules/file.h"
#include "modules/globals.h"

/******************************************************************************
******************************************************************************/

int main(int argc, char *argv[])
{
	ASSERT(argc > 1)
	file_init_stdin(argv[1]);

	const char arrang = 96 + (int) file_keyword(stdin, "arrang", 1.0, 3.0, 1.0);

	int n = 0;
	FILE *input = NULL;

	while ((input = pes_multipole_file(arrang, n, "rb", true)) != NULL)
	{
		pes_multipole m;
		pes_multipole_read(&m, input);

		FILE *output = pes_multipole_file(arrang, n, "w", true);

		fprintf(output, "# R = %06f\n", m.R);
		fprintf(output, "# File created at %s\n", time_stamp());

		for (int p = 0; p < m.grid_size; ++p)
		{
			fprintf(output, "%06f\t", m.r_min + as_double(p)*m.r_step);

			for (int lambda = m.lambda_min; lambda <= m.lambda_max; lambda += m.lambda_step)
			{
/*
 *				NOTE: "% -8e\t" print numbers left-justified with an invisible
 *				plus sign, if any, 8 digits wide in scientific notation + tab.
 */
				fprintf(output, " % -8e\t", m.value[lambda][p]);
			}

			fprintf(output, "\n");
		}

		pes_multipole_free(&m);
		fclose(output);
		fclose(input);
		input = NULL;
		++n;
	}

	return EXIT_SUCCESS;
}
