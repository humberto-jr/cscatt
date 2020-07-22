#include "modules/file.h"
#include "modules/globals.h"

#include "basis_config.h"

int main(int argc, char *argv[])
{
	ASSERT(argc > 1)
	file_init_stdin(argv[1]);

	const int J = (int) file_keyword(stdin, "J", 0.0, INF, 0.0);

	const char arrang = 96 + (int) file_keyword(stdin, "arrang", 1.0, 3.0, 1.0);

	int ch = 0;
	while (check_basis_file(arrang, ch, J) == true)
	{
		scatt_basis basis;
		read_basis_file(arrang, ch, J, &basis);

		char filename[MAX_LINE_LENGTH];
		sprintf(filename, BASIS_DATAFILE_FORMAT, arrang, ch, J);

		FILE *output = file_open_output(filename, false, false);

		printf("# Writing %s\n", filename);

		fprintf(output, "# v = %d\n", basis.v);
		fprintf(output, "# j = %d\n", basis.j);
		fprintf(output, "# l = %d\n", basis.l);

		fprintf(output, "# Component  = %d\n", basis.state);
		fprintf(output, "# Eigenvalue = % -8e\n", basis.energy);
		fprintf(output, "# File created on %s\n", time_stamp());

		fprintf(output, "#\n");

		for (int n = 0; n < basis.grid_size; ++n)
		{
			const double r = basis.r_min + as_double(n)*basis.grid_step;

/*
 *			NOTE: "% -8e\t" print numbers left-justified with an invisible
 *			plus sign, if any, 8 digits wide in scientific notation + tab.
 */

			fprintf(output, "%06f\t % -8e\t\n", r, basis.wavef[n]);
		}

		free(basis.wavef);
		fclose(output);
		++ch;
	}

	return EXIT_SUCCESS;
}
