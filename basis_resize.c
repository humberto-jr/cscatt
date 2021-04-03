#include "modules/fgh.h"
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

	const int J_min
		= read_int_keyword(stdin, "J_min", 0, 10000, 0);

	const int J_max
		= read_int_keyword(stdin, "J_max", J_min, 10000, J_min);

	const int J_step
		= read_int_keyword(stdin, "J_step", 1, 10000, 1);

/*
 *	Vibrational grid:
 */

	const int grid_size
		= read_int_keyword(stdin, "rovib_grid_size", 1, 10000, 1000);

	const double r_min
		= read_dbl_keyword(stdin, "r_min", 0.0, INF, 0.5);

	const double r_max
		= read_dbl_keyword(stdin, "r_max", r_min, INF, r_min + 30.0);

	const double r_step
		= (r_max - r_min)/as_double(grid_size);

/*
 *	Arrangement (a = 1, b = 2, c = 3):
 */

	const char arrang = 96 + read_int_keyword(stdin, "arrang", 1, 3, 1);

/*
 *	Directory to load basis functions from:
 */

	char *dir = read_str_keyword(stdin, "basis_dir", ".");

/*
 *	Resize the basis functions for all J:
 */

	for (int J = J_min; J <= J_max; J += J_step)
	{
		const int max_channel = fgh_basis_count(dir, arrang, J);

		for (int ch = 0; ch < max_channel; ++ch)
		{
			fgh_basis old;

			FILE *input = fgh_basis_file(dir, arrang, ch, J, "rb", true);

			fgh_basis_read(&old, input);

			file_close(&input);

			ASSERT(r_min >= old.r_min)
			ASSERT(r_max <= old.r_max)

			FILE *output = fgh_basis_file(dir, arrang, ch, J, "wb", true);

			file_write(&old.v, sizeof(int), 1, output);

			file_write(&old.j, sizeof(int), 1, output);

			file_write(&old.l, sizeof(int), 1, output);

			file_write(&old.n, sizeof(int), 1, output);

			file_write(&r_min, sizeof(double), 1, output);

			file_write(&r_max, sizeof(double), 1, output);

			file_write(&r_step, sizeof(double), 1, output);

			file_write(&old.eigenval, sizeof(double), 1, output);

			file_write(&grid_size, sizeof(int), 1, output);

			for (int n = 0; n < grid_size; ++n)
			{
				const double r = r_min + as_double(n)*r_step;

				const double wavef
					= fgh_interpolation(old.grid_size, old.r_min, old.r_max, r, old.eigenvec);

				file_write(&wavef, sizeof(double), 1, output);
			}

			file_close(&output);
			free(old.eigenvec);
		}
	}

	free(dir);
	return EXIT_SUCCESS;
}
