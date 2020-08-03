#include "modules/file.h"
#include "modules/matrix.h"
#include "modules/globals.h"

#include "utils.h"

/******************************************************************************

 Function fgh_eigenvec(): interpolate the eigenvector a of the Hamiltonian
 built by dvr_fgh(), using Eq. (4.1) of Ref. [3], and return its amplitude
 value at a given new r.

******************************************************************************/

double fgh_eigenvec(const int grid_size,
                    const double r_min,
                    const double r_max,
                    const double r_new,
                    const double eigenvec[])
{
	ASSERT(grid_size > 0)
	ASSERT(eigenvec != NULL)

	const double r_step = (r_max - r_min)/as_double(grid_size);

	double result = 0.0;
	for (int n = 0; n < grid_size; ++n)
	{
		const double r_old = r_min + as_double(n)*r_step;
		const double param = M_PI*(r_new - r_old)/r_step;

		result += (fabs(param) > 1.0E-7? eigenvec[n]*sinc(param) : eigenvec[n]);
	}

	return result;
}

/******************************************************************************

******************************************************************************/

int main(int argc, char *argv[])
{
	ASSERT(argc > 1)
	file_init_stdin(argv[1]);

/*
 *	Total angular momentum, J:
 */

	const int J = (int) file_keyword(stdin, "J", 0.0, INF, 0.0);

/*
 *	Vibrational grid:
 */

	const int grid_size = (int) file_keyword(stdin, "rovib_grid_size", 1.0, INF, 500.0);

	const double r_min = file_keyword(stdin, "r_min", 0.0, INF, 0.5);

	const double r_max = file_keyword(stdin, "r_max", r_min, INF, r_min + 30.0);

	const double r_step = (r_max - r_min)/as_double(grid_size);

/*
 *	Arrangement (a = 1, b = 2, c = 3):
 */

	const char arrang = 96 + (int) file_keyword(stdin, "arrang", 1.0, 3.0, 1.0);

	const int n_max = basis_count(arrang, J);

	for (int n = 0; n < n_max; ++n)
	{
		basis old;
		basis_read(arrang, n, J, &old);

		FILE *output = basis_file(arrang, n, J, "wb");

		file_write(&old.v, sizeof(int), 1, output);

		file_write(&old.j, sizeof(int), 1, output);

		file_write(&old.l, sizeof(int), 1, output);

		file_write(&old.n, sizeof(int), 1, output);

		file_write(&r_min, sizeof(double), 1, output);

		file_write(&r_max, sizeof(double), 1, output);

		file_write(&r_step, sizeof(double), 1, output);

		file_write(&old.eigenval, sizeof(double), 1, output);

		file_write(&grid_size, sizeof(int), 1, output);

		for (int m = 0; m < grid_size; ++m)
		{
			const double r = r_min + as_double(m)*r_step;

			const double wavef = fgh_eigenvec(old.grid_size,
			                                  old.r_min, old.r_max, r, old.eigenvec);

			file_write(&wavef, sizeof(double), 1, output);
		}

		fclose(output);
		free(old.eigenvec);
	}

	return EXIT_SUCCESS;
}
