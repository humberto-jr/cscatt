#include "modules/fgh.h"
#include "modules/file.h"
#include "modules/matrix.h"
#include "modules/globals.h"

#define FORMAT "# %5zu  %5zu   %5zu   %5zu   %5zu   %+5d   %5zu     % -8e  % -8e  % -8e\n"

/******************************************************************************

 Function parity(): return the parity p = (-1)^l for a given l; either p = +1
 or p = -1.

******************************************************************************/

inline static int parity(const size_t l)
{
	return (int) pow(-1.0, l);
}

/******************************************************************************

 Function print_level(): prints the quantum numbers and eigenvalues for a given
 channel and total angular momentum J.

******************************************************************************/

void print_level(const fgh_basis *b, const size_t ch, const size_t J)
{
	printf(FORMAT, J, ch, b->v, b->j, b->l, parity(b->j + b->l), b->n,
	       b->eigenval, b->eigenval*219474.63137054, b->eigenval*27.211385);
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

	const size_t J_min = read_int_keyword(stdin, "J_min", 0, 10000, 0);

	const size_t J_max = read_int_keyword(stdin, "J_max", J_min, 10000, J_min);

	const size_t J_step = read_int_keyword(stdin, "J_step", 1, 10000, 1);

/*
 *	Vibrational grid:
 */

	const size_t n_max = read_int_keyword(stdin, "rovib_grid_size", 1, 10000, 1000);

	const double r_min = read_dbl_keyword(stdin, "r_min", 0.0, INF, 0.5);

	const double r_max = read_dbl_keyword(stdin, "r_max", r_min, INF, r_min + 30.0);

	const double r_step = (r_max - r_min)/as_double(n_max);

/*
 *	Arrangement (a = 1, b = 2, c = 3) and eigenvalue shift:
 */

	const char arrang = 96 + read_int_keyword(stdin, "arrang", 1, 3, 1);

	const double shift = read_dbl_keyword(stdin, "energy_shift", -INF, INF, 0.0);

/*
 *	Directory to load basis functions from:
 */

	char *dir = read_str_keyword(stdin, "basis_dir", ".");

/*
 *	Resize the basis functions for all J:
 */

	printf("#     J     ch.      v       j       l       p    comp.       E (a.u.)       E (cm-1)         E (eV)   \n");
	printf("# -----------------------------------------------------------------------------------------------------\n");

	for (size_t J = J_min; J <= J_max; J += J_step)
	{
		const size_t max_channel = fgh_basis_count(dir, arrang, J);

		for (size_t ch = 0; ch < max_channel; ++ch)
		{
			fgh_basis old, new;
			fgh_basis_load(&old, dir, arrang, ch, J);

			ASSERT(r_min >= old.r_min)
			ASSERT(r_max <= old.r_max)

			new.v = old.v;
			new.j = old.j;
			new.l = old.l;
			new.n = old.n;
			new.r_min = r_min;
			new.r_max = r_max;
			new.r_step = r_step;
			new.grid_size = n_max;
			new.eigenval = old.eigenval + shift;
			new.eigenvec = allocate(n_max, sizeof(double), false);

			for (size_t n = 0; n < n_max; ++n)
			{
				const double r = r_min + as_double(n)*r_step;

				new.eigenvec[n]
					= fgh_interpolation(old.grid_size, old.eigenvec, old.r_min, old.r_max, r);
			}

			print_level(&new, ch, J);

			fgh_basis_save(&new, ".", arrang, ch, J);

			free(old.eigenvec);
			free(new.eigenvec);
		}
	}

	free(dir);
	return EXIT_SUCCESS;
}
