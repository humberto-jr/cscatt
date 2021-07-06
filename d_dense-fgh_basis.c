#include "modules/pes.h"
#include "modules/fgh.h"
#include "modules/file.h"
#include "modules/matrix.h"
#include "modules/globals.h"

#define FORMAT "# %5zu   %5zu   %5zu   %5zu   %5zu   %+5d     % -8e  % -8e  % -8e\n"

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
	printf(FORMAT, J, ch, b->v, b->j, b->l, parity(b->j + b->l),
	       b->eigenval, b->eigenval*219474.63137054, b->eigenval*27.211385);
}

/******************************************************************************
******************************************************************************/

int main(int argc, char *argv[])
{
	ASSERT(argc > 1)

	matrix_init_gpu();

	file_init_stdin(argv[1]);

/*
 *	Total angular momentum, J:
 */

	const size_t J_min = read_int_keyword(stdin, "J_min", 0, 10000, 0);

	const size_t J_max = read_int_keyword(stdin, "J_max", J_min, 10000, J_min);

	const size_t J_step = read_int_keyword(stdin, "J_step", 1, 10000, 1);

	const int J_parity = read_int_keyword(stdin, "J_parity", -1, 1, 0);

/*
 *	Vibrational quantum numbers, v:
 */

	const size_t v_min = read_int_keyword(stdin, "v_min", 0, 10000, 0);

	const size_t v_max = read_int_keyword(stdin, "v_max", v_min, 10000, v_min);

	const size_t v_step = read_int_keyword(stdin, "v_step", 1, 10000, 1);

/*
 *	Rotational quantum numbers, j:
 */

	const size_t j_min = read_int_keyword(stdin, "j_min", 0, 10000, 0);

	const size_t j_max = read_int_keyword(stdin, "j_max", j_min, 10000, j_min);

	const size_t j_step = read_int_keyword(stdin, "j_step", 1, 10000, 1);

/*
 *	Vibrational grid:
 */

	const size_t n_max = read_int_keyword(stdin, "rovib_grid_size", v_max + 1, 1000000, 1000);

	const double r_min = read_dbl_keyword(stdin, "r_min", 0.0, INF, 0.5);

	const double r_max = read_dbl_keyword(stdin, "r_max", r_min, INF, r_min + 30.0);

	const double r_step = (r_max - r_min)/as_double(n_max);

/*
 *	Arrangement (a = 1, b = 2, c = 3), atomic masses and PES:
 */

	const char arrang = 96 + read_int_keyword(stdin, "arrang", 1, 3, 1);

	pes_init_mass(stdin, 'a');
	pes_init_mass(stdin, 'b');
	pes_init_mass(stdin, 'c');

	pes_init();

	double (*pec)(const size_t, const double), mass = 0.0;

	switch (arrang)
	{
		case 'a':
			pec = pes_bc;
			mass = pes_mass_bc();
		break;

		case 'b':
			pec = pes_ac;
			mass = pes_mass_ac();
		break;

		case 'c':
			pec = pes_ab;
			mass = pes_mass_ab();
		break;

		default:
			PRINT_ERROR("invalid arrangement %c\n", arrang)
			exit(EXIT_FAILURE);
	}

	ASSERT(mass != 0.0)

/*
 *	Directory to store all basis functions:
 */

	char *dir = read_str_keyword(stdin, "basis_dir", ".");

/*
 *	Resolve the diatomic eigenvalue for each j-case and sort results as scatt. channels:
 */

	printf("# Reduced mass = %f a.u., basis dir. = %s, PES name = %s\n", mass, dir, pes_name());
	printf("#     J      ch.      v       j       l       p        E (a.u.)       E (cm-1)        E (eV)   \n");
	printf("# ---------------------------------------------------------------------------------------------\n");

	size_t *ch_counter = allocate(J_max + 1, sizeof(size_t), true);

	fgh_basis basis =
	{
		.v = 0,
		.j = 0,
		.l = 0,
		.n = 0,
		.r_min = r_min,
		.r_max = r_max,
		.r_step = r_step,
		.grid_size = n_max,
		.eigenval = 0.0,
		.eigenvec = NULL
	};

	for (basis.j = j_min; basis.j <= j_max; basis.j += j_step)
	{
		double *pot_energy = allocate(n_max, sizeof(double), false);

		for (size_t n = 0; n < n_max; ++n)
			pot_energy[n] = pec(basis.j, r_min + as_double(n)*r_step);

		matrix *fgh = fgh_dense_single_channel(n_max, r_step, pot_energy, mass);

		free(pot_energy);

		double *eigenval = matrix_symm_eigen(fgh, 'v');

		for (basis.v = v_min; basis.v <= v_max; basis.v += v_step)
		{
			basis.eigenval = eigenval[basis.v];
			basis.eigenvec = fgh_eigenvec(fgh, basis.v, r_step);

			for (size_t J = J_min; J <= J_max; J += J_step)
			{
				for (basis.l = abs(J - basis.j); basis.l <= (J + basis.j); ++basis.l)
				{
					if (parity(basis.j + basis.l) != J_parity && J_parity != 0) continue;

					print_level(&basis, ch_counter[J], J);

					fgh_basis_save(&basis, dir, arrang, ch_counter[J], J);

					ch_counter[J] += 1;
				}
			}

			free(basis.eigenvec);
		}

		matrix_free(fgh);
		free(eigenval);
	}

	free(dir);
	free(ch_counter);

	matrix_end_gpu();

	return EXIT_SUCCESS;
}
