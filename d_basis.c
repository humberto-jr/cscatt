#include "modules/pes.h"
#include "modules/fgh.h"
#include "modules/file.h"
#include "modules/matrix.h"
#include "modules/globals.h"

#define FORMAT "# %5d   %5d   %5d   %5d   %5d   %+5d     % -8e  % -8e  % -8e\n"

/******************************************************************************

 Function parity(): return the parity p = (-1)^l for a given l; either p = +1
 or p = -1.

******************************************************************************/

inline static int parity(const int l)
{
	return (int) pow(-1.0, l);
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

	const int J_min = (int) file_keyword(stdin, "J_min", 0.0, INF, 0.0);

	const int J_max = (int) file_keyword(stdin, "J_max", 0.0, INF, 0.0);

	const int J_step = (int) file_keyword(stdin, "J_step", 1.0, INF, 0.0);

	const int J_parity = (int) file_keyword(stdin, "parity", -1.0, 1.0, 0.0);

	ASSERT(J_max >= J_min)

/*
 *	Vibrational quantum numbers, v:
 */

	const int v_min = (int) file_keyword(stdin, "v_min", 0.0, INF, 0.0);

	const int v_max = (int) file_keyword(stdin, "v_max", 0.0, INF, 0.0);

	const int v_step = (int) file_keyword(stdin, "v_step", 1.0, INF, 1.0);

	ASSERT(v_max >= v_min)

/*
 *	Rotational quantum numbers, j:
 */

	const int j_min = (int) file_keyword(stdin, "j_min", 0.0, INF, 0.0);

	const int j_max = (int) file_keyword(stdin, "j_max", 0.0, INF, 0.0);

	const int j_step = (int) file_keyword(stdin, "j_step", 1.0, INF, 1.0);

	ASSERT(j_max >= j_min)

/*
 *	Vibrational grid:
 */

	const int n_max = (int) file_keyword(stdin, "rovib_grid_size", 1.0, INF, 500.0);

	const double r_min = file_keyword(stdin, "r_min", 0.0, INF, 0.5);

	const double r_max = file_keyword(stdin, "r_max", r_min, INF, r_min + 30.0);

	const double r_step = (r_max - r_min)/as_double(n_max);

	ASSERT(n_max >= v_max + 1)

/*
 *	Arrangement (a = 1, b = 2, c = 3), atomic masses and PES:
 */

	const char arrang = 96 + (int) file_keyword(stdin, "arrang", 1.0, 3.0, 1.0);

	pes_init_mass(stdin, 'a');
	pes_init_mass(stdin, 'b');
	pes_init_mass(stdin, 'c');

	pes_init();

	double (*pec)(const int, const double), mass = 0.0;

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
 *	Resolve the diatomic eigenvalue for each j-case and sort results as scatt. channels:
 */

	printf("# Reduced mass = %f a.u.\n", mass);
	printf("#     J     Ch.       v       j       l       p        E (a.u.)       E (cm-1)        E (eV)   \n");
	printf("# ---------------------------------------------------------------------------------------------\n");

	int *ch_counter = allocate(J_max + 1, sizeof(int), true);

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

/*
 *	Step 1: loop over rotational states j of the diatom and solve the diatomic eigenvalue problem
 *	using the Fourier grid Hamiltonian discrete variable representation (FGH-DVR) method.
 */

	for (basis.j = j_min; basis.j <= j_max; basis.j += j_step)
	{
		double *pot_energy = allocate(n_max, sizeof(double), false);

		for (int n = 0; n < n_max; ++n)
			pot_energy[n] = pec(basis.j, r_min + as_double(n)*r_step);

		matrix *fgh
			= fgh_dense_single_channel(n_max, r_step, pot_energy, mass);

		free(pot_energy);

		double *eigenval = matrix_symm_eigen(fgh, 'v');

/*
 *		Step 2: loop over the vibrational states v of the diatom, solutions of step 2, and select
 *		only those of interest.
 */

		for (basis.v = v_min; basis.v <= v_max; basis.v += v_step)
		{
			basis.eigenval = eigenval[basis.v];
			basis.eigenvec = fgh_eigenvec(fgh, basis.v, r_step);

/*
 *			Step 3: loop over all partial waves l of the atom around the diatom given by the
 *			respective J and j.
 */

			for (int J = J_min; J <= J_max; J += J_step)
			{
				for (basis.l = abs(J - basis.j); basis.l <= (J + basis.j); ++basis.l)
				{
					if (parity(basis.j + basis.l) != J_parity && J_parity != 0) continue;

					printf(FORMAT, J, ch_counter[J], basis.v, basis.j, basis.l, parity(basis.j + basis.l),
					       basis.eigenval, basis.eigenval*219474.63137054, basis.eigenval*27.211385);

/*
 *					Step 4: save each basis function |vjl> in the disk and increment the counter of
 *					channels.
 */

					FILE *output = fgh_basis_file(arrang, ch_counter[J], J, "wb", false);

					fgh_basis_write(&basis, output);

					file_close(&output);

					ch_counter[J] += 1;
				}
			}

			free(basis.eigenvec);
		}

		matrix_free(fgh);
		free(eigenval);
	}

	free(ch_counter);

	matrix_end_gpu();

	return EXIT_SUCCESS;
}
