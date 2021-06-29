#include "modules/pes.h"
#include "modules/fgh.h"
#include "modules/file.h"
#include "modules/matrix.h"
#include "modules/globals.h"

#if !defined(COUPLING_MATRIX_FILE_FORMAT)
	#define COUPLING_MATRIX_FILE_FORMAT "cmatrix_arrang=%c_n=%zu_J=%zu.bin"
#endif

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

 Function load_cmatrix(): read the coupling matrix from the disk for the n-th
 grid point index, arrangement and total angular momentum J.

******************************************************************************/

inline static matrix *load_cmatrix(const char arrang, const size_t n, const size_t J)
{
	char filename[MAX_LINE_LENGTH];
	sprintf(filename, COUPLING_MATRIX_FILE_FORMAT, arrang, n, J);

	return matrix_load(filename);
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

	const int J_parity = read_int_keyword(stdin, "parity", -1.0, 1.0, 0.0);

/*
 *	Vibrational quantum numbers, v:
 */

	const size_t v_min = 	read_int_keyword(stdin, "v_min", 0, 10000, 0);

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

	size_t n_max = read_int_keyword(stdin, "scatt_grid_size", 1, 1000000, 100);

	double R_min = read_dbl_keyword(stdin, "R_min", 0.0, INF, 0.5);

	double R_max = read_dbl_keyword(stdin, "R_max", R_min, INF, R_min + 100.0);

	const double R_step = (R_max - R_min)/as_double(n_max);

/*
 *	Redefine the vibrational grid if working on a smaller range of pre-computed
 *	coupling matrices:
 */

	const size_t n_min = read_int_keyword(stdin, "n_min", 0, n_max, 0);

	n_max = read_int_keyword(stdin, "n_max", n_min, n_max, n_max);

	R_min = R_min + as_double(n_min)*R_step;

	R_max = R_min + as_double(n_max - 1)*R_step;

	ASSERT(n_max >= v_max + 1)

/*
 *	Arrangement (a = 1, b = 2, c = 3) and atomic masses:
 */

	const char arrang = 96 + read_int_keyword(stdin, "arrang", 1, 3, 1);

	pes_init_mass(stdin, 'a');
	pes_init_mass(stdin, 'b');
	pes_init_mass(stdin, 'c');

	const double mass = pes_mass_abc(arrang);

	ASSERT(mass != 0.0)

/*
 *	Resolve the atom-diatom eigenvalues for each j-case and sort results as scatt. channels:
 */

	printf("# Reduced mass = %f a.u., n = [%zu, %zu]\n", mass, n_min, n_max);
	printf("#     J     Ch.      v       j       l       p    Comp.       E (a.u.)       E (cm-1)         E (eV)   \n");
	printf("# -----------------------------------------------------------------------------------------------------\n");

/*
 *	Step 1: loop over rotational states j of the atom-diatom and solve a multichannel
 *	eigenvalue problem using the Fourier grid Hamiltonian discrete variable
 *	representation (FGH-DVR) method.
 */

	size_t *ch_counter = allocate(J_max + 1, sizeof(size_t), true);

	fgh_basis basis =
	{
		.v = 0,
		.j = 0,
		.l = 0,
		.n = 0,
		.r_min = R_min,
		.r_max = R_max,
		.r_step = R_step,
		.grid_size = n_max - n_min,
		.eigenval = 0.0,
		.eigenvec = NULL
	};

	for (basis.j = j_min; basis.j <= j_max; basis.j += j_step)
	{
		tensor *pot_energy = allocate(n_max, sizeof(tensor), true);

		size_t n_counter = 0;
		for (size_t n = n_min; n < n_max; ++n)
		{
			pot_energy[n_counter].value = load_cmatrix(arrang, n, basis.j);

			if (n_counter > 0)
			{
				ASSERT(matrix_rows(pot_energy[n_counter - 1].value) == matrix_rows(pot_energy[n_counter].value))
			}

			++n_counter;
		}

		ASSERT(n_counter == basis.grid_size)

		const size_t max_state = matrix_rows(pot_energy[0].value);

		matrix *fgh
			= fgh_dense_multi_channel(max_state, n_counter, R_step, pot_energy, mass);

		for (size_t n = 0; n < n_counter; ++n)
			matrix_free(pot_energy[n].value);

		free(pot_energy);

		double *eigenval = matrix_symm_eigen(fgh, 'v');

		for (basis.v = v_min; basis.v <= v_max; basis.v += v_step)
		{
			basis.eigenval = eigenval[basis.v];

			for (basis.n = 0; basis.n < max_state; ++basis.n)
			{
				basis.eigenvec = fgh_multi_channel_eigenvec(fgh, R_step, max_state, basis.v, basis.n);

				for (size_t J = J_min; J <= J_max; J += J_step)
				{
					for (basis.l = abs(J - basis.j); basis.l <= (J + basis.j); ++basis.l)
					{
						if (parity(basis.j + basis.l) != J_parity && J_parity != 0) continue;

						printf(FORMAT, J, ch_counter[J], basis.v, basis.j, basis.l, parity(basis.j + basis.l),
						       basis.n, basis.eigenval, basis.eigenval*219474.63137054, basis.eigenval*27.211385);

						fgh_basis_save(&basis, ".", arrang, ch_counter[J], J);

						ch_counter[J] += 1;
					}
				}

				free(basis.eigenvec);
			}
		}

		matrix_free(fgh);
		free(eigenval);
	}

	free(ch_counter);

	matrix_end_gpu();

	return EXIT_SUCCESS;
}
