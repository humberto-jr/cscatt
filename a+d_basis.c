#include "modules/pes.h"
#include "modules/file.h"
#include "modules/matrix.h"
#include "modules/mpi_lib.h"
#include "modules/globals.h"

#include "utils.h"

#define FORMAT "# %5d   %5d   %5d   %5d   %5d   %+5d   %5d     % -8e  % -8e  % -8e\n"

/******************************************************************************

 Function parity(): return the parity p = (-1)^l for a given l; either p = +1
 or p = -1.

******************************************************************************/

inline static int parity(const int l)
{
	return (int) pow(-1.0, l);
}

/******************************************************************************

 Function fgh_multich_matrix(): the same as fgh_matrix(), except that the
 potential energy is the one of a problem with max_ch channels within grid_size
 points. For details see Eq. (19a), Eq. (19b) and Eq. (19c) of Ref. [2].

******************************************************************************/

mpi_matrix *fgh_multich_matrix(const int max_ch,
                               const int grid_size,
                               const double grid_step,
                               const tensor pot_energy[],
                               const double mass)
{
	ASSERT(max_ch > 0)
	ASSERT(pot_energy != NULL)

	const int size = grid_size*max_ch;

	const int non_zeros[2] = {1, 1};

	mpi_matrix *result = mpi_matrix_alloc(size, size, non_zeros);

	const double box_length = as_double(grid_size - 1)*grid_step;

	const double factor = (M_PI*M_PI)/(mass*box_length*box_length);

	const double nn_term = factor*as_double(grid_size*grid_size + 2)/6.0;

	int row = 0;
	for (int p = 0; p < max_ch; ++p)
	{
		for (int n = 0; n < grid_size; ++n)
		{
			int col = 0;
			for (int q = 0; q < max_ch; ++q)
			{
				for (int m = 0; m < grid_size; ++m)
				{
					double pnqm = 0.0;

					if (n == m && p == q)
					{
						pnqm = nn_term
						     + matrix_get(pot_energy[n].value, p, p);
					}

					if (n != m && p == q)
					{
						const double nm = as_double((n + 1) - (m + 1));

						const double nm_term = sin(nm*M_PI/as_double(grid_size));

						pnqm = pow(-1.0, nm)*factor/pow(nm_term, 2);
					}

					if (n == m && p != q)
					{
						pnqm = matrix_get(pot_energy[n].value, p, q);
					}

					if (pnqm != 0.0) mpi_matrix_set(result, row, col, pnqm);
					++col;
				}
			}

			++row;
		}
	}

	mpi_barrier();
	mpi_matrix_build(result);

	return result;
}

int main(int argc, char *argv[])
{
	mpi_init(argc, argv);

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

	int n_max = (int) file_keyword(stdin, "scatt_grid_size", 1.0, INF, 500.0);

	double R_min = file_keyword(stdin, "R_min", 0.0, INF, 0.5);

	double R_max = file_keyword(stdin, "R_max", R_min, INF, R_min + 100.0);

	const double R_step = (R_max - R_min)/as_double(n_max);

/*
 *	Redefine the vibrational grid if working on a smaller range of pre-computed
 *	data:
 */

	const int n_min = (int) file_keyword(stdin, "n_min", 0.0, as_double(n_max), 0.0);

	n_max = (int) file_keyword(stdin, "n_max", as_double(n_min), as_double(n_max), as_double(n_max));

	R_min = R_min + as_double(n_min)*R_step;

	R_max = R_min + as_double(n_max - 1)*R_step;

	ASSERT(n_max >= v_max + 1)

/*
 *	Arrangement (a = 1, b = 2, c = 3), atomic masses:
 */

	const char arrang = 96 + (int) file_keyword(stdin, "arrang", 1.0, 3.0, 1.0);

	pes_init_mass(stdin, 'a');
	pes_init_mass(stdin, 'b');
	pes_init_mass(stdin, 'c');

	const double mass = pes_mass_abc(arrang);

	ASSERT(mass != 0.0)

/*
 *	Number of iterations used in the eigensolver:
 */

	const int max_step = (int) file_keyword(stdin, "max_step", 1.0, INF, 300.0);

	const double tol = file_keyword(stdin, "tolerance", -INF, INF, 1.0E-8);

/*
 *	Resolve the triatomic eigenvalues for each j-case and sort results as scatt. channels:
 */

	if (mpi_rank() == 0)
	{
		printf("# Reduced mass = %f a.u., n = [%d, %d]\n", mass, n_min, n_max);
		printf("#     J      Ch.      v       j       l       p    Comp.      E (a.u.)       E (cm-1)         E (eV)   \n");
		printf("# -----------------------------------------------------------------------------------------------------\n");
	}

/*
 *	Step 1: loop over rotational states j of the triatom and solve a multichannel
 *	eigenvalue problem using the Fourier grid Hamiltonian discrete variable
 *	representation (FGH-DVR) method.
 */

	int *ch_counter = allocate(J_max + 1, sizeof(int), true);

	for (int j = j_min; j <= j_max; j += j_step)
	{
		tensor *pot_energy = allocate(n_max, sizeof(tensor), true);

		int n_counter = 0;
		for (int n = n_min; n < n_max; ++n)
		{
			pot_energy[n_counter].value = coupling_read(arrang, n, false, j);

			if (n_counter > 0)
			{
				ASSERT(matrix_row(pot_energy[n_counter - 1].value) == matrix_row(pot_energy[n_counter].value))
			}

			++n_counter;
		}

/*
 *		NOTE: the total number of diatomic states used to expand the triatomic
 *		eigenvectos is named max_state and it is equal ch_counter from dbasis driver.
 */

		const int max_state = matrix_row(pot_energy[0].value);

		mpi_matrix *fgh
			= fgh_multich_matrix(max_state, n_counter, R_step, pot_energy, mass);

		for (int n = 0; n < n_counter; ++n)
			matrix_free(pot_energy[n].value);

		free(pot_energy);

		const int info
			= mpi_matrix_sparse_eigen(fgh, v_max + 1, max_step, tol, false);

		if (info < (v_max + 1))
		{
			PRINT_ERROR("only %d/%d solutions converged for j = %d\n", info, v_max + 1, j)
			exit(EXIT_FAILURE);
		}

/*
 *		Step 2: loop over the vibrational states v of the triatom, solutions of
 *		step 1, and select only those of interest.
 */

		for (int v = v_min; v <= v_max; v += v_step)
		{
			double eigenval = 0.0;
			mpi_vector *eigenvec = mpi_matrix_eigenpair(fgh, v, &eigenval);
/*
 *			Step 3: loop over all partial waves l of the atom around the triatom
 *			given by the respective J and j.
 */
			for (int J = J_min; J <= J_max; J += J_step)
			{
				for (int l = abs(J - j); l <= (J + j); ++l)
				{
					if (parity(j + l) != J_parity && J_parity != 0) continue;

					for (int n = 0; n < max_state; ++n)
					{
/*
 *						Step 4: save each basis function |vjln> in the disk and increment
 *						the counter of atom-triatom channels.
 */
						FILE *output = NULL;

						if (eigenval != 0.0)
						{
							/* NOTE: only one CPU shall have the v-th non-zero eigenvalue. */

							printf(FORMAT, J, ch_counter[J], v, j, l, parity(j + l), n,
							       eigenval, eigenval*219474.63137054, eigenval*27.211385);
						}

						if (mpi_rank() == 0)
						{
							output = basis_file(arrang, ch_counter[J], J, "wb", false);

							file_write(&v, sizeof(int), 1, output);
							file_write(&j, sizeof(int), 1, output);
							file_write(&l, sizeof(int), 1, output);
							file_write(&n, sizeof(int), 1, output);

							file_write(&R_min, sizeof(double), 1, output);
							file_write(&R_max, sizeof(double), 1, output);
							file_write(&R_step, sizeof(double), 1, output);

							file_write(&eigenval, sizeof(double), 1, output);

							file_write(&n_counter, sizeof(int), 1, output);
						}

						const int start = n*n_counter;
						const int end = start + n_counter;

						/* NOTE: all pieces of the v-th eigenvector, per CPU, should collapse on the same file written by CPU 0. */
						mpi_vector_write(eigenvec, start, end, output);

						if (mpi_rank() == 0) file_close(&output);

						ch_counter[J] += 1;
					}
				}
			}

			mpi_vector_free(eigenvec);
		}

		mpi_matrix_free(fgh);
	}

	free(ch_counter);

	mpi_end();
	return EXIT_SUCCESS;
}
