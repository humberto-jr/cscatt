#include "modules/pes.h"
#include "modules/fgh.h"
#include "modules/file.h"
#include "modules/matrix.h"
#include "modules/mpi_lib.h"
#include "modules/globals.h"

#define FORMAT "# %5d   %5d   %5d   %5d   %5d   %+5d   %5d     % -8e  % -8e  % -8e\n"

#if !defined(COUPLING_MATRIX_FILE_FORMAT)
	#define COUPLING_MATRIX_FILE_FORMAT "cmatrix_arrang=%c_n=%d_J=%d.bin"
#endif

/******************************************************************************

 Function parity(): return the parity p = (-1)^l for a given l; either p = +1
 or p = -1.

******************************************************************************/

inline static int parity(const int l)
{
	return (int) pow(-1.0, l);
}

/******************************************************************************

 Function read_coupling(): load from the disk a set of multipole terms for the
 n-th grid point index, arrangement and total angular momentum, J.

******************************************************************************/

matrix *read_coupling(const char arrang, const int n, const int J)
{
	char filename[MAX_LINE_LENGTH];
	sprintf(filename, COUPLING_MATRIX_FILE_FORMAT, arrang, n, J);

	return matrix_load(filename);
}

/******************************************************************************
******************************************************************************/

int main(int argc, char *argv[])
{
	mpi_init(argc, argv);

	file_init_stdin(argv[1]);

/*
 *	Total angular momentum, J:
 */

	const int J_min = read_int_keyword(stdin, "J_min", 0, 10000, 0);

	const int J_max = read_int_keyword(stdin, "J_max", J_min, 10000, J_min);

	const int J_step = read_int_keyword(stdin, "J_step", 1, 10000, 1);

	const int J_parity = read_int_keyword(stdin, "parity", -1, 1, 0);

/*
 *	Vibrational quantum numbers, v:
 */

	const int v_min = read_int_keyword(stdin, "v_min", 0, 10000, 0);

	const int v_max = read_int_keyword(stdin, "v_max", v_min, 10000, v_min);

	const int v_step = read_int_keyword(stdin, "v_step", 1, 10000, 1);

/*
 *	Rotational quantum numbers, j:
 */

	const int j_min = read_int_keyword(stdin, "j_min", 0, 10000, 0);

	const int j_max = read_int_keyword(stdin, "j_max", j_min, 10000, j_min);

	const int j_step = read_int_keyword(stdin, "j_step", 1, 10000, 1);

/*
 *	Vibrational grid, R:
 */

	int n_max = read_int_keyword(stdin, "scatt_grid_size", 1, 1000000, 500);

	double R_min = read_dbl_keyword(stdin, "R_min", 0.0, INF, 0.5);

	double R_max = read_dbl_keyword(stdin, "R_max", R_min, INF, R_min + 30.0);

	const double R_step = (R_max - R_min)/as_double(n_max);

/*
 *	Redefine the vibrational grid if working on a smaller range of pre-computed
 *	data, i.e. n_min > 0 and n_max < scatt_grid_size.
 */

	const int n_min = read_int_keyword(stdin, "n_min", 0, n_max, 0);

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
 *	Number of iterations and tolerance (tol) used by the eigensolver:
 */

	const int max_step = read_int_keyword(stdin, "max_step", 1, 1000000, 20000);

	const double tol = read_dbl_keyword(stdin, "tolerance", -INF, INF, 1.0E-6);

/*
 *	Resolve the triatomic eigenvalues for each j-case and sort results as scatt. channels:
 */

	if (mpi_rank() == 0)
	{
		printf("# Reduced mass = %f a.u., n = [%d, %d]\n", mass, n_min, n_max);
		printf("#     J      Ch.      v       j       l       p    Comp.      E (a.u.)       E (cm-1)         E (eV)   \n");
		printf("# -----------------------------------------------------------------------------------------------------\n");
	}

	int *ch_counter = allocate(J_max + 1, sizeof(int), true);

	for (int j = j_min; j <= j_max; j += j_step)
	{
		tensor *pot_energy = allocate(n_max, sizeof(tensor), true);

		int n_counter = 0;
		for (int n = n_min; n < n_max; ++n)
		{
			pot_energy[n_counter].value = read_coupling(arrang, n, j);

			if (n_counter > 0)
			{
				ASSERT(matrix_rows(pot_energy[n_counter - 1].value) == matrix_rows(pot_energy[n_counter].value))
			}

			++n_counter;
		}

/*
 *		NOTE: the total number of diatomic states used to expand the triatomic
 *		eigenvectos is named max_state and it is equal ch_counter from dbasis driver.
 */

		const int max_state
			= matrix_rows(pot_energy[0].value);

		mpi_matrix *fgh
			= fgh_sparse_multi_channel(max_state, n_counter, R_step, pot_energy, mass);

		for (int n = 0; n < n_counter; ++n)
			matrix_free(pot_energy[n].value);

		free(pot_energy);

		const int info
			= mpi_matrix_sparse_eigen(fgh, v_max + 1, max_step, tol, false);

		if (info < (v_max + 1) && mpi_rank() == 0)
		{
			PRINT_ERROR("only %d/%d solutions converged for j = %d\n", info, v_max + 1, j)
			exit(EXIT_FAILURE);
		}

		for (int v = v_min; v <= v_max; v += v_step)
		{
			double eigenval = 0.0;
			mpi_vector *eigenvec = mpi_matrix_eigenpair(fgh, v, &eigenval);

			for (int J = J_min; J <= J_max; J += J_step)
			{
				for (int l = abs(J - j); l <= (J + j); ++l)
				{
					if (parity(j + l) != J_parity && J_parity != 0) continue;

					for (int n = 0; n < max_state; ++n)
					{
						FILE *output = NULL;

						if (eigenval != 0.0)
						{
							/* NOTE: only one CPU shall have the v-th non-zero eigenvalue. */

							printf(FORMAT, J, ch_counter[J], v, j, l, parity(j + l), n,
							       eigenval, eigenval*219474.63137054, eigenval*27.211385);
						}

						if (mpi_rank() == 0)
						{
							output = fgh_basis_file(".", arrang, ch_counter[J], J, "wb", false);

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

						/* NOTE: all pieces of the v-th eigenvector, per CPU, should collapse into the same file written by CPU 0. */
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
