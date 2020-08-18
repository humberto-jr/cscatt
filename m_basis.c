#include "modules/pes.h"
#include "modules/file.h"
#include "modules/matrix.h"
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

 Function fgh_matrix_product(): the same as fgh_matrix(), except that the
 potential energy is the one of a problem with max_ch channels within grid_size
 points. For details see Eq. (19a), Eq. (19b) and Eq. (19c) of Ref. [2].

******************************************************************************/

struct fgh_params
{
	const tensor *pot_energy;
	const int grid_size, max_ch;
	const double grid_step, mass;
};

void fgh_matrix_product(const int size,
                        const int from,
                        const int to,
                        double vector[],
                        void *params)
{
	ASSERT(vector != NULL)
	ASSERT(params != NULL)
	ASSERT(size == fgh->max_ch*fgh->grid_size)

	struct fgh_params *fgh = (struct fgh_params *) params;

	const double box_length = as_double(fgh->grid_size - 1)*fgh->grid_step;

	const double factor = (M_PI*M_PI)/(fgh->mass*box_length*box_length);

	const double nn_term = factor*as_double(fgh->grid_size*fgh->grid_size + 2)/6.0;

	int row = 0;
	for (int p = 0; p < fgh->max_ch; ++p)
	{
		for (int n = 0; n < fgh->grid_size; ++n)
		{
			double sum = 0.0;

			int col = 0;
			for (int q = 0; q < fgh->max_ch; ++q)
			{
				for (int m = 0; m < fgh->grid_size; ++m)
				{
					double pnqm = 0.0;

					if (n == m && p == q)
					{
						pnqm = nn_term + matrix_get(fgh->pot_energy[n].value, p, p);
					}

					if (n != m && p == q)
					{
						const double nm = as_double((n + 1) - (m + 1));

						const double nm_term = sin(nm*M_PI/as_double(fgh->grid_size));

						pnqm = pow(-1.0, nm)*factor/pow(nm_term, 2);
					}

					if (n == m && p != q)
					{
						pnqm = matrix_get(fgh->pot_energy[n].value, p, q);
					}

					sum += *(vector + from + col)*pnqm;
					++col;
				}
			}

			*(vector + to + row) = sum;
			++row;
		}
	}
}

/******************************************************************************

 Function fgh_multich_matrix(): the same as fgh_matrix(), except that the
 potential energy is the one of a problem with max_ch channels within grid_size
 points. For details see Eq. (19a), Eq. (19b) and Eq. (19c) of Ref. [2].

******************************************************************************/

matrix *fgh_multich_matrix(const int max_ch,
                           const int grid_size,
                           const double grid_step,
                           const tensor pot_energy[],
                           const double mass)
{
	ASSERT(max_ch > 0)
	ASSERT(pot_energy != NULL)

	matrix *result = matrix_alloc(grid_size*max_ch, grid_size*max_ch, false);

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
						pnqm = nn_term + matrix_get(pot_energy[n].value, p, p);
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

					matrix_set(result, row, col, pnqm);
					++col;
				}
			}

			++row;
		}
	}

	return result;
}

/******************************************************************************

 Function norm(): normalize to unity a radial multichannel eigenvector composed
 of max_state components and grid_size points, using a 1/3-Simpson quadrature
 rule.

******************************************************************************/

void norm(const int max_state, const int grid_size,
          const double grid_step, const bool use_omp, double eigenvec[])
{
	ASSERT(max_state > 0)
	ASSERT(grid_size > 0)
	ASSERT(eigenvec != NULL)

	double total_sum = 0.0;

	#pragma omp parallel for default(none) shared(eigenvec) reduction(+:total_sum) if(use_omp)
	for (int m = 0; m < max_state; ++m)
	{
		const int n_min = m*grid_size;
		const int n_max = n_min + (grid_size%2 == 0? grid_size : grid_size - 1);

		double sum
			= eigenvec[n_min]*eigenvec[n_min] + eigenvec[n_max - 1]*eigenvec[n_max - 1];

		for (int n = n_min + 1; n < (n_max - 2); n += 2)
		{
			sum += 4.0*eigenvec[n]*eigenvec[n];
			sum += 2.0*eigenvec[n + 1]*eigenvec[n + 1];
		}

		sum = grid_step*sum/3.0;
		total_sum += sum;
	}

	total_sum = sqrt(total_sum);

	#pragma omp parallel for default(none) shared(eigenvec, total_sum) if(use_omp)
	for (int n = 0; n < max_state*grid_size; ++n)
	{
		eigenvec[n] = eigenvec[n]/total_sum;
	}
}

int main(int argc, char *argv[])
{
	ASSERT(argc > 1)
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
 *	Resolve the triatomic eigenvalues for each j-case and sort results as scatt. channels:
 */

	printf("# Reduced mass = %f a.u., n = [%d, %d]\n", mass, n_min, n_max);
	printf("#     J     Ch.      v       j       l       p    Comp.       E (a.u.)       E (cm-1)         E (eV)   \n");
	printf("# -----------------------------------------------------------------------------------------------------\n");

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

		matrix *eigenvec
			= fgh_multich_matrix(max_state, n_counter, R_step, pot_energy, mass);

		for (int n = 0; n < n_counter; ++n)
		{
			matrix_free(pot_energy[n].value);
		}

		free(pot_energy);

		double *eigenval = matrix_symm_eigen(eigenvec, 'v');

/*
 *		Step 2: loop over the vibrational states v of the triatom, solutions of
 *		step 1, and select only those of interest.
 */

		for (int v = v_min; v <= v_max; v += v_step)
		{
			double *wavef = matrix_raw_col(eigenvec, v, false);

			norm(max_state, n_counter, R_step, (n_counter >= 2000), wavef);

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
						printf(FORMAT, J, ch_counter[J], v, j, l, parity(j + l), n,
						       eigenval[v], eigenval[v]*219474.63137054, eigenval[v]*27.211385);

/*
 *						Step 4: save each basis function |vjln> in the disk and increment
 *						the counter of atom-triatom channels.
 */

						FILE *output = basis_file(arrang, ch_counter[J], J, "wb");

						file_write(&v, sizeof(int), 1, output);
						file_write(&j, sizeof(int), 1, output);
						file_write(&l, sizeof(int), 1, output);
						file_write(&n, sizeof(int), 1, output);

						file_write(&R_min, sizeof(double), 1, output);
						file_write(&R_max, sizeof(double), 1, output);
						file_write(&R_step, sizeof(double), 1, output);

						file_write(&eigenval[v], sizeof(double), 1, output);

						file_write(&n_counter, sizeof(int), 1, output);
						file_write(wavef + n*n_counter, sizeof(double), n_counter, output);

						file_close(&output);
						ch_counter[J] += 1;
					}
				}
			}

			free(wavef);
		}

		matrix_free(eigenvec);
		free(eigenval);
	}

	free(ch_counter);

	return EXIT_SUCCESS;
}
