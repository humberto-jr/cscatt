/******************************************************************************

 About
 -----

 This module defines wrappers for popular implementations of the BLAS library,
 including those using GPUs, and selectable during compilation.

******************************************************************************/

#include "fgh.h"
#include "mpi_lb.h"
#include "matrix.h"
#include "globals.h"

/******************************************************************************

 Function fgh_matrix(): return the discrete variable representation (DVR) of a
 Fourier grid Hamiltonian (FGH) subjected to a given single channel potential
 energy and reduced mass. As shown in Eq. (6a) and (6b) of Ref. [1].

 NOTE: The m-th eigenvector is the m-th column in a grid_size-by-grid_size row-
 major matrix: eigenvec[n*grid_size + m], where n = m = [0, grid_size).

******************************************************************************/

matrix *fgh_dense_single_ch(const int grid_size, const double grid_step,
                            const double pot_energy[], const double mass)
{
	ASSERT(pot_energy != NULL)

	matrix *result = matrix_alloc(grid_size, grid_size, false);

	const double box_length = as_double(grid_size - 1)*grid_step;

	const double factor = (M_PI*M_PI)/(mass*box_length*box_length);

	const double nn_term = factor*as_double(grid_size*grid_size + 2)/6.0;

	for (int n = 0; n < grid_size; ++n)
	{
		matrix_diag_set(result, n, nn_term + pot_energy[n]);

		for (int m = (n + 1); m < grid_size; ++m)
		{
			const double nm = as_double((n + 1) - (m + 1));

			const double nm_term = sin(nm*M_PI/as_double(grid_size));

			matrix_symm_set(result, n, m, pow(-1.0, nm)*factor/pow(nm_term, 2));
		}
	}

	return result;
}

/******************************************************************************

 Function fgh_multich_matrix(): the same as fgh_matrix(), except that the
 potential energy is the one of a problem with max_ch channels within grid_size
 points. For details see Eq. (19a), Eq. (19b) and Eq. (19c) of Ref. [2].

******************************************************************************/

mpi_matrix *fgh_sparse_multi_ch(const int max_ch,
                                const int grid_size,
                                const double grid_step,
                                const tensor pot_energy[],
                                const double mass)
{
	ASSERT(max_ch > 0)
	ASSERT(pot_energy != NULL)

	const int size = grid_size*max_ch;

	const int non_zeros[2] = {grid_size + max_ch, grid_size + max_ch};

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

/******************************************************************************

 Function fgh_norm(): normalize to unity a radial eigenvector of grid_size
 points, using a 1/3-Simpson quadrature rule.

******************************************************************************/

void fgh_norm(const int grid_size,
              const double grid_step, const bool use_omp, double eigenvec[])
{
	ASSERT(grid_size > 0)
	ASSERT(eigenvec != NULL)

	const int n_max
		= (grid_size%2 == 0? grid_size : grid_size - 1);

	double sum
		= eigenvec[0]*eigenvec[0] + eigenvec[n_max - 1]*eigenvec[n_max - 1];

	#pragma omp parallel for default(none) shared(eigenvec) reduction(+:sum) if(use_omp)
	for (int n = 1; n < (n_max - 2); n += 2)
	{
		sum += 4.0*eigenvec[n]*eigenvec[n];
		sum += 2.0*eigenvec[n + 1]*eigenvec[n + 1];
	}

	sum = grid_step*sum/3.0;
	sum = sqrt(sum);

	#pragma omp parallel for default(none) shared(eigenvec, sum) if(use_omp)
	for (int n = 0; n < grid_size; ++n)
	{
		eigenvec[n] = eigenvec[n]/sum;
	}
}

/******************************************************************************

 Function basis_read(): loads from the disk the basis function for the n-th
 channel for a given arrangement and total angular momentum, J.

******************************************************************************/

void fgh_basis_load(const char arrang, const int n, const int J, fgh_basis *b)
{
	ASSERT(b != NULL)

	FILE *input = basis_file(arrang, n, J, "rb");

	file_read(&b->v, sizeof(int), 1, input, 0);

	file_read(&b->j, sizeof(int), 1, input, 0);

	file_read(&b->l, sizeof(int), 1, input, 0);

	file_read(&b->n, sizeof(int), 1, input, 0);

	file_read(&b->r_min, sizeof(double), 1, input, 0);

	file_read(&b->r_max, sizeof(double), 1, input, 0);

	file_read(&b->r_step, sizeof(double), 1, input, 0);

	file_read(&b->eigenval, sizeof(double), 1, input, 0);

	file_read(&b->grid_size, sizeof(int), 1, input, 0);

	b->eigenvec = allocate(b->grid_size, sizeof(double), false);

	file_read(b->eigenvec, sizeof(double), b->grid_size, input, 0);

	fclose(input);
}
