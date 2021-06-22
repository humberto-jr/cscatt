/******************************************************************************

 About
 -----

 This module defines wrappers for popular implementations of the BLAS library,
 including those using GPUs, and selectable during compilation.

******************************************************************************/

#include "globals.h"
#include "mpi_lib.h"
#include "matrix.h"
#include "file.h"
#include "math.h"
#include "fgh.h"

#if !defined(FGH_BASIS_FORMAT)
	#define FGH_BASIS_FORMAT "%s/basis_arrang=%c_ch=%zu_J=%zu.%s"
#endif

/******************************************************************************

 Function fgh_matrix(): return the discrete variable representation (DVR) of a
 Fourier grid Hamiltonian (FGH) subjected to a given single channel potential
 energy and reduced mass. As shown in Eq. (6a) and (6b) of Ref. [1].

 NOTE: The m-th eigenvector is the m-th column in a grid_size-by-grid_size row-
 major matrix: eigenvec[n*grid_size + m], where n = m = [0, grid_size).

******************************************************************************/

matrix *fgh_dense_single_channel(const size_t grid_size,
                                 const double grid_step,
                                 const double pot_energy[], const double mass)
{
	ASSERT(pot_energy != NULL)

	matrix *result = matrix_alloc(grid_size, grid_size, false);

	const double box_length = as_double(grid_size - 1)*grid_step;

	const double factor = (M_PI*M_PI)/(mass*box_length*box_length);

	const double nn_term = factor*as_double(grid_size*grid_size + 2)/6.0;

	for (size_t n = 0; n < grid_size; ++n)
	{
		matrix_set_diag(result, n, nn_term + pot_energy[n]);

		for (size_t m = (n + 1); m < grid_size; ++m)
		{
			const double nm = as_double((n + 1) - (m + 1));

			const double nm_term = sin(nm*M_PI/as_double(grid_size));

			matrix_set_symm(result, n, m, pow(-1.0, nm)*factor/pow(nm_term, 2));
		}
	}

	return result;
}

/******************************************************************************

 Function dvr_multich_fgh(): the same as dvr_fgh(), except that the potential
 energy is the one of a problem with max_ch channels within grid_size points.
 For details see Eq. (19a), Eq. (19b) and Eq. (19c) of Ref. [2].

******************************************************************************/

matrix *fgh_dense_multi_channel(const size_t max_state,
                                const size_t grid_size,
                                const double grid_step,
                                const tensor pot_energy[], const double mass)
{
	ASSERT(max_state > 0)
	ASSERT(pot_energy != NULL)

	matrix *result = matrix_alloc(grid_size*max_state, grid_size*max_state, false);

	const double box_length = as_double(grid_size - 1)*grid_step;

	const double factor = (M_PI*M_PI)/(mass*box_length*box_length);

	const double nn_term = factor*as_double(grid_size*grid_size + 2)/6.0;

	size_t row = 0;
	for (size_t p = 0; p < max_state; ++p)
	{
		for (size_t n = 0; n < grid_size; ++n)
		{
			size_t col = 0;
			for (size_t q = 0; q < max_state; ++q)
			{
				for (size_t m = 0; m < grid_size; ++m)
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

 Function fgh_multich_matrix(): the same as fgh_matrix(), except that the
 potential energy is the one of a problem with max_ch channels within grid_size
 points. For details see Eq. (19a), Eq. (19b) and Eq. (19c) of Ref. [2].

******************************************************************************/

mpi_matrix *fgh_sparse_multi_channel(const size_t max_state,
                                     const size_t grid_size,
                                     const double grid_step,
                                     const tensor pot_energy[], const double mass)
{
	ASSERT(max_state > 0)
	ASSERT(pot_energy != NULL)

	const size_t size = grid_size*max_state;

	const int non_zeros[2] = {1, 1};

	mpi_matrix *result = mpi_matrix_alloc(size, size, non_zeros);

	const double box_length = as_double(grid_size - 1)*grid_step;

	const double factor = (M_PI*M_PI)/(mass*box_length*box_length);

	const double nn_term = factor*as_double(grid_size*grid_size + 2)/6.0;

	size_t row = 0;
	for (size_t p = 0; p < max_state; ++p)
	{
		for (size_t n = 0; n < grid_size; ++n)
		{
			size_t col = 0;
			for (size_t q = 0; q < max_state; ++q)
			{
				for (size_t m = 0; m < grid_size; ++m)
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

 Function dvr_multich_fgh_comp(): return the n-th component of a multichannel
 eigenvector v from the Hamiltonian built by dvr_multich_fgh(). On entry, the
 matrix is expected to be properly diagonalized and its columns the respective
 eigenvectors.

******************************************************************************/
/*
matrix *dvr_multich_fgh_comp(matrix *fgh,
                             const int max_ch, const int n, const int v)
{
	ASSERT(n < max_ch)
	const int grid_size = matrix_rows(fgh)/max_ch;

	return matrix_get_block(fgh, n*grid_size, n*grid_size + (grid_size - 1), v, v);
}*/

/******************************************************************************

 Function fgh_interpolation(): interpolate the eigenvector of the Hamiltonian
 built by fgh_dense_single_channel(), using Eq. (4.1) of Ref. [3], and return
 its amplitude value at a given new r.

******************************************************************************/

double fgh_interpolation(const size_t grid_size,
                         const double eigenvec[],
                         const double r_min,
                         const double r_max,
                         const double r_new)
{
	ASSERT(eigenvec != NULL)

	const double r_step = (r_max - r_min)/as_double(grid_size);

	double result = 0.0;
	for (size_t n = 0; n < grid_size; ++n)
	{
		const double r_old = r_min + as_double(n)*r_step;
		const double param = M_PI*(r_new - r_old)/r_step;

		result += (fabs(param) > 1.0E-7? eigenvec[n]*math_sinc(param) : eigenvec[n]);
	}

	return result;
}

/******************************************************************************

 Function dvr_fgh_product(): the same of dvr_fgh_wavef() but for the product of
 two eigenvectors, a and b.

******************************************************************************/
/*
double dvr_fgh_product(const matrix *fgh, const int a, const int b,
                       const double r_min, const double r_max, const double r_new)
{
	const int n_max = matrix_rows(fgh);

	ASSERT(a > -1)
	ASSERT(a < n_max)

	ASSERT(b > -1)
	ASSERT(b < n_max)

	const double r_step = (r_max - r_min)/as_double(n_max);
	double result_a = 0.0, result_b = 0.0;

	for (int n = 0; n < n_max; ++n)
	{
		const double wavef_a = matrix_get(fgh, n, a);
		const double wavef_b = matrix_get(fgh, n, b);

		const double r_old = r_min + as_double(n)*r_step;
		const double param = M_PI*(r_new - r_old)/r_step;

		if (fabs(param) > 1.0E-7)
		{
			result_a += wavef_a*sinc(param);
			result_b += wavef_b*sinc(param);
		}
		else
		{
			result_a += wavef_a;
			result_b += wavef_b;
		}
	}

	return result_a*result_b;
}*/

/******************************************************************************

 Function fgh_eigenvec(): normalize to unity the eigenvector v of a Hamiltonian
 built by fgh_dense_single_channel(). On entry, the matrix is expected to be
 properly diagonalized with its columns being the respective eigenvectors.

******************************************************************************/

double *fgh_eigenvec(const matrix *fgh, const size_t v, const double grid_step)
{
	ASSERT(fgh != NULL)

	const size_t grid_size = matrix_rows(fgh);
	const size_t n_max = (grid_size%2 == 0? grid_size : grid_size - 1);

	double *eigenvec = NULL;

	if (matrix_using_magma())
		eigenvec = matrix_get_raw_row(fgh, v);
	else
		eigenvec = matrix_get_raw_col(fgh, v);

	ASSERT(eigenvec != NULL)

	double sum
		= eigenvec[0]*eigenvec[0] + eigenvec[n_max - 1]*eigenvec[n_max - 1];

	for (size_t n = 1; n < (n_max - 2); n += 2)
		sum += 4.0*eigenvec[n]*eigenvec[n] + 2.0*eigenvec[n + 1]*eigenvec[n + 1];

	sum = grid_step*sum/3.0;

	const double norm = 1.0/sqrt(sum);

	for (size_t n = 0; n < grid_size; ++n)
		eigenvec[n] = norm*eigenvec[n];

	return eigenvec;
}

/******************************************************************************

 Function dvr_fgh_norm(): the same of dvr_fgh_norm() but for eigenvectos from
 those matrices built by dvr_multich_fgh().

******************************************************************************/

/*
void dvr_multich_fgh_norm(matrix *fgh, const int max_ch,
                          const int v, const double grid_step, const bool use_omp)
{
	const int n_max
	= ((matrix_row(fgh)/max_ch)%2 == 0? matrix_row(fgh)/max_ch : matrix_row(fgh)/max_ch - 1);

	double sum = 0.0;

	#pragma omp parallel for default(none) shared(fgh) reduction(+:sum) if(use_omp)
	for(int m = 0; m < max_ch; ++m)
	{
		matrix *wavef_m = dvr_multich_fgh_comp(fgh, max_ch, m, v);

		sum += matrix_get_pow(wavef_m, 0, 0, 2.0)
		     + matrix_get_pow(wavef_m, n_max - 1, 0, 2.0);

		for (int n = 1; n < (n_max - 2); n += 2)
		{
			sum += 4.0*matrix_get_pow(wavef_m, n, 0, 2.0)
			     + 2.0*matrix_get_pow(wavef_m, n + 1, 0, 2.0);
		}

		matrix_free(wavef_m);
	}

	sum = grid_step*sum/3.0;
	matrix_col_scale(fgh, v, 1.0/sqrt(sum), use_omp);
}*/

/******************************************************************************

 Function fgh_basis_count(): counts how many FGH basis functions are available
 in the disk for a given arrangement with total angular momentum J.

******************************************************************************/

size_t fgh_basis_count(const char dir[], const char arrang, const size_t J)
{
	char filename[MAX_LINE_LENGTH];

	size_t counter = 0;
	sprintf(filename, FGH_BASIS_FORMAT, dir, arrang, counter, J, "bin");

	while (file_exist(filename))
	{
		++counter;
		sprintf(filename, FGH_BASIS_FORMAT, dir, arrang, counter, J, "bin");
	}

	return counter;
}

/******************************************************************************

 Function fgh_basis_file(): opens the file for the n-th channel and a given
 arrangement with total angular momentum J. Where, mode is the file access mode
 of fopen() from the C library.

 NOTE: mode = "wb" for write + binary format and "rb" for read + binary format.
 Extension used is .bin otherwise .dat is used assuming text mode ("w" or "r").

******************************************************************************/

FILE *fgh_basis_file(const char dir[], const char arrang, const size_t n,
                     const size_t J, const char mode[], const bool verbose)
{
	char filename[MAX_LINE_LENGTH], ext[4];

	if (strlen(mode) > 1 && mode[1] == 'b')
		sprintf(ext, "%s", "bin");
	else
		sprintf(ext, "%s", "dat");

	sprintf(filename, FGH_BASIS_FORMAT, dir, arrang, n, J, ext);

	FILE *stream = fopen(filename, mode);

	if (verbose)
	{
		if (stream != NULL && mode[0] == 'r') printf("# Reading %s\n", filename);
		if (stream != NULL && mode[0] == 'w') printf("# Writing %s\n", filename);
	}

	return stream;
}

/******************************************************************************

 Function fgh_basis_read(): saves in the disk the FGH basis function for the
 n-th channel and a given arrangement with total angular momentum J.

******************************************************************************/

void fgh_basis_write(const fgh_basis *b, FILE *output)
{
	ASSERT(b != NULL)

	file_write(&b->v, sizeof(size_t), 1, output);
	file_write(&b->j, sizeof(size_t), 1, output);
	file_write(&b->l, sizeof(size_t), 1, output);
	file_write(&b->n, sizeof(size_t), 1, output);

	file_write(&b->r_min, sizeof(double), 1, output);
	file_write(&b->r_max, sizeof(double), 1, output);
	file_write(&b->r_step, sizeof(double), 1, output);

	file_write(&b->eigenval, sizeof(double), 1, output);

	file_write(&b->grid_size, sizeof(size_t), 1, output);

	file_write(b->eigenvec, sizeof(double), b->grid_size, output);
}

/******************************************************************************

 Function fgh_basis_read(): loads from the disk the FGH basis function for the
 n-th channel and a given arrangement with total angular momentum J.

******************************************************************************/

void fgh_basis_read(fgh_basis *b, FILE *input)
{
	ASSERT(b != NULL)

	file_read(&b->v, sizeof(size_t), 1, input, 0);
	file_read(&b->j, sizeof(size_t), 1, input, 0);
	file_read(&b->l, sizeof(size_t), 1, input, 0);
	file_read(&b->n, sizeof(size_t), 1, input, 0);

	file_read(&b->r_min, sizeof(double), 1, input, 0);
	file_read(&b->r_max, sizeof(double), 1, input, 0);
	file_read(&b->r_step, sizeof(double), 1, input, 0);

	file_read(&b->eigenval, sizeof(double), 1, input, 0);

	file_read(&b->grid_size, sizeof(size_t), 1, input, 0);

	b->eigenvec = allocate(b->grid_size, sizeof(double), false);

	file_read(b->eigenvec, sizeof(double), b->grid_size, input, 0);
}
