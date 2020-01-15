/******************************************************************************

 About
 -----

 This module is a collection of routines applying algorithms based on the
 discrete variable representation (DVR) for atomic and molecular systems.


 References
 ----------

 [1] O. Dulieu et al. J. Chem. Phys. 103 (1) (1995)
     doi: http://dx.doi.org/10.1063/1.469622

 [2] M. Monnerville et al. J. Chem. Phys. 101, 7580 (1994)
     doi: http://dx.doi.org/10.1063/1.468252

 [3] V. Kokoouline et al. J. Chem. Phys. 110, 20 (1999)
     doi:

******************************************************************************/

#include "matrix.h"
#include "dvr.h"

/******************************************************************************

 Function dvr_fgh(): return the discrete variable representation (DVR) of the
 Fourier grid Hamiltonian (FGH) subjected to a given single channel potential
 energy and reduced mass. As shown in Eq. (6a) and (6b) of Ref. [1].

 NOTE: The m-th eigenvector is the m-th column in a grid_size-by-grid_size row-
 major matrix: eigenvec[n*grid_size + m], where n = m = [0, grid_size).

******************************************************************************/

matrix *dvr_fgh(const int grid_size,
                const double grid_step,
                const double pot_energy[],
                const double mass)
{
	ASSERT(pot_energy != NULL)

	matrix *result
		= matrix_alloc(grid_size, grid_size, false);

	const double box_length
		= as_double(grid_size - 1)*grid_step;

	const double factor
		= (M_PI*M_PI)/(mass*box_length*box_length);

	const double nn_term
		= factor*as_double(grid_size*grid_size + 2)/6.0;

	for (int n = 0; n < grid_size; ++n)
	{
		matrix_diag_set(result, n, nn_term + pot_energy[n]);

		for (int m = (n + 1); m < grid_size; ++m)
		{
			const double nm
				= as_double((n + 1) - (m + 1));

			const double nm_term
				= sin(nm*M_PI/as_double(grid_size));

			matrix_symm_set(result, n, m, pow(-1.0, nm)*factor/pow(nm_term, 2));
		}
	}

	return result;
}

/******************************************************************************

 Function dvr_multich_fgh(): the same as dvr_fgh(), except that the potential
 energy is the one of a problem with max_ch channels within grid_size points.
 For details see Eq. (19a), Eq. (19b) and Eq. (19c) of Ref. [2].

******************************************************************************/

matrix *dvr_multich_fgh(const int max_ch,
                        const int grid_size,
                        const double grid_step,
                        const tensor pot_energy[],
                        const double mass)
{
	ASSERT(max_ch > 0)
	ASSERT(pot_energy != NULL)

	matrix *result
		= matrix_alloc(grid_size*max_ch, grid_size*max_ch, false);

	const double box_length
		= as_double(grid_size - 1)*grid_step;

	const double factor
		= (M_PI*M_PI)/(mass*box_length*box_length);

	const double nn_term
		= factor*as_double(grid_size*grid_size + 2)/6.0;

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
						const double nm
							= as_double((n + 1) - (m + 1));

						const double nm_term
							= sin(nm*M_PI/as_double(grid_size));

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

 Function dvr_multich_fgh_comp(): return the n-th component of a multichannel
 eigenvector v from the Hamiltonian built by dvr_multich_fgh(). On entry, the
 matrix is expected to be properly diagonalized and its columns the respective
 eigenvectors.

******************************************************************************/

matrix *dvr_multich_fgh_comp(matrix *fgh,
                             const int max_ch, const int n, const int v)
{
	ASSERT(n < max_ch)
	const int grid_size = matrix_row(fgh)/max_ch;

	return matrix_get_block(fgh, n*grid_size, n*grid_size + (grid_size - 1), v, v);
}

/******************************************************************************

 Function dvr_fgh_wavef(): interpolate the eigenvector a of the Hamiltonian
 built by dvr_fgh(), using Eq. (4.1) of Ref. [3], and return its amplitude
 value at a given new r.

******************************************************************************/

double dvr_fgh_wavef(const matrix *fgh, const int a,
                     const double r_min, const double r_max, const double r_new)
{
	const int n_max = matrix_row(fgh);

	ASSERT(a > -1)
	ASSERT(a < n_max)

	const double r_step = (r_max - r_min)/as_double(n_max);

	double result = 0.0;

	for (int n = 0; n < n_max; ++n)
	{
		const double wavef_a = matrix_get(fgh, n, a);

		const double r_old = r_min + as_double(n)*r_step;
		const double param = M_PI*(r_new - r_old)/r_step;

		result += (fabs(param) > 1.0E-7? wavef_a*sinc(param) : wavef_a);
	}

	return result;
}

/******************************************************************************

 Function dvr_fgh_product(): the same of dvr_fgh_wavef() but for the product of
 two eigenvectors, a and b.

******************************************************************************/

double dvr_fgh_product(const matrix *fgh, const int a, const int b,
                       const double r_min, const double r_max, const double r_new)
{
	const int n_max = matrix_row(fgh);

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
}

/******************************************************************************

 Function dvr_fgh_norm(): normalize to unity the eigenvector a of a Hamiltonian
 built by dvr_fgh(). On entry, the matrix is expected to be properly
 diagonalized and its columns the respective eigenvectors.

******************************************************************************/

void dvr_fgh_norm(matrix *fgh,
                  const int a, const double grid_step, const bool use_omp)
{
	const int n_max
		= (matrix_row(fgh)%2 == 0? matrix_row(fgh) : matrix_row(fgh) - 1);

	double sum = matrix_get_pow(fgh, 0, a, 2.0)
	           + matrix_get_pow(fgh, n_max - 1, a, 2.0);

	#pragma omp parallel for default(none) shared(fgh) reduction(+:sum) if(use_omp)
	for (int n = 1; n < (n_max - 2); n += 2)
	{
		sum += 4.0*matrix_get_pow(fgh, n, a, 2.0)
		     + 2.0*matrix_get_pow(fgh, n + 1, a, 2.0);
	}

	sum = grid_step*sum/3.0;
	matrix_col_scale(fgh, a, 1.0/sqrt(sum), use_omp);
}

/******************************************************************************

 Function dvr_fgh_norm(): the same of dvr_fgh_norm() but for eigenvectos from
 those matrices built by dvr_multich_fgh().

******************************************************************************/

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
}
