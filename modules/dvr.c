#include <omp.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>

#include "globals.h"
#include "matrix.h"
#include "dvr.h"

struct dvr
{
	int max_ch;
	matrix *eigenvec;
	double *eigenval, r_min, r_max;
};

/******************************************************************************

 Function fgh_dvr(): return all eigenvalues/eigenvectors of a discrete variable
 representation (DVR) of the Fourier grid Hamiltonian (FGH), H = T + V, subject
 to a given single channel potential energy and reduced mass. As shown in Eq.
 (6a) and (6b) of Ref. [3].

 NOTE: The m-th eigenvector is the m-th column in a grid_size-by-grid_size row-
 major matrix: eigenvec[n*grid_size + m], where n = m = [0, grid_size).

******************************************************************************/

dvr *dvr_fgh_alloc(const int grid_size,
                   const double grid_step,
                   const double pot_energy[],
                   const double mass,
                   const double r_min,
                   const double r_max)
{
	ASSERT(pot_energy != NULL)

	dvr *result = calloc(1, sizeof(dvr));
	ASSERT(result != NULL)

	result->eigenvec = matrix_alloc(1, grid_size, grid_size, false);
	ASSERT(result->eigenvec != NULL)

	register const double box_length
		= as_double(grid_size - 1)*grid_step;

	register const double factor
		= (M_PI*M_PI)/(mass*box_length*box_length);

	register const double nn_term
		= factor*as_double(grid_size*grid_size + 2)/6.0;

	for (int n = 0; n < grid_size; ++n)
	{
		matrix_set(result->eigenvec, n, n, nn_term + pot_energy[n]);

		for (int m = (n + 1); m < grid_size; ++m)
		{
			register const double nm = as_double((n + 1) - (m + 1));
			register const double nm_term = sin(nm*M_PI/as_double(grid_size));

			matrix_set(result->eigenvec, n, m, pow(-1.0, nm)*factor/pow(nm_term, 2));
			matrix_set(result->eigenvec, m, n, pow(-1.0, nm)*factor/pow(nm_term, 2));
		}
	}

	result->eigenval = matrix_symm_eigen(result->eigenvec, 'v');
	result->r_min = r_min;
	result->r_max = r_max;
	result->max_ch = 1;

	return result;
}

/******************************************************************************

 Function multich_fgh_dvr(): the same as fgh_dvr(), except that the potential
 energy is the one of a problem with max_ch channels within grid_size points:
 pot_energy[p].data[n*max_ch + m] where, p = [0, grid_size) and n = m = [0,
 max_ch). For details see Eq. (19a), Eq. (19b) and Eq. (19c) of Ref. [12].

 NOTE: in the multichannel problem there are (grid_size*max_ch) eigenvalues
 and the eigenvectors are stored in a (grid_size*max_ch)-by-(grid_size*max_ch)
 row-major matrix. Thus, each multichannel eigenvector will have one component
 in every channel of the problem following the order in which the channels (or
 blocks) are arranged in the original potential matrix representation.

******************************************************************************/

dvr *dvr_multich_fgh_alloc(const int max_ch,
                           const int grid_size,
                           const double grid_step,
                           const matrix pot_energy[],
                           const double mass)
{
	ASSERT(pot_energy != NULL)

	dvr *result = calloc(1, sizeof(dvr));
	ASSERT(result != NULL)

	result->eigenvec = matrix_alloc(1, grid_size*max_ch, grid_size*max_ch, false);
	ASSERT(result->eigenvec != NULL)

	register const double box_length
		= as_double(grid_size - 1)*grid_step;

	register const double factor
		= (M_PI*M_PI)/(mass*box_length*box_length);

	register const double nn_term
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
					register double nmpq = 0.0;

					if (n == m && p == q)
					{
						nmpq = nn_term + matrix_get(&pot_energy[n], p, p);
					}

					if (n != m && p == q)
					{
						register const double nm
							= as_double((n + 1) - (m + 1));

						register const double nm_term
							= sin(nm*M_PI/as_double(grid_size));

						nmpq = pow(-1.0, nm)*factor/pow(nm_term, 2);
					}

					if (n == m && p != q)
					{
						nmpq = matrix_get(&pot_energy[n], p, q);
					}

					matrix_set(result->eigenvec, row, col, nmpq);
					++col;
				}
			}

			++row;
		}
	}

	result->eigenval = matrix_symm_eigen(result->eigenvec, 'v');
	result->max_ch = max_ch;
	result->r_min = r_min;
	result->r_max = r_max;

	return result;
}

/******************************************************************************

 Function fgh_dvr_wavef(): interpolate one eigenvector computed by fgh_dvr()
 using Eq. (4.1) of Ref. [5] and return its amplitude value at a given new r.

******************************************************************************/

double dvr_fgh_eigenvec(const dvr *g, const int v, const double r_new)
{
	register double result = 0.0;

	for (int n = 0; n < g->grid_size; ++n)
	{
		register const double wavef
			= matrix_get(g->eigenvec, n, v);

		register const double r_old
			= g->r_min + as_double(n)*g->r_step;

		register const double param
			= M_PI*(r_new - r_old)/g->r_step;

		result += (fabs(param) > 1.0E-7? wavef*sinc(param) : wavef);
	}

	return result;
}

/******************************************************************************

 Function fgh_dvr_product(): interpolate two eigenvectors computed by fgh_dvr()
 using Eq. (4.1) of Ref. [5] and return the product of their amplitude value at
 a given new r.

******************************************************************************/

double fgh_dvr_product(const int grid_size,
                       const double eigenvec_a[],
                       const double eigenvec_b[],
                       const double r_min,
                       const double r_max,
                       const double r_new)
{
	ASSERT(grid_size > 0)
	ASSERT(eigenvec_a != NULL)
	ASSERT(eigenvec_b != NULL)

	register double a = 0.0, b = 0.0;

	register const double r_inc
		= (r_max - r_min)/as_double(grid_size);

	for (int n = 0; n < grid_size; ++n)
	{
		register const double r_old
			= r_min + as_double(n)*r_inc;

		register const double param
			= M_PI*(r_new - r_old)/r_inc;

		if (fabs(param) > 1.0e-7)
		{
			a += eigenvec_a[n]*sinc(param);
			b += eigenvec_b[n]*sinc(param);
		}
		else
		{
			a += eigenvec_a[n];
			b += eigenvec_b[n];
		}
	}

	return a*b;
}
