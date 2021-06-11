/******************************************************************************

 About
 -----

 This module is a collection of routines applying the algorithms published by
 B. R. Johnson on regarding computation of bound and scattering wavefunctions
 (either single or multichannel).


 References
 ----------

 [1] B. R. Johnson J. Chem. Phys. 69, 4678 (1978)
     doi: https://doi.org/10.1063/1.436421

 [2] B. R. Johnson J. Comp. Phys. 13, 445-449 (1973)
     doi: https://doi.org/10.1016/0021-9991(73)90049-1

 [3] B. R. Johnson J. Chem. Phys. 67, 4086 (1977)
     doi: https://doi.org/10.1063/1.435384

 [4] B. R. Johnson. J. Comp. Phys. 13, 445 (1973)
     doi: https://doi.org/10.1016/0021-9991(73)90049-1

******************************************************************************/

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>

#include "math.h"
#include "matrix.h"
#include "johnson.h"

/******************************************************************************

 Function centr_term(): returns the actual centrifugal potential term at x for
 a given angular momentum l, provided the mass under consideration.

******************************************************************************/

inline static double centr_term(const int l,
                                const double mass, const double x)
{
	return as_double(l*(l + 1))/(2.0*mass*x*x);
}

/******************************************************************************

 Function johnson_riccati_bessel(): returns the Riccati-Bessel function J(l, x)
 if type = 'j', or N(l, x) if type = 'n', as defined in Eq. (14) and (15) of
 Ref. [2]. These are made upon spherical Bessel functions, j(l, x) and
 y(l, x).

******************************************************************************/

double johnson_riccati_bessel(const char type, const int l,
                              const double wavenum, const double x)
{
	switch (type)
	{
		case 'j': return  (wavenum*x)*gsl_sf_bessel_jl(l, wavenum*x)/sqrt(wavenum);
		case 'n': return -(wavenum*x)*gsl_sf_bessel_yl(l, wavenum*x)/sqrt(wavenum);

		default: return 0.0;
	}
}

/******************************************************************************

 Function johnson_modif_spher_bessel(): returns the modified spherical Bessel
 functions J(l, x) if type = 'j', or N(l, x) if type = 'n', as defined in Eq.
 (16) and (17) of Ref. [2]. These are made upon regular and irregular modified
 spherical Bessel functions, i(l, x) and k(l, x).

******************************************************************************/

double johnson_modif_spher_bessel(const char type, const int l,
                                  const double wavenum, const double x)
{
	switch (type)
	{
		case 'j': return sqrt(wavenum*x)*gsl_sf_bessel_il_scaled(l, wavenum*x);
		case 'n': return sqrt(wavenum*x)*gsl_sf_bessel_kl_scaled(l, wavenum*x);

		default: return 0.0;
	}
}

/******************************************************************************

 Function johnson_jcp77_numerov(): use the method of B. R. Johnson, Ref. [3],
 to construct the single channel wavefunction at a given trial energy. The
 return pointer contains the unormalized amplitude at grid_size points.

 NOTE: on exit, a negative number of nodes implies that no matching point has
 been found (see Ref. [3]).

******************************************************************************/

double *johnson_jcp77_numerov(const size_t grid_size,
                              const double grid_step,
                              const double pot_energy[],
                              const double trial_energy,
                              const double mass, double *error, size_t *nodes)
{
	ASSERT(pot_energy != NULL)

	double *T = allocate(grid_size, sizeof(double), false);
	double *F = allocate(grid_size, sizeof(double), false);

	double *inward_R  = allocate(grid_size, sizeof(double), false);
	double *outward_R = allocate(grid_size, sizeof(double), false);

/*
 *	Solve Eq. (32), (35) and (46) for each T, U and R (inward) coefficient:
 */

	inward_R[grid_size - 1] = 1.0E20;

	T[grid_size - 1] =
		-pow(grid_step, 2)*2.0*mass*(trial_energy - pot_energy[grid_size - 1])/12.0;

	size_t M = 0;
	for (size_t n = (grid_size - 2); n < grid_size; --n)
	{
		T[n] = -pow(grid_step, 2)*2.0*mass*(trial_energy - pot_energy[n])/12.0;

		inward_R[n] = (2.0 + 10.0*T[n])/(1.0 - T[n]) - 1.0/inward_R[n + 1];

/*
 *		Check if the matching point, M, is found:
 */

		if (inward_R[n] <= -1.0)
		{
			M = n;
			break;
		}
	}

	if (M == 0)
	{
		free(T);
		free(F);
		free(inward_R);
		free(outward_R);

		*error =  1.0E20;
		*nodes = -1;
		return NULL;
	}

/*
 *	Solve Eq. (32), (35) and (37) for each T, U and R (outward) coefficient:
 */

	outward_R[0] = 1.0E20;

	T[0] = -pow(grid_step, 2)*2.0*mass*(trial_energy - pot_energy[0])/12.0;

	*nodes = 0;
	for (size_t n = 1; n < M; ++n)
	{
		T[n] = -pow(grid_step, 2)*2.0*mass*(trial_energy - pot_energy[n])/12.0;

		outward_R[n] = (2.0 + 10.0*T[n])/(1.0 - T[n]) - 1.0/outward_R[n - 1];

/*
 *		Count the number of nodes in the wavefunction:
 */

		if (outward_R[n] < 0.0) *nodes += 1;
	}

/*
 *	Solve Eq. (36) and (45) for each F coefficient:
 */

	F[M] = 1.0;

	for (size_t n = (M - 1); n < M; --n)
		F[n] = F[n + 1]/outward_R[n];

	for (size_t n = (M + 1); n < grid_size; ++n)
		F[n] = F[n - 1]/inward_R[n];

/*
 *	Solve Eq. (33) for the wavefunction (unnormalized):
 */

	double *wavef = allocate(grid_size, sizeof(double), false);

	for (size_t n = 0; n < grid_size; ++n)
		wavef[n] = F[n]/(1.0 - T[n]);

/*
 *	Solve Eq. (48) for the difference, D (named as error here), between the
 *	left and right integrations:
 */

	*error = (0.5 - T[M + 1])*(1.0/inward_R[M + 1] - outward_R[M])/(1.0 - T[M + 1])
	       - (0.5 - T[M - 1])*(inward_R[M] - 1.0/outward_R[M - 1])/(1.0 - T[M - 1]);

	*error *= (1.0 - T[M]);

	free(T);
	free(F);
	free(inward_R);
	free(outward_R);

	return wavef;
}

/******************************************************************************

 Function johnson_jcp78_numerov(): use the method of B. R. Johnson, Ref. [1],
 to propagate the ratio matrix of a multichannel wavefunction from the radial
 grid point (n - 1) to n at a given total energy.

 NOTE: the potential matrix is used as workspace and its content is destroyed.

******************************************************************************/

void johnson_jcp78_numerov(const double grid_step, matrix *pot_energy,
                           const double tot_energy, const double mass, matrix *ratio)
{
	ASSERT(ratio != NULL)
	ASSERT(pot_energy != NULL)

	if (!matrix_is_null(ratio)) matrix_inverse(ratio);

/*
 *	NOTE: From Eq. (2) and (17) of Ref. [1] the following numerical
 *	factor is defined in atomic units:
 */

	const double factor = -grid_step*grid_step*2.0*mass/12.0;

/*
 *	Resolve Eq. (23) of Ref. [1] with Eq. (2) and (17) plugged in:
 */

	matrix *w = pot_energy;

	for (size_t n = 0; n < matrix_rows(pot_energy); ++n)
	{
		double t = factor*(tot_energy - matrix_get(pot_energy, n, n));
		matrix_set(w, n, n, 1.0 - t);

		for (size_t m = (n + 1); m < matrix_cols(pot_energy); ++m)
		{
			t = factor*(0.0 - matrix_get(pot_energy, n, m));
			matrix_set(w, n, m, 0.0 - t);

			t = factor*(0.0 - matrix_get(pot_energy, m, n));
			matrix_set(w, m, n, 0.0 - t);
		}
	}

/*
 *	Solve Eq. (22) and (24) of Ref. [1]:
 */

	matrix_inverse(w);

	for (size_t n = 0; n < matrix_rows(pot_energy); ++n)
	{
		double u = 12.0*matrix_get(w, n, n) - 10.0;
		matrix_set(ratio, n, n, u - matrix_get(ratio, n, n));

		for (size_t m = (n + 1); m < matrix_cols(pot_energy); ++m)
		{
			u = 12.0*matrix_get(w, n, m);
			matrix_set(ratio, n, m, u - matrix_get(ratio, n, m));

			u = 12.0*matrix_get(w, m, n);
			matrix_set(ratio, m, n, u - matrix_get(ratio, m, n));
		}
	}

	w = NULL;
}

/******************************************************************************

 Function johnson_jcp73_logd(): use the algorithm of B. R. Johnson, Ref. [4],
 in order to propagate the multichannel log derivative matrix, Y, from the
 radial grid point (n - 1) to n driven by the interaction potential V.

 NOTE: matrix V is defined as Eq. (2) of Ref. [4].

******************************************************************************/

void johnson_jcp73_logd(const int n, const int grid_size,
                        const double grid_step, const matrix *pot_energy, matrix *y)
{
	ASSERT(n > 0)
	ASSERT(y != NULL)
	ASSERT(grid_size >= n)
	ASSERT(pot_energy != NULL)
	ASSERT(GSL_IS_ODD(grid_size) == 1)

	const int max_ch = matrix_rows(pot_energy);

	/* NOTE: a = (I + hY)^-1 from the left term of Eq. (6). */
	matrix *a = matrix_alloc(max_ch, max_ch, false);

	matrix *eq6_right_term = matrix_alloc(max_ch, max_ch, false);

/*
 *	Solve Eq. (7) of Ref. [4]:
 */

	if (GSL_IS_EVEN(n) == 1)
	{
		const double weight = (n == 0? 1.0 : 2.0);
		const double factor = grid_step*weight/3.0;

		for (int p = 0; p < max_ch; ++p)
		{
			matrix_set(a, p, p, 1.0 + matrix_get(y, p, p)*grid_step);

			matrix_set(eq6_right_term, p, p, matrix_get(pot_energy, p, p)*factor);

			for (int q = (p + 1); q < max_ch; ++q)
			{
				matrix_set(a, p, q, matrix_get(y, p, q)*grid_step);
				matrix_set(a, q, p, matrix_get(y, q, p)*grid_step);

				matrix_set(eq6_right_term, p, q, matrix_get(pot_energy, p, q)*factor);
				matrix_set(eq6_right_term, q, p, matrix_get(pot_energy, q, p)*factor);
			}
		}
	}
	else
	{
		const double weight = (n == grid_size? 1.0 : 4.0);
		const double factor = grid_step*grid_step/6.0;

		/* NOTE: b = [I + (h^2/6)V]^-1 from the right term of Eq. (6). */
		matrix *b = matrix_alloc(max_ch, max_ch, false);

		for (int p = 0; p < max_ch; ++p)
		{
			matrix_set(a, p, p, 1.0 + matrix_get(y, p, p)*grid_step);

			matrix_set(b, p, p, 1.0 + matrix_get(pot_energy, p, p)*factor);

			for (int q = (p + 1); q < max_ch; ++q)
			{
				matrix_set(a, p, q, matrix_get(y, p, q)*grid_step);
				matrix_set(a, q, p, matrix_get(y, q, p)*grid_step);

				matrix_set(b, p, q, matrix_get(pot_energy, p, q)*factor);
				matrix_set(b, q, p, matrix_get(pot_energy, q, p)*factor);
			}
		}

		matrix_inverse(b);
		matrix_multiply(grid_step*weight/3.0, b, pot_energy, 0.0, eq6_right_term);

		matrix_free(b);
	}

	matrix *eq6_left_term = matrix_alloc(max_ch, max_ch, false);

	matrix_inverse(a);
	matrix_multiply(1.0, a, y, 0.0, eq6_left_term);

	matrix_free(a);

/*
 *	Solve Eq. (6) of Ref. [4]:
 */

	matrix_sub(1.0, eq6_left_term, 1.0, eq6_right_term, y);

	matrix_free(eq6_left_term);
	matrix_free(eq6_right_term);
}

/******************************************************************************

 Function johnson_kmatrix(): build the augmented reactant matrix K as described
 by Eq. (A19) of Ref. [1], from the ratio matrix of the multichannel scattering
 wavefunction. Where, l and level are the asymptotic angular momentum and eigen
 value of each channel, respectively, as R -> inf.

 NOTE: both l and level are vectors with the same size of the rows (or columns)
 of the ratio matrix, i.e. total number of channels.

******************************************************************************/

matrix *johnson_kmatrix(const int l[],
                        const double grid_step,
                        const double tot_energy,
                        const double mass,
                        const double level[],
                        const matrix *ratio,
                        const double R)
{
	ASSERT(l != NULL)
	ASSERT(ratio != NULL)
	ASSERT(level != NULL)

	const int max_ch = matrix_rows(ratio);

/*
 *	NOTE: from Eq. (2) and (17) of Ref. [1] the following numerical factor is
 *	defined in atomic units:
 */

	const double factor
		= -grid_step*grid_step*2.0*mass/12.0;

/*
 *	Step 1: build the j(R) and n(R) diagonal matrices defined by Eq. (A16) and
 *	(A17) of Ref. [1]; where, Eq. (23) is also used:
 */

	matrix *n = matrix_alloc(max_ch, max_ch, true);
	matrix *j = matrix_alloc(max_ch, max_ch, true);

	for (int i = 0; i < max_ch; ++i)
	{
		if (level[i] < tot_energy)
		{
			const double wavenum
				= sqrt(2.0*mass*(tot_energy - level[i]));

			const double w
				= 1.0 - factor*(tot_energy - centr_term(l[i], mass, R));

			matrix_set(n, i, i, w*johnson_riccati_bessel('n', l[i], wavenum, R));
			matrix_set(j, i, i, w*johnson_riccati_bessel('j', l[i], wavenum, R));
		}
		else
		{
			matrix_set(n, i, i, 1.0);
			matrix_set(j, i, i, 1.0);
		}
	}

/*
 *	Step 2: build the product rn and rj (r := ratio) at the grid point R, as
 *	shown in Eq. (A19) of Ref. [1]:
 */

	matrix *rn = matrix_alloc(max_ch, max_ch, false);
	matrix *rj = matrix_alloc(max_ch, max_ch, false);

	matrix_multiply(1.0, ratio, n, 0.0, rn);
	matrix_multiply(1.0, ratio, j, 0.0, rj);

	matrix_free(j);
	matrix_free(n);

/*
 *	Step 3: subtract j(R + grid_step) and n(R + grid_step) diagonal matrices
 *	from the rn(R) and rj(R) products, following Eq. (A19) of Ref. [1]:
 */

	for (int i = 0; i < max_ch; ++i)
	{
		const double wavenum
			= sqrt(2.0*mass*(tot_energy - level[i]));

		const double w_prime
			= 1.0 - factor*(tot_energy - centr_term(l[i], mass, R + grid_step));

		const double n_prime
			= w_prime*johnson_riccati_bessel('n', l[i], wavenum, R + grid_step);

		const double j_prime
			= w_prime*johnson_riccati_bessel('j', l[i], wavenum, R + grid_step);

		if (level[i] < tot_energy)
		{
			matrix_decr(rn, i, i, n_prime);
			matrix_decr(rj, i, i, j_prime);
		}
		else
		{
			const double w
				= 1.0 - factor*(tot_energy - centr_term(l[i], mass, R));

			matrix_decr(rn, i, i, n_prime/(w*johnson_modif_spher_bessel('n', l[i], wavenum, R)));
			matrix_decr(rj, i, i, j_prime/(w*johnson_modif_spher_bessel('j', l[i], wavenum, R)));
		}
	}

/*
 *	Step 4: resolve Eq. (A19) of Ref. [1] for the K matrix:
 */

	matrix *k = matrix_alloc(max_ch, max_ch, false);

	matrix_inverse(rn);
	matrix_multiply(-1.0, rn, rj, 0.0, k);

	matrix_free(rn);
	matrix_free(rj);
	return k;
}

/******************************************************************************

 Function johnson_smatrix(): return both real and imaginary parts of a
 scattering matrix, S, defined as

 S := (I + iK)inv(I - iK)
    = (I - KK)inv(I + KK) + i(2K)inv(I + KK)

 Where, I is the unit matrix and K is the open-open block of a reactant matrix.

 For details see Eq. (20) of Ref. [2].

******************************************************************************/

smatrix *johnson_smatrix(const matrix *k)
{
	const int max_ch = matrix_rows(k);

	matrix *a = matrix_alloc(max_ch, max_ch, false);
	matrix *b = matrix_alloc(max_ch, max_ch, false);
	matrix *c = matrix_alloc(max_ch, max_ch, false);

/*
 *	Resolve A = KK, B = (I + A) and C = (I - A):
 */

	matrix_multiply(1.0, k, k, 0.0, a);

	for (int n = 0; n < max_ch; ++n)
	{
		matrix_set_diag(b, n, 1.0 + matrix_get(a, n, n));
		matrix_set_diag(c, n, 1.0 - matrix_get(a, n, n));

		for (int m = (n + 1); m < max_ch; ++m)
		{
			matrix_set_symm(b, n, m, 0.0 + matrix_get(a, n, m));
			matrix_set_symm(c, n, m, 0.0 - matrix_get(a, n, m));
		}
	}

/*
 *	Build the S matrix, Re(S) = C*inv(B) and Im(S) = 2*K*inv(B):
 */

	matrix_inverse(b);

	smatrix *s = calloc(1, sizeof(smatrix));

	s->re_part = matrix_alloc(max_ch, max_ch, false);
	s->im_part = matrix_alloc(max_ch, max_ch, false);

	matrix_multiply(1.0, c, b, 0.0, s->re_part);
	matrix_multiply(2.0, k, b, 0.0, s->im_part);

	matrix_free(a);
	matrix_free(b);
	matrix_free(c);
	return s;
}
