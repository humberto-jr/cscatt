/******************************************************************************

 About
 -----

 This module is a collection of special math functions for problems of AMO and
 quantum physics as well as routines to perform numerical integration, either
 using quadratures or statistical methods.

******************************************************************************/

#include "math.h"
#include "gsl_lib.h"

//#include "arpack/arpack.h"

static int workspace_size = 5000;
static double abs_error = 1.0E-6;

/******************************************************************************

 Function math_legendre_poly(): returns a Legendre polynomial at x in [-1, 1].
 Where, l is positive.

******************************************************************************/

double math_legendre_poly(const int l, const double x)
{
	ASSERT(l >= 0)
	ASSERT(fabs(x) <= 1.0)

	return gsl_sf_legendre_Pl(l, x);
}

/******************************************************************************

 Function math_sphe_harmonics(): returns the spherical harmonics y for angular
 momentum l and projection m at (theta, phi).

******************************************************************************/

double math_sphe_harmonics(const int l, const int m,
                           const double theta, const double phi)
{
	ASSERT(l >= abs(m))

	const double x
		= cos(theta*M_PI/180.0);

	const double m_phase
		= (m > 0? pow(-1.0, m) : 1.0);

	const double theta_wavef
		= gsl_sf_legendre_sphPlm(l, abs(m), x);

	const double phi_wavef
		= exp(as_double(m)*phi*M_PI/180.0)/sqrt(2.0*M_PI);

	/* NOTE: see Eq. 1.43 (pag. 8) of Angular Momentum by Richard N. Zare. */
	return m_phase*theta_wavef*phi_wavef;
}

/******************************************************************************

 Function math_wigner_3j(): returns the Wigner 3j-symbol of coupling a and b to
 result in c, where d, e and f are the projections of a, b, c, respectively.

******************************************************************************/

double math_wigner_3j(const int a, const int b, const int c,
                      const int d, const int e, const int f)
{
	return gsl_sf_coupling_3j(2*a, 2*b, 2*c, 2*d, 2*e, 2*f);
}

/******************************************************************************

 Function math_wigner_6j(): similarly to math_wigner_3j(), returns the Wigner
 6j-symbol.

******************************************************************************/

double math_wigner_6j(const int a, const int b, const int c,
                      const int d, const int e, const int f)
{
	return gsl_sf_coupling_6j(2*a, 2*b, 2*c, 2*d, 2*e, 2*f);
}

/******************************************************************************

 Function math_wigner_9j(): similarly to math_wigner_6j(), returns the Wigner
 9j-symbol.

******************************************************************************/

double math_wigner_9j(const int a, const int b, const int c,
                      const int d, const int e, const int f,
                      const int g, const int h, const int i)

{
	return gsl_sf_coupling_9j(2*a, 2*b, 2*c, 2*d, 2*e, 2*f, 2*g, 2*h, 2*i);
}

/******************************************************************************

 Function math_clebsch_gordan(): returns the Clebsch–Gordan coefficients that
 arise in the angular momentum coupling of j1 and j2 to result in j3. Where, m
 is the respective projection.

******************************************************************************/

double math_clebsch_gordan(const int j1, const int j2, const int j3,
                           const int m1, const int m2, const int m3)
{
	return pow(-1.0, j1 - j2 + m3)*sqrt(2*j3 + 1)*math_wigner_3j(j1, j2, j3, m1, m2, -m3);
}

/******************************************************************************

 Function math_racah_coeff():

******************************************************************************/

double math_racah_coeff(const int J,
                        const int j1,
                        const int j2,
                        const int l1,
                        const int l2,
                        const int lambda)
{
	return pow(-1.0, j1 + l1 + j2 + l2)*math_wigner_6j(j1, l1, J, l2, j2, lambda);
}

/******************************************************************************

 Function math_percival_seaton(): return the so-called Percival & Seaton term
 often used in the definition of the atom-diatom collisional coupling matrix.
 Where,

 j1     = diatomic rotational angular momentum quantum number for channel 1
 j2     = diatomic rotational angular momentum quantum number for channel 2
 l1     = atom-diatom orbital angular momentum quantum number for channel 1
 l2     = atom-diatom orbital angular momentum quantum number for channel 2
 J      = total angular momentum of the problem
 lambda = an integer parameter (0, 1, 2, ...)

 and also depends parametrically on the electronic spin multiplicity (1 for
 singlet, 2 for doublet and 3 for triplet).

 See Ref. [1] for more, in particular, Eq. (A1).

******************************************************************************/

double math_percival_seaton(const int J,
                            const int j1,
                            const int j2,
                            const int l1,
                            const int l2,
                            const int lambda)
{
	double result = pow(-1.0, j1 + j2 - J);

	result *= math_wigner_3j(l1, l2, lambda, 0, 0, 0);
	result *= math_wigner_3j(j1, j2, lambda, 0, 0, 0);
	result *= math_wigner_6j(j1, j2, lambda, l2, l1, J);
	result *= sqrt(as_double(2*j1 + 1)*as_double(2*j2 + 1)*as_double(2*l1 + 1)*as_double(2*l2 + 1));

	return result;
}

/******************************************************************************

 Function math_wigner_d():

******************************************************************************/

double *math_wigner_d(const int k,
                      const int m,
                      const int j_max,
                      const double beta)
{
	ASSERT(m >= 0)
	ASSERT(k >= m)
	ASSERT(j_max > -1)
	ASSERT(j_max < 46340)

	double dkm[3], *result = allocate(j_max + 1, sizeof(double), true);

	/* NOTE: from degree to radian */
	const double x = beta*M_PI/180.0;

	/* Eq. (20) */
	const double seed_c = pow(cos(x/2.0), k + m);

	/* Eq. (20) */
	const double seed_s = pow(sin(x/2.0), k - m);

	/* Eq. (21) */
	const double seed_e = sqrt(factorial(2*k)/(factorial(k + m)*factorial(k - m)));

	const double t = 2.0*pow(sin(x/2.0), 2);

	/* Eq. (18) */
	dkm[0] = seed_c*seed_s*seed_e;

	result[k/2] = dkm[0];

	if (j_max <= k) return result;

	int j = k + 2;
	dkm[1] = dkm[0]*sqrt(as_double(j - 1)/as_double((j + m)*(j - m)))*(as_double(j - m) - as_double(j)*t);

	result[j/2] = dkm[1];

	const double i = as_double(k*m);

	for (j = (k + 4); j <= j_max; j += 2)
	{
		const double h = as_double(j*(j - 2));

		const double g = as_double(j + k);

		const double f = as_double(j - k);

		const double e = as_double(j + m);

		const double d = as_double(j - m);

		const double c = 1.0/(as_double(j - 2)*sqrt(g*f*e*d));

		const double b = c*as_double(j)*sqrt((g - 2.0)*(f - 2.0)*(e - 2.0)*(d - 2.0));

		const double a = c*as_double(2*j - 2)*((h - i) - h*t);

		dkm[2] = a*dkm[1] - b*dkm[0];

		result[j/2] = dkm[2];

		dkm[0] = dkm[1];
		dkm[1] = dkm[2];
	}

	return result;
}

/******************************************************************************

 Function math_integral_yyy():

******************************************************************************/

double math_integral_yyy(const int j1, const int m1,
                         const int j2, const int m2,
                         const int j3, const int m3)
{
	const double a = as_double(2*j1 + 1);

	const double b = as_double(2*j2 + 1);

	const double c = as_double(2*j3 + 1);

	const double d = sqrt(a*b*c/(4.0*M_PI));

	const double e = math_wigner_3j(j1, j2, j3, 0, 0, 0);

	const double f = math_wigner_3j(j1, j2, j3, m1, m2, m3);

/*	const double a = as_double(2*j1 + 1);

	const double b = as_double(2*j2 + 1);

	const double c = as_double(2*j3 + 1);

	const double d = sqrt(a*b/(4.0*M_PI*c));

	const double e = math_clebsch_gordan(j1, j2, j3, 0, 0, 0);

	const double f = math_clebsch_gordan(j1, j2, j3, m1, m2, m3);*/

	return d*e*f;
}

/******************************************************************************

 Function math_set_error(): sets the absolute error for the QAG kind of methods
 (1e-6 by default). If a given abs. error cannot be reached, an error message
 is printed in the C stderr.

******************************************************************************/

void math_set_error(const double error)
{
	ASSERT(error > 0.0)
	abs_error = error;
}

/******************************************************************************

 Function math_set_workspace(): sets the size of auxiliary workspace used by
 integrators (5000 by default).

******************************************************************************/

void math_set_workspace(const int size)
{
	ASSERT(size > 0)
	workspace_size = size;
}

/******************************************************************************

 Function math_no_gsl_handler(): sets the GSL error handler off.

******************************************************************************/

void math_no_gsl_handler()
{
	gsl_set_error_handler_off();
}

/******************************************************************************

 Function math_simpson(): return the integral of f(x) from a to b, using the
 1/3-Simpson quadrature rule. Where, values of f are evaluated in a grid of n
 points and params points to a struct of parameters that f may depend on (NULL
 if not needed).

******************************************************************************/

double math_simpson(const int n,
                    const double a,
                    const double b,
                    void *params,
                    const bool use_omp,
                    double (*f)(double x, void *params))
{
	ASSERT(n%2 == 0)

	const double grid_step = (b - a)/as_double(n);
	double sum = f(a, params) + f(b, params);

	#pragma omp parallel for default(none) shared(f, params) reduction(+:sum) if(use_omp)
	for (int m = 1; m < (n - 2); m += 2)
	{
		sum += 4.0*f(a + as_double(m)*grid_step, params);
		sum += 2.0*f(a + as_double(m + 1)*grid_step, params);
	}

	sum += 4.0*f(a + as_double(n - 1)*grid_step, params);

	return grid_step*sum/3.0;
}

/******************************************************************************

 Function math_2nd_simpson(): the same as math_simpson() but using Simpson's
 second rule, i.e. 3/8-Simpson quadrature.

******************************************************************************/

double math_2nd_simpson(const int n,
                        const double a,
                        const double b,
                        void *params,
                        const bool use_omp,
                        double (*f)(double x, void *params))
{
	ASSERT(n%3 == 0)

	const double grid_step = (b - a)/as_double(n);
	double sum = f(a, params) + f(b, params);

	#pragma omp parallel for default(none) shared(f, params) reduction(+:sum) if(use_omp)
	for (int m = 1; m < (n - 2); m += 3)
	{
		sum += 3.0*f(a + as_double(m)*grid_step, params);
		sum += 3.0*f(a + as_double(m + 1)*grid_step, params);
		sum += 2.0*f(a + as_double(m + 2)*grid_step, params);
	}

	sum += 3.0*f(a + as_double(n - 1)*grid_step, params);
	sum += 3.0*f(a + as_double(n - 2)*grid_step, params);

	return grid_step*3.0*sum/8.0;
}

/******************************************************************************

 Function math_qag(): return the integral of f(x) from a to b, using a 61 point
 Gauss-Kronrod rule in a QAG framework. Where, params points to a struct of
 parameters that f may depend on (NULL if not needed).

******************************************************************************/

double math_qag(const double a,
                const double b,
                void *params,
                double (*f)(double x, void *params))
{
	gsl_function gsl_f =
	{
		.params = params,
		.function = f
	};

	gsl_integration_workspace *work
		= gsl_integration_workspace_alloc(workspace_size);

	ASSERT(work != NULL)
	double result = 0.0, error = 0.0;

	const int info = gsl_integration_qag(&gsl_f, a, b, abs_error, 0.0,
	                                     workspace_size, GSL_INTEG_GAUSS61, work, &result, &error);

	gsl_integration_workspace_free(work);

	if (info != 0)
	{
		PRINT_ERROR("%s; error = %f\n", gsl_strerror(info), error)
	}

	return result;
}

/******************************************************************************

 Function math_qags(): the same as integral_qag() but using a smaller order in
 the Gauss-Kronrod rule (21) and assuming that f may be singular.

******************************************************************************/

double math_qags(const double a,
                 const double b,
                 void *params,
                 double (*f)(double x, void *params))
{
	gsl_function gsl_f =
	{
		.params = params,
		.function = f
	};

	gsl_integration_workspace *work
		= gsl_integration_workspace_alloc(workspace_size);

	ASSERT(work != NULL)
	double result = 0.0, error = 0.0;

	const int info = gsl_integration_qags(&gsl_f, a, b, abs_error, 0.0,
	                                      workspace_size, work, &result, &error);

	gsl_integration_workspace_free(work);

	if (info != 0)
	{
		PRINT_ERROR("%s; error = %f\n", gsl_strerror(info), error)
	}

	return result;
}

/******************************************************************************

 Function math_plain_mcarlo(): return the n-dimensional integral of f(x) from
 a[0, 1, ..., n] to b[0, 1, ..., n], using a plain Monte Carlo algorithm for a
 given number of f calls.

******************************************************************************/

double math_plain_mcarlo(const int n,
                         const int calls,
                         const double a[],
                         const double b[],
                         void *params,
                         double (*f)(double x[], size_t n, void *params))
{
	gsl_monte_function gsl_f =
	{
		.params = params,
		.dim = n,
		.f = f
	};

	gsl_monte_plain_state *work = gsl_monte_plain_alloc(n);
	ASSERT(work != NULL)

	gsl_rng *r = gsl_rng_alloc(gsl_rng_ranlxd2);
	ASSERT(r != NULL)

	double result = 0.0, error = 0.0;

	const int info = gsl_monte_plain_integrate(&gsl_f, a, b, n,
	                                           calls, r, work, &result, &error);

	gsl_rng_free(r);
	gsl_monte_plain_free(work);

	if (info != 0)
	{
		PRINT_ERROR("%s; error = %f\n", gsl_strerror(info), error)
	}

	return result;
}

/******************************************************************************

 Function math_vegas_mcarlo(): the same as math_plain_mcarlo() but using the
 VEGAS algorithm.

******************************************************************************/

double math_vegas_mcarlo(const int n,
                         const int calls,
                         double a[],
                         double b[],
                         void *params,
                         double (*f)(double x[], size_t n, void *params))
{
	gsl_monte_function gsl_f =
	{
		.params = params,
		.dim = n,
		.f = f
	};

	gsl_monte_vegas_state *work = gsl_monte_vegas_alloc(n);
	ASSERT(work != NULL)

	gsl_rng *r = gsl_rng_alloc(gsl_rng_ranlxd2);
	ASSERT(r != NULL)

	double result = 0.0, error = 0.0;

	const int info = gsl_monte_vegas_integrate(&gsl_f, a, b, n,
	                                           calls, r, work, &result, &error);

	gsl_rng_free(r);
	gsl_monte_vegas_free(work);

	if (info != 0)
	{
		PRINT_ERROR("%s; error = %f\n", gsl_strerror(info), error)
	}

	return result;
}

/******************************************************************************

 Function math_miser_mcarlo(): the same as math_plain_mcarlo() but using the
 MISER algorithm.

******************************************************************************/

double math_miser_mcarlo(const int n,
                         const int calls,
                         const double a[],
                         const double b[],
                         void *params,
                         double (*f)(double x[], size_t n, void *params))
{
	gsl_monte_function gsl_f =
	{
		.params = params,
		.dim = n,
		.f = f
	};

	gsl_monte_miser_state *work = gsl_monte_miser_alloc(n);
	ASSERT(work != NULL)

	gsl_rng *r = gsl_rng_alloc(gsl_rng_ranlxd2);
	ASSERT(r != NULL)

	double result = 0.0, error = 0.0;

	const int info = gsl_monte_miser_integrate(&gsl_f, a, b, n,
	                                           calls, r, work, &result, &error);

	gsl_rng_free(r);
	gsl_monte_miser_free(work);

	if (info != 0)
	{
		PRINT_ERROR("%s; error = %f\n", gsl_strerror(info), error)
	}

	return result;
}

/******************************************************************************

 Function math_lanczos():

******************************************************************************/
/*
void math_lanczos(math_lanczos_setup *s)
{
	ASSERT(s != NULL)

	ASSERT(s->n > 0)
	ASSERT(s->n_max >= s->n)
	ASSERT(s->max_step > 0)
*/
/*
 *	NOTE: this driver uses the same naming convention used by ARPACK library
 *	except for the following: N = n_max, nev = n, bmat = "I" and which = "SM".
 *	Likewise, for dseupd() routine, rmat = true, howmny = "A" and sigma = 0.
 */
/*
	a_int ido = 0;
	a_int info = 1;
	a_int ncv = 2*s->n + 10;
	a_int lworkl = 3*ncv*ncv + 8*ncv;
	a_int ipntr[14] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	a_int iparam[11] = {1, 1, (a_int) s->max_step, 1, 0, 0, 1, 0, 0, 0, 0};

	double *v = allocate(s->n_max*ncv, sizeof(double), false);
	double *workd = allocate(3*s->n_max, sizeof(double), false);
	double *workl = allocate(lworkl, sizeof(double), false);

	if (s->start_vector == NULL)
	{
		info = 0;
		s->start_vector = allocate(s->n_max, sizeof(double), false);

		for (int p = 0; p < s->n_max; ++p)
		{
			s->start_vector[p] = (double) rand()/RAND_MAX;
		}
	}

	while (ido != (a_int) 99)
	{
		dsaupd_c(&ido, "I", s->n_max, "SM", s->n, abs_error, s->start_vector,
		         ncv, v, s->n_max, iparam, ipntr, workd, workl, lworkl, &info);

		double *x = workd + ((int) ipntr[0]) - 1;
		double *y = workd + ((int) ipntr[1]) - 1;

		s->product(s->n_max, x, y, s->params);
	}

	if (info != (a_int) 0)
	{
		PRINT_ERROR("dsaupd_c() failed with error code %d\n", (int) info)
		exit(EXIT_FAILURE);
	}

	a_int *select = allocate(ncv, sizeof(a_int), false);

	if (s->eigenval == NULL)
	{
		s->eigenval = allocate(s->n, sizeof(double), false);
	}

	if (s->eigenvec == NULL)
	{
		s->eigenvec = allocate(s->n_max*s->n, sizeof(double), false);
	}

	dseupd_c(1, "A", select, s->eigenval, s->eigenvec, s->n_max, 0.0, "I",
	         s->n_max, "SM", s->n, abs_error, s->start_vector, ncv, v,
	         s->n_max, iparam, ipntr, workd, workl, lworkl, &info);

	free(v);
	free(workd);
	free(workl);
	free(select);

	if (info != (a_int) 0)
	{
		PRINT_ERROR("dseupd_c() failed with error code %d\n", (int) info)
		exit(EXIT_FAILURE);
	}
}*/
