/******************************************************************************

 About
 -----

 This module is a collection of routines to perform numerical integration,
 either using quadratures or statistical methods.

******************************************************************************/

#include "integral.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_integration.h>

static int workspace_size = 5000;
static double abs_error = 1.0e-6;

/******************************************************************************

 Function integral_set_error(): sets the absolute error for the QAG kind of
 methods (1e-6 by default). If a given abs. error cannot be reached, an error
 message is printed in the C stderr.

******************************************************************************/

void integral_set_error(const double error)
{
	ASSERT(error > 0.0)
	abs_error = error;
}

/******************************************************************************

 Function integral_set_workspace(): sets the size of auxiliary workspace used
 during the integrations (5000 by default).

******************************************************************************/

void integral_set_workspace(const int size)
{
	ASSERT(size > 0)
	workspace_size = size;
}

/******************************************************************************

 Function integral_simpson(): return the integral of f(x) using the 1/3-Simpson
 quadrature rule. Where, values of f in a grid-mesh of grid_size x-points is
 given.

******************************************************************************/

double integral_simpson(const int grid_size,
                        const double grid_step,
                        const bool use_omp,
                        const double f[])
{
	ASSERT(f != NULL)
	ASSERT(grid_size > 6)

	const int n_max = (grid_size%2 == 0? grid_size : grid_size - 1);

	double sum = f[0] + f[n_max - 1];

	#pragma omp parallel for default(none) shared(f) reduction(+:sum) if(use_omp)
	for (int n = 1; n < (n_max - 2); n += 2)
	{
		sum += 4.0*f[n] + 2.0*f[n + 1];
	}

	return grid_step*sum/3.0;
}

/******************************************************************************

 Function integral_simpson_2nd(): the same as integral_simpson() but using
 Simpson's second rule, i.e. 3/8-Simpson quadrature.

******************************************************************************/

double integral_simpson_2nd(const int grid_size,
                            const double grid_step,
                            const bool use_omp,
                            const double f[])
{
	ASSERT(f != NULL)
	ASSERT(grid_size > 12)

	const int n_max = (grid_size%3 == 0? grid_size : grid_size - 2);

	double sum = f[0] + f[n_max - 1];

	#pragma omp parallel for default(none) shared(f) reduction(+:sum) if(use_omp)
	for (int n = 1; n < (n_max - 2); n += 3)
	{
		sum += 3.0*f[n] + 3.0*f[n + 1] + 2.0*f[n + 2];
	}

	return grid_step*3.0*sum/8.0;
}

/******************************************************************************

 Function integral_qag(): return the integral of f(x) from a to b, using a 61
 point Gauss-Kronrod rule in a QAG framework. Where, params points to a struct
 of parameters that f may depend on (NULL if not needed).

******************************************************************************/

double integral_qag(const double a,
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

 Function integral_qags(): the same as integral_qag() but using a smaller order
 in the Gauss-Kronrod rule (21) and assuming that f may be singular.

******************************************************************************/

double integral_qags(const double a,
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

 Function integral_plain_mcarlo(): return the n-dimensional integral of f(x)
 from a[0, 1, ..., n] to b[0, 1, ..., n], using a plain Monte Carlo algorithm
 for a given number of f calls.

******************************************************************************/

double integral_plain_mcarlo(const int n,
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

 Function integral_vegas_mcarlo(): the same as integral_plain_mcarlo() but
 using the VEGAS algorithm.

******************************************************************************/

double integral_vegas_mcarlo(const int n,
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

 Function integral_miser_mcarlo(): the same as integral_plain_mcarlo() but
 using the MISER algorithm.

******************************************************************************/

double integral_miser_mcarlo(const int n,
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

 Function example_a(): is an auxiliary function that returns the 3-dimensional
 integrand

 I(x, y, z) = A/[1 - cos(x)*cos(y)*cos(z)]; where, A = 1/pi*pi*pi.

 Useful to test Monte Carlo methods from a = (0, 0, 0) to b = (pi, pi, pi). The
 answer is 1.3932039296856768591842462603255.

******************************************************************************/

static inline double example_a(double x[], size_t n, void *params)
{
	ASSERT(n > 0)
	ASSERT(params == NULL)

	const double A = 1.0/(M_PI*M_PI*M_PI);

	return A/(1.0 - cos(x[0])*cos(x[1])*cos(x[2]));
}

#define EXAMPLEA_ANSWER 1.3932039296856768591842462603255

/******************************************************************************

 Function example_b(): is an auxiliary function that returns the integrand

 I(x) = exp(-x)*cos(x).

 Useful to test all integral methods from a = 0 to b = inf. The answer is 1/2.

******************************************************************************/

static inline double example_b(double x, void *params)
{
	ASSERT(params == NULL)
	return exp(-x)*cos(x);
}

#define EXAMPLEB_ANSWER 0.5

/******************************************************************************

 Function integral_benchmark(): return the error and wall time for a given
 method of integration.

******************************************************************************/

void integral_benchmark(const integral_method type,
                        const int n_max, double *error, double *wtime)
{
	switch (type)
	{
		case simpson_1st:
		{
			const double dx = 100.0/as_double(n_max);
			double *f = allocate(n_max, sizeof(double), false);

			for (int n = 0; n < n_max; ++n)
			{
				f[n] = example_b(as_double(n)*dx, NULL);
			}

			const double start_time = wall_time();
			const double result = integral_simpson(n_max, dx, false, f);
			const double end_time = wall_time();

			free(f);

			*wtime = end_time - start_time;
			*error = result - EXAMPLEB_ANSWER;
		}
		break;

		case simpson_2nd:
		{
			const double dx = 100.0/as_double(n_max);
			double *f = allocate(n_max, sizeof(double), false);

			for (int n = 0; n < n_max; ++n)
			{
				f[n] = example_b(as_double(n)*dx, NULL);
			}

			const double start_time = wall_time();
			const double result = integral_simpson_2nd(n_max, dx, false, f);
			const double end_time = wall_time();

			free(f);

			*wtime = end_time - start_time;
			*error = result - EXAMPLEB_ANSWER;
		}
		break;

		case plain_monte_carlo:
		{
			const double a[3] = {0.0, 0.0, 0.0};
			const double b[3] = {M_PI, M_PI, M_PI};

			const double start_time = wall_time();
			const double result = integral_plain_mcarlo(3, n_max, a, b, NULL, example_a);
			const double end_time = wall_time();

			*wtime = end_time - start_time;
			*error = result - EXAMPLEA_ANSWER;
		}
		break;

		case vegas_monte_carlo:
		{
			double a[3] = {0.0, 0.0, 0.0};
			double b[3] = {M_PI, M_PI, M_PI};

			const double start_time = wall_time();
			const double result = integral_vegas_mcarlo(3, n_max, a, b, NULL, example_a);
			const double end_time = wall_time();

			*wtime = end_time - start_time;
			*error = result - EXAMPLEA_ANSWER;
		}
		break;

		case miser_monte_carlo:
		{
			const double a[3] = {0.0, 0.0, 0.0};
			const double b[3] = {M_PI, M_PI, M_PI};

			const double start_time = wall_time();
			const double result = integral_miser_mcarlo(3, n_max, a, b, NULL, example_a);
			const double end_time = wall_time();

			*wtime = end_time - start_time;
			*error = result - EXAMPLEA_ANSWER;
		}
		break;

		default:
			PRINT_ERROR("invalid method %d\n", type)
			exit(EXIT_FAILURE);
	}
}
