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
******************************************************************************/

void integral_set_error(const double error)
{
	ASSERT(error > 0.0)
	abs_error = error;
}

/******************************************************************************
******************************************************************************/

void integral_set_workspace(const int size)
{
	ASSERT(size > 0)
	workspace_size = size;
}

/******************************************************************************
******************************************************************************/

double integral_qag(const double x_min,
                    const double x_max, void *params, qag_integrand *f)
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

	const int info = gsl_integration_qag(&gsl_f, x_min, x_max, abs_error, 0.0,
	                                     workspace_size, GSL_INTEG_GAUSS61, work, &result, &error);

	gsl_integration_workspace_free(work);

	if (info != 0)
	{
		PRINT_ERROR("%s; error = %f\n", gsl_strerror(info), error)
	}

	return result;
}

/******************************************************************************
******************************************************************************/

double integral_qags(const double x_min,
                     const double x_max, void *params, qag_integrand *f)
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

	const int info = gsl_integration_qags(&gsl_f, x_min, x_max, abs_error, 0.0,
	                                      workspace_size, work, &result, &error);

	gsl_integration_workspace_free(work);

	if (info != 0)
	{
		PRINT_ERROR("%s; error = %f\n", gsl_strerror(info), error)
	}

	return result;
}

/******************************************************************************
******************************************************************************/

double integral_mc_plain(const int dim,
                         const double x_min[],
                         const double x_max[],
                         const int max_call,
                         void *params,
                         mc_integrand *f)
{
	gsl_monte_function gsl_f =
	{
		.params = params,
		.dim = dim,
		.f = f
	};

	gsl_monte_plain_state *work
		= gsl_monte_plain_alloc(workspace_size);

	ASSERT(work != NULL)

	gsl_rng *r
		= gsl_rng_alloc(gsl_rng_ranlxd2);

	ASSERT(r != NULL)

	double result = 0.0, error = 0.0;

	const int info = gsl_monte_plain_integrate(&gsl_f, x_min, x_max, dim,
	                                           max_call, r, work, &result, &error);

	gsl_rng_free(r);
	gsl_monte_plain_free(work);

	if (info != 0)
	{
		PRINT_ERROR("%s; error = %f\n", gsl_strerror(info), error)
	}

	return result;
}

/******************************************************************************
******************************************************************************/

double integral_mc_vegas(const int dim,
                         double x_min[],
                         double x_max[],
                         const int max_call,
                         void *params,
                         mc_integrand *f)
{
	gsl_monte_function gsl_f =
	{
		.params = params,
		.dim = dim,
		.f = f
	};

	gsl_monte_vegas_state *work
		= gsl_monte_vegas_alloc(workspace_size);

	ASSERT(work != NULL)

	gsl_rng *r
		= gsl_rng_alloc(gsl_rng_ranlxd2);

	ASSERT(r != NULL)

	double result = 0.0, error = 0.0;

	const int info = gsl_monte_vegas_integrate(&gsl_f, x_min, x_max, dim,
	                                           max_call, r, work, &result, &error);

	gsl_rng_free(r);
	gsl_monte_vegas_free(work);

	if (info != 0)
	{
		PRINT_ERROR("%s; error = %f\n", gsl_strerror(info), error)
	}

	return result;
}

/******************************************************************************
******************************************************************************/

double integral_mc_miser(const int dim,
                         const double x_min[],
                         const double x_max[],
                         const int max_call,
                         void *params,
                         mc_integrand *f)
{
	gsl_monte_function gsl_f =
	{
		.params = params,
		.dim = dim,
		.f = f
	};

	gsl_monte_miser_state *work
		= gsl_monte_miser_alloc(workspace_size);

	ASSERT(work != NULL)

	gsl_rng *r
		= gsl_rng_alloc(gsl_rng_ranlxd2);

	ASSERT(r != NULL)

	double result = 0.0, error = 0.0;

	const int info = gsl_monte_miser_integrate(&gsl_f, x_min, x_max, dim,
	                                           max_call, r, work, &result, &error);

	gsl_rng_free(r);
	gsl_monte_miser_free(work);

	if (info != 0)
	{
		PRINT_ERROR("%s; error = %f\n", gsl_strerror(info), error)
	}

	return result;
}
