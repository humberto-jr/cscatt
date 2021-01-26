/******************************************************************************

 About
 -----

 This module defines the opaque type spline, for interpolations of numerical
 functions f = f(x), and few general purpose routines needed to manipulate it.

******************************************************************************/

#include "spline.h"
#include "gsl_lib.h"
#include "globals.h"

struct spline
{
	int grid_size;
	double *x, *f;
	gsl_interp *data;
	gsl_interp_accel *state;
};

/******************************************************************************

 Function spline_alloc(): allocate resources for a spline interpolation of f(x)
 for a given set of tabulated values of f and x. Where, type is one of: 'a' for
 Akima spline, 'c' for cubic spline or 's' for Steffen spline.

******************************************************************************/

spline *spline_alloc(const int grid_size,
                     const double f[], const double x[], const char type)
{
	ASSERT(f != NULL)
	ASSERT(x != NULL)
	ASSERT(grid_size > 0)

	spline *s = allocate(1, sizeof(struct spline), true);

	switch (type)
	{
		case 'a':
			s->data = gsl_interp_alloc(gsl_interp_akima, grid_size);
			break;

		case 'c':
			s->data = gsl_interp_alloc(gsl_interp_cspline, grid_size);
			break;

		case 's':
			s->data = gsl_interp_alloc(gsl_interp_steffen, grid_size);
			break;

		default:
			PRINT_ERROR("invalid type %c\n", type)
			exit(EXIT_FAILURE);
	}

	ASSERT(s->data != NULL)

	s->state = gsl_interp_accel_alloc();

	ASSERT(s->state != NULL)

	s->x = allocate(grid_size, sizeof(double), false);
	s->f = allocate(grid_size, sizeof(double), false);

	for (int n = 0; n < grid_size; ++n)
	{
		s->x[n] = x[n];
		s->f[n] = f[n];
	}

	s->grid_size = grid_size;

	const int info = gsl_interp_init(s->data, s->x, s->f, s->grid_size);

	if (info != 0)
	{
		PRINT_ERROR("%s\n", gsl_strerror(info))
		exit(EXIT_FAILURE);
	}

	return s;
}

/******************************************************************************

 Function spline_free(): release resources allocated by spline_alloc().

******************************************************************************/

void spline_free(spline *s)
{
	ASSERT(s != NULL)

	gsl_interp_accel_free(s->state);
	gsl_interp_free(s->data);
	free(s->f);
	free(s->x);
	free(s);
}

/******************************************************************************

 Function spline_value(): return the interpolated value of f, in s, for a given
 point x.

******************************************************************************/

double spline_value(const spline *s, const double x)
{
	ASSERT(s != NULL)

	return gsl_interp_eval(s->data, s->x, s->f, x, s->state);
}

/******************************************************************************

 Function spline_derivative(): return the first (order = 1) or second (order =
 2) derivative of a spline interpolated function, s, for a given point x.

******************************************************************************/

double spline_derivative(const spline *s, const int order, const double x)
{
	ASSERT(s != NULL)

	switch (order)
	{
		case 1: return gsl_interp_eval_deriv(s->data, s->x, s->f, x, s->state);
		case 2: return gsl_interp_eval_deriv2(s->data, s->x, s->f, x, s->state);

		default:
			PRINT_ERROR("invalid order %d\n", order)
			exit(EXIT_FAILURE);
	}
}

/******************************************************************************

 Function spline_integral(): return the numerical integral result of a spline
 interpolated function, s, over the range [a, b].

******************************************************************************/

double spline_integral(const spline *s, const double a, const double b)
{
	ASSERT(s != NULL)

	return gsl_interp_eval_integ(s->data, s->x, s->f, a, b, s->state);
}

/******************************************************************************

 Function spline_about(): print in a given output file the conditions in which
 the module was compiled.

******************************************************************************/

void spline_about(FILE *output)
{
	ASSERT(output != NULL)

	fprintf(output, "# build date  = %s\n", __DATE__);
	fprintf(output, "# source code = %s\n", __FILE__);
}
