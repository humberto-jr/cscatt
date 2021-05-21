/******************************************************************************

 About
 -----

 This module defines the opaque type spline, for interpolations of numerical
 functions f = f(x), and few general purpose routines needed to manipulate it.

******************************************************************************/

#include "spline.h"
#include "gsl_lib.h"

struct spline
{
	size_t grid_size;
	const double *x, *f;

	gsl_interp *data;
	gsl_interp_accel *state;
};

/******************************************************************************

 Function spline_alloc(): allocate resources for a spline interpolation of f(x)
 for a given set of tabulated values of f and x. Where, type is one of: 'a' for
 Akima spline, 'c' for cubic spline or 's' for Steffen spline.

 NOTE: internal pointers pointing to x and f are used, thus they should not be
 freed while the spline is used.

******************************************************************************/

spline *spline_alloc(const size_t grid_size,
                     const double x[], const double f[], const char type)
{
	ASSERT(x != NULL)
	ASSERT(f != NULL)
	ASSERT(grid_size > 6)

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

	s->x = x;
	s->f = f;
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

double spline_derivative(const spline *s, const size_t order, const double x)
{
	ASSERT(s != NULL)

	switch (order)
	{
		case 1: return gsl_interp_eval_deriv(s->data, s->x, s->f, x, s->state);
		case 2: return gsl_interp_eval_deriv2(s->data, s->x, s->f, x, s->state);

		default:
			PRINT_ERROR("invalid order %d\n", (int) order)
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
