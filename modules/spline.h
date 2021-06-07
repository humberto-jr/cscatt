#if !defined(SPLINE_HEADER)
	#define SPLINE_HEADER
	#include "globals.h"

	typedef struct spline spline;

	struct spline_handle
	{
		spline *s;
		double *x, *f;
	};

	typedef struct spline_handle spline_handle;

	spline *spline_alloc(const size_t grid_size,
	                     const double x[], const double f[], const char type);

	void spline_free(spline *s);

	double spline_value(const spline *s, const double x);

	double spline_derivative(const spline *s, const size_t order, const double x);

	double spline_integral(const spline *s, const double a, const double b);

	void spline_about(FILE *output);
#endif
