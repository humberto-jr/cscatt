#if !defined(SPLINE_HEADER)
	#define SPLINE_HEADER

	typedef struct spline spline;

	spline *spline_alloc(const int grid_size,
	                     const double f[], const double x[], const char type);

	void spline_free(spline *s);

	double spline_value(const spline *s, const double x);

	double spline_derivative(const spline *s, const int order, const double x);

	double spline_integral(const spline *s, const double a, const double b);
#endif
