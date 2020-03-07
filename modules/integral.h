#if !defined(INTEGRAL_HEADER)
	#define INTEGRAL_HEADER
	#include "globals.h"

	void integral_set_error(const double error);

	void integral_set_workspace(const int size);

	double integral_qag(const double a,
	                    const double b,
	                    void *params,
	                    double (*f)(double x, void *params));

	double integral_qags(const double a,
	                     const double b,
	                     void *params,
	                     double (*f)(double x, void *params));

	double integral_plain_mcarlo(const int n,
	                             const int calls,
	                             const double a[],
	                             const double b[],
	                             void *params,
	                             double (*f)(double x[], size_t n, void *params));

	double integral_vegas_mcarlo(const int n,
	                             const int calls,
	                             double a[],
	                             double b[],
	                             void *params,
	                             double (*f)(double x[], size_t n, void *params));

	double integral_miser_mcarlo(const int n,
	                             const int calls,
	                             const double a[],
	                             const double b[],
	                             void *params,
	                             double (*f)(double x[], size_t n, void *params));

	void integral_mcarlo_benchmark(const char type,
	                               const int calls, double *error, double *wtime);
#endif
