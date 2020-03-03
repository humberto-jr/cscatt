#if !defined(INTEGRAL_HEADER)
	#define INTEGRAL_HEADER
	#include "globals.h"

	void integral_set_error(const double error);

	void integral_set_workspace(const int size);

	double integral_qag(const double x_min,
	                    const double x_max,
	                    void *params,
	                    double (*integrand)(double x, void *params));

	double integral_qags(const double x_min,
	                     const double x_max,
	                     void *params,
	                     double (*integrand)(double x, void *params));

	double integral_mc_plain(const int n,
	                         const double x_min[],
	                         const double x_max[],
	                         const int max_call,
	                         void *params,
	                         double (*integrand)(double x[], size_t n, void *params));

	double integral_mc_vegas(const int n,
	                         double x_min[],
	                         double x_max[],
	                         const int max_call,
	                         void *params,
	                         double (*integrand)(double x[], size_t n, void *params));

	double integral_mc_miser(const int n,
	                         const double x_min[],
	                         const double x_max[],
	                         const int max_call,
	                         void *params,
	                         double (*integrand)(double x[], size_t n, void *params));
#endif
