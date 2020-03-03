#if !defined(INTEGRAL_HEADER)
	#define INTEGRAL_HEADER
	#include "globals.h"

	typedef double (*qag_integrand)(double x, void *params);

	typedef double (*mc_integrand)(double x[], int dim, void *params);

	void integral_set_error(const double error);

	void integral_set_workspace(const int size);

	double integral_qag(const double x_min,
	                    const double x_max, void *params, qag_integrand *f);

	double integral_qags(const double x_min,
	                     const double x_max, void *params, qag_integrand *f);

	double integral_mc_plain(const int dim,
	                         const double x_min[],
	                         const double x_max[],
	                         const int max_call,
	                         void *params,
	                         mc_integrand *f);

	double integral_mc_vegas(const int dim,
	                         double x_min[],
	                         double x_max[],
	                         const int max_call,
	                         void *params,
	                         mc_integrand *f);

	double integral_mc_miser(const int dim,
	                         const double x_min[],
	                         const double x_max[],
	                         const int max_call,
	                         void *params,
	                         mc_integrand *f);
#endif
