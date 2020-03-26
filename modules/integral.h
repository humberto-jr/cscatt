#if !defined(INTEGRAL_HEADER)
	#define INTEGRAL_HEADER
	#include "globals.h"

	enum integral_method
	{
		qag,
		qags,
		simpson_1st,
		simpson_2nd,
		plain_monte_carlo,
		vegas_monte_carlo,
		miser_monte_carlo
	};

	typedef enum integral_method integral_method;

	void integral_set_error(const double error);

	void integral_set_workspace(const int size);

	double integral_simpson(const int n,
	                        const double a,
	                        const double b,
	                        void *params,
	                        const bool use_omp,
	                        double (*f)(double x, void *params));

	double integral_simpson_2nd(const int n,
	                            const double a,
	                            const double b,
	                            void *params,
	                            const bool use_omp,
	                            double (*f)(double x, void *params));

	double integral_tab_simpson(const int grid_size,
	                            const double grid_step,
	                            const bool use_omp,
	                            const double f[]);

	double integral_tab_simpson_2nd(const int grid_size,
	                                const double grid_step,
	                                const bool use_omp,
	                                const double f[]);

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

	void integral_benchmark(const integral_method type,
	                        const int n_max, double *error, double *wtime);
#endif
