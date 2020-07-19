#if !defined(MATH_HEADER)
	#define MATH_HEADER
	#include "globals.h"

	double math_legendre_poly(const int l, const double x);

	double math_sphe_harmonics(const int l, const int m,
	                           const double theta, const double phi);

	double math_wigner_3j(const int a, const int b, const int c,
	                      const int d, const int e, const int f);

	double math_wigner_6j(const int a, const int b, const int c,
	                      const int d, const int e, const int f);

	double math_wigner_9j(const int a, const int b, const int c,
	                      const int d, const int e, const int f,
	                      const int g, const int h, const int i);

	double math_clebsch_gordan(const int j1, const int j2, const int j3,
	                           const int m1, const int m2, const int m3);

	double *math_wigner_d(const int k,
	                      const int m,
	                      const int j_max,
	                      const double beta);

	void math_set_workspace(const int size);

	void math_set_error(const double error);

	void math_no_gsl_handler();

	double math_simpson(const int n,
	                    const double a,
	                    const double b,
	                    void *params,
	                    const bool use_omp,
	                    double (*f)(double x, void *params));

	double math_2nd_simpson(const int n,
	                        const double a,
	                        const double b,
	                        void *params,
	                        const bool use_omp,
	                        double (*f)(double x, void *params));

	double math_qag(const double a,
	                const double b,
	                void *params,
	                double (*f)(double x, void *params));

	double math_qags(const double a,
	                 const double b,
	                 void *params,
	                 double (*f)(double x, void *params));

	double math_plain_mcarlo(const int n,
	                         const int calls,
	                         const double a[],
	                         const double b[],
	                         void *params,
	                         double (*f)(double x[], size_t n, void *params));

	double math_vegas_mcarlo(const int n,
	                         const int calls,
	                         double a[],
	                         double b[],
	                         void *params,
	                         double (*f)(double x[], size_t n, void *params));

	double math_miser_mcarlo(const int n,
	                         const int calls,
	                         const double a[],
	                         const double b[],
	                         void *params,
	                         double (*f)(double x[], size_t n, void *params));
#endif
