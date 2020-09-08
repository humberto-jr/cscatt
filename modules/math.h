#if !defined(MATH_HEADER)
	#define MATH_HEADER
	#include "globals.h"

	struct math_lanczos_setup
	{
		void *params;
		int n, n_max, max_step;
		double *start_vector, *eigenval, *eigenvec;
		void (*product)(const int, const double *, double *, void *);
	};

	typedef struct math_lanczos_setup math_lanczos_setup;

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

	double math_racah_coeff(const int J,
	                        const int j1,
	                        const int j2,
	                        const int l1,
	                        const int l2,
	                        const int lambda);

	double math_percival_seaton(const int J,
	                            const int j1,
	                            const int j2,
	                            const int l1,
	                            const int l2,
	                            const int lambda);

	double *math_wigner_d(const int k,
	                      const int m,
	                      const int j_max,
	                      const double beta);

	double math_integral_yyy(const int j1, const int m1,
	                         const int j2, const int m2,
	                         const int j3, const int m3);

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

	void math_lanczos(math_lanczos_setup *s);

	/******************************************************************************

	 Function math_side_c() uses the law of cosines to resolve the side c opposite
	 to the interior angle C (in degrees) of a triangle with sides a, b and c.

	******************************************************************************/

	inline static double math_side_c(const double side_a,
	                                 const double side_b, const double angle_c)
	{
		return sqrt(side_a*side_a + side_b*side_b - 2.0*side_a*side_b*cos(angle_c*M_PI/180.0));
	}

	/******************************************************************************

	 Function math_angle_c() uses the law of cosines to resolve the angle C (in
	 degrees) opposite to the side c of a triangle with sides a, b and c.

	******************************************************************************/

	inline static double math_angle_c(const double side_a,
	                                  const double side_b, const double side_c)
	{
		return acos((side_c*side_c - side_a*side_a - side_b*side_b)/(-2.0*side_a*side_b))*180.0/M_PI;
	}
#endif
