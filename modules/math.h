#if !defined(MATH_HEADER)
	#define MATH_HEADER
	#include "globals.h"

	#if !defined(MATH_PI)
		#define MATH_PI 3.14159265358979323846264338328
	#endif

	#if !defined(MATH_SQRT_2PI)
		#define MATH_SQRT_2PI 2.5066282746310005024157652848110452530069867406099383166299
	#endif

	/******************************************************************************

	 Type math_xyz: represents a set of Cartesian coordinates.

	******************************************************************************/

	struct math_xyz
	{
		double x, y, z;
	};

	typedef struct math_xyz math_xyz;

	double math_distance(const math_xyz *a, const math_xyz *b);

	double math_dot_product(const math_xyz *a, const math_xyz *b);

	double math_length(const math_xyz *a);

	void math_spherical_coor(const math_xyz *a,
	                         double *rho, double *theta, double *phi);

	void math_cartesian_coor(math_xyz *a,
	                         const double rho, const double theta, const double phi);

	void math_euler_rotation(math_xyz *a,
	                         const double alpha,
	                         const double beta,
	                         const double gamma);

	double math_legendre_poly(const size_t l, const double x);

	double math_assoc_legendre_poly(const size_t l,
	                                const size_t m, const double x);

	double complex math_sphe_harmonics(const size_t l, const int m,
	                                   const double theta, const double phi);

	double math_wigner_3j(const int a, const int b, const int c,
	                      const int d, const int e, const int f);

	double math_wigner_6j(const int a, const int b, const int c,
	                      const int d, const int e, const int f);

	double math_wigner_9j(const int a, const int b, const int c,
	                      const int d, const int e, const int f,
	                      const int g, const int h, const int i);

	double math_sphe_bessel(const char type, const size_t l, const double x);

	double math_clebsch_gordan(const int j1, const int j2, const int j3,
	                           const int m1, const int m2, const int m3);

	double math_racah_coeff(const int a,
	                        const int b,
	                        const int c,
	                        const int d,
	                        const int e,
	                        const int f);

	double math_percival_seaton(const int J,
	                            const int j1,
	                            const int j2,
	                            const int l1,
	                            const int l2,
	                            const int lambda);

	double math_gaunt(const int k,
	                  const int j1,
	                  const int j2,
	                  const int lambda);

	double *math_wigner_d(const double k_in,
	                      const double m_in,
	                      const double j_in,
	                      const double beta);

	double math_wigner_d2(const double n,
	                      const double m,
	                      const double j,
	                      const double beta);

	double math_integral_yyy(const int j1, const int j2, const int j3,
	                         const int m1, const int m2, const int m3);

	void math_set_workspace(const int size);

	void math_set_error(const double error);

	void math_no_gsl_handler();

	double math_sinc(const double x);

	double math_sigmoid(const double x);

	double math_side_c(const double side_a,
	                   const double side_b, const double angle_c);

	double math_angle_c(const double side_a,
	                    const double side_b, const double side_c);

	double math_simpson(const double a,
	                    const double b,
	                    const size_t n_max,
	                    const void *params,
	                    const bool use_omp,
	                    double (*f)(const double x, const void *params));

	double math_simpson_array(const double a,
	                          const double b,
	                          const size_t n_max,
	                          const bool use_omp,
	                          const double array[]);

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

	size_t *math_bubble_sort(const size_t n_max, const double x[]);

	double math_factorial(const int n);

	double math_gauss_legendre(const double a,
	                           const double b,
	                           const int order,
	                           void *params,
	                           double (*f)(const double x, void *params));

	void math_about(FILE *output);
#endif
