#if !defined(GSL_LIB_HEADER)
	#define GSL_LIB_HEADER
	#include "c_lib.h"

	#if !defined(HAVE_INLINE)
		#define HAVE_INLINE
	#endif

	#if !defined(__GNUC__) && !defined(GSL_C99_INLINE)
		#define GSL_C99_INLINE
	#endif

	/* error handling */
	#include <gsl/gsl_errno.h>

	/* basic functions */
	#include <gsl/gsl_math.h>

	/* special functions */
	#include <gsl/gsl_sf_bessel.h>
	#include <gsl/gsl_sf_legendre.h>
	#include <gsl/gsl_sf_coupling.h>

	/* random numbers generation */
	#include <gsl/gsl_rng.h>

	/* matrices and vectors */
	#include <gsl/gsl_vector.h>
	#include <gsl/gsl_matrix.h>

	/* blas and linear algebra */
	#include <gsl/gsl_cblas.h>
	#include <gsl/gsl_eigen.h>
	#include <gsl/gsl_linalg.h>

	/* Monte Carlo integration */
	#include <gsl/gsl_monte_plain.h>
	#include <gsl/gsl_monte_miser.h>
	#include <gsl/gsl_monte_vegas.h>

	/* interpolation */
	#include <gsl/gsl_interp.h>
	#include <gsl/gsl_spline.h>

	/* numerical integration */
	#include <gsl/gsl_integration.h>
#endif
