#if !defined(SPHERICAL_HEADER)
	#define SPHERICAL_HEADER
	#include "globals.h"
	#include "cartesian.h"

	/* TODO: external lib. headers shall go to .c implementation files. */
	#include <gsl/gsl_sf_legendre.h>

	/******************************************************************************

	 Type spherical: represents a set of spherical coordinates.

	******************************************************************************/

	struct spherical
	{
		double rho, theta, phi;
	};

	typedef struct spherical spherical;

	/******************************************************************************

	 Function spherical_to_cartesian(): resolves a set of Cartesian coordinates, b,
	 from spherical ones, a.

	******************************************************************************/

	inline static double spherical_to_cartesian(const spherical *a, cartesian *b)
	{
		b->x = a->rho*sin(a->theta*M_PI/180.0)*cos(a->phi*M_PI/180.0);
		b->y = a->rho*sin(a->theta*M_PI/180.0)*sin(a->phi*M_PI/180.0);
		b->z = a->rho*cos(a->theta*M_PI/180.0);
	}

	/******************************************************************************

	 Function spherical_harmonics(): returns the spherical harmonics for angular
	 momentum l, with projection m, at x.

	******************************************************************************/

	inline static double spherical_harmonics(const int l,
	                                         const int m, const double x)
	{
		ASSERT(m >= 0)
		ASSERT(l >= m)
		ASSERT(fabs(x) <= 1.0)

		return pow(-1.0, m)*gsl_sf_legendre_sphPlm(l, m, x);
	}
#endif
