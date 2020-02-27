#if !defined(SPHERICAL_HEADER)
	#define SPHERICAL_HEADER
	#include "globals.h"
	#include "cartesian.h"

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
#endif
