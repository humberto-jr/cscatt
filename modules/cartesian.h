#if !defined(CARTESIAN_HEADER)
	#define CARTESIAN_HEADER
	#include "globals.h"
	#include "spherical.h"

	/******************************************************************************

	 Type cartesian: represents a set of Cartesian coordinates.

	******************************************************************************/

	struct cartesian
	{
		double x, y, z;
	};

	typedef struct cartesian cartesian;

	/******************************************************************************

	 Function cartesian_distance(): resolves the distance between two vectors, a =
	 (x, y, z) and b = (x', y', z'), in Cartesian coordinates.

	******************************************************************************/

	inline static double cartesian_distance(const cartesian *a, const cartesian *b)
	{
		return sqrt(pow(a->x - b->x, 2) + pow(a->y - b->y, 2) + pow(a->z - b->z, 2));
	}

	/******************************************************************************

	 Function cartesian_dot_prod(): resolves the dot product between two vectors, a
	 = (x, y, z) and b = (x', y', z'), in Cartesian coordinates.

	******************************************************************************/

	inline static double cartesian_dot_prod(const cartesian *a, const cartesian *b)
	{
		return a->x*b->x + a->y*b->y + a->z*b->z;
	}

	/******************************************************************************

	 Function cartesian_length(): resolves the length of a vector a = (x, y, z) in
	 Cartesian coordinates.

	******************************************************************************/

	inline static double cartesian_length(const cartesian *a)
	{
		return sqrt(a->x*a->x + a->y*a->y + a->z*a->z);
	}

	/******************************************************************************

	 Function cartesian_to_spherical(): resolves a set of spherical coordinates, b,
	 from Cartesian ones, a.

	 NOTE: outputted angles in degree.

	******************************************************************************/

	inline static void cartesian_to_spherical(const cartesian *a, spherical *b)
	{
		b->rho = cartesian_length(a);

		if (b->rho == 0.0)
		{
			b->theta = 0.0;
			b->phi = 0.0;
			return;
		}

		b->theta = acos(a->z/b->rho)*180.0/M_PI;

		if (a->x == 0.0 && a->y == 0.0)
			b->phi = 0.0;

		else if (a->x == 0.0 && a->y > 0.0)
			b->phi = 90.0;

		else if (a->x == 0.0 && a->y < 0.0)
			b->phi = 270.0;

		else
			b->phi = atan(a->y/a->x)*180.0/M_PI;
	}
#endif
