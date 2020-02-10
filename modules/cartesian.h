#if !defined(CARTESIAN_HEADER)
	#define CARTESIAN_HEADER
	#include "globals.h"

	/******************************************************************************

	 Type cartesian: represent one set of Cartesian coordinates.

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

	 Function cartesian_to_jacobi(): resolves space-fixed Jacobi coordinates from
	 Cartesian coordinates.

	******************************************************************************/

	inline static void cartesian_to_jacobi(const cartesian *a, jacobi_sfixed *b)
	{
		const double c = sqrt(a.x*a.x + a.y*a.y + a.z*a.z);

		if (c == 0.0)
		{
			b->theta = 0.0;
			b->phi = 90.0;
			b->r = 0.0;
			return;
		}

		b->r = c;
		b->theta = acos(a.z/c);

		if (a.x > 0.0)
			b->phi = atan(a.y/a.x);
		else if (a.x < 0.0 && a.y >= 0.0)
			b->phi = atan(a.y/a.x) + M_PI;
		else if (a.x < 0.0 && a.y < 0.0)
			b->phi = atan(a.y/a.x) - M_PI;
		else if (a.x == 0.0 && a.y > 0.0)
			b->phi = M_PI/2.0;
		else if (a.x == 0.0 && a.y < 0.0)
			b->phi = -M_PI/2.0;
		else if (a.x == 0.0 && a.y == 0.0)
			b->phi = -M_PI/2.0;

		b->theta *= 180.0/M_PI;
		b->phi *= 180.0/M_PI;
	}
#endif
