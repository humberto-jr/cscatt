#if !defined(MISC_HEADER)
	#define MISC_HEADER
	#include <math.h>

	#include "globals.h"

	/******************************************************************************

	 Function as_double(): a simple interface to a very common operation of casting
	 an integer, n, as a double precision number.

	******************************************************************************/

	inline static double as_double(const int n)
	{
		return (double) n;
	}

	/******************************************************************************

	 Function side_c(): uses the law of cosines to resolve the side c opposite to
	 the interior angle C (in degree) of a triangle with sides a, b and c.

	******************************************************************************/

	inline static double side_c(const double side_a,
	                            const double side_b, const double angle_c)
	{
		return sqrt(side_a*side_a + side_b*side_b - 2.0*side_a*side_b*cos(angle_c*M_PI/180.0));
	}

	/******************************************************************************

	 Function angle_c(): uses the law of cosines to resolve the angle C (in degree)
	 opposite to the side c of a triangle with sides a, b and c.

	******************************************************************************/

	inline static double angle_c(const double side_a,
	                             const double side_b, const double side_c)
	{
		return acos((side_c*side_c - side_a*side_a - side_b*side_b)/(-2.0*side_a*side_b))*180.0/M_PI;
	}

	/******************************************************************************

	 Function centr_term(): returns the actual centrifugal potential term at x for
	 a given angular momentum l, provided the mass under consideration.

	******************************************************************************/

	inline static double centr_term(const int l,
	                                const double mass, const double x)
	{
		return as_double(l*(l + 1))/(2.0*mass*x*x);
	}

	/******************************************************************************

	 Function sinc(): returns sin(x)/x for a given x.

	******************************************************************************/

	inline static double sinc(const double x)
	{
		return sin(x)/x;
	}

	/******************************************************************************

	 Function parity(): return the parity p = (-1)^l for a given l; either p = +1
	 or p = -1.

	******************************************************************************/

	inline static int parity(const int l)
	{
		return (int) pow(-1.0, l);
	}

	/******************************************************************************

	 Function sigmoid(): return the sigmoid function at x.

	******************************************************************************/

	inline static double sigmoid(const double x)
	{
		return 1.0/(1.0 + exp(-x));
	}
#endif
