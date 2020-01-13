#if !defined(PHYS_HEADER)
	#define PHYS_HEADER
	#include "globals.h"

	double phys_wigner_3j(const int a, const int b, const int c,
	                      const int d, const int e, const int f);

	double phys_wigner_6j(const int a, const int b, const int c,
	                      const int d, const int e, const int f);

	double phys_percival_seaton(const int spin_mult, const int j1, const int j2,
	                            const int l1, const int l2, const int lambda, const int J);


	/******************************************************************************

	 Function phys_centr_term(): returns the actual centrifugal potential term at x
	 for a given angular momentum l, provided the mass under consideration.

	******************************************************************************/

	inline static double phys_centr_term(const int l,
	                                     const double mass, const double x)
	{
		return as_double(l*(l + 1))/(2.0*mass*x*x);
	}

	/******************************************************************************

	 Function phys_parity(): return the parity p = (-1)^l for a given l; either p =
	 +1 or p = -1.

	******************************************************************************/

	inline static int phys_parity(const int l)
	{
		return (int) pow(-1.0, l);
	}
#endif
