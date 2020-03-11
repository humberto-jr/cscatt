#if !defined(SPHERICAL_HEADER)
	#define SPHERICAL_HEADER
	#include "globals.h"

	/******************************************************************************

	 Type spherical: represents a set of spherical coordinates.

	******************************************************************************/

	struct spherical
	{
		double rho, theta, phi;
	};

	typedef struct spherical spherical;

	double spherical_harmonics(const int l,
	                           const int m,
	                           const spherical *r);

	double spherical_harmonics_coupl(const int n,
	                                 const int l[],
	                                 const int m[],
	                                 const spherical r[]);
#endif
