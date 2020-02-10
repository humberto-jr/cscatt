#if !defined(JACOBI_HEADER)
	#define JACOBI_HEADER
	#include "globals.h"

	/******************************************************************************

	 Type jacobi: represent one set of Jacobi coordinates for triatomic systems
	 at a given arrangement indexed either by 'a' (A + BC), 'b' (B + CA) or 'c'
	 (C + AB).

	******************************************************************************/

	struct jacobi
	{
		char arrang;
		double r, R, theta;
	};

	typedef struct jacobi jacobi;

	struct jacobi_sf
	{
		char arrang;
		double r, theta, phi, mass_a, mass_b;
	};

	typedef struct jacobi_sf jacobi_sf;

	inline static void jacobi_to_cartesian(const jacobi_sf *a, cartesian *b)
	{
		b->x = a->r*sin(a->theta*M_PI/180.0)*cos(a->phi*M_PI/180.0);
		b->y = a->r*sin(a->theta*M_PI/180.0)*sin(a->phi*M_PI/180.0);
		b->z = a->r*cos(a->theta*M_PI/180.0);
	}
#endif
