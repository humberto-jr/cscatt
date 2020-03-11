#if !defined(JACOBI_HEADER)
	#define JACOBI_HEADER
	#include "globals.h"
	#include "cartesian.h"
	#include "spherical.h"

	#define MAX_JACOBI_VECTOR 3

	/******************************************************************************

	 Type jacobi_sf: represent one set of Jacobi coordinates for triatomic systems
	 at a given arrangement indexed either by 'a' (A + BC), 'b' (B + CA) or 'c'
	 (C + AB).

	******************************************************************************/

	struct jacobi_sf
	{
		int size;
		char arrang;
		double mass[MAX_JACOBI_VECTOR];
		spherical r[MAX_JACOBI_VECTOR];
	};

	typedef struct jacobi_sf jacobi_sf;

	void jacobi_to_cartesian(const jacobi_sf *set,
	                         cartesian com[], cartesian atom[]);

	void jacobi_to_internuc(const jacobi_sf *set, double r[]);
#endif
