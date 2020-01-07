#if !defined(MASS_HEADER)
	#define MASS_HEADER
	#include <stdio.h>

	/******************************************************************************

	 Type mass_case: is an enumeration of several types of masses for a triatomic
	 system (A, B and C).

	 NOTE: See function mass().

	******************************************************************************/

	enum mass_case
	{
		atom_a,
		atom_b,
		atom_c,
		pair_ab,
		pair_ac,
		pair_bc,
		arrang_a,
		arrang_b,
		arrang_c,
		total
	};

	void mass_init(FILE *input);

	double mass(const enum mass_case m);
#endif
