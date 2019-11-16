#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "globals.h"
#include "coor.h"
#include "pes.h"

#if !defined(EXTERNAL_PES_NAME)
	#define EXTERNAL_PES_NAME pes_noname

	inline static double EXTERNAL_PES_NAME(const double *r_ab,
	                                       const double *r_bc,
	                                       const double *r_ac)
	{
		ASSERT(r_ab != NULL)
		ASSERT(r_bc != NULL)
		ASSERT(r_ac != NULL)

		PRINT_ERROR("%s\n", "no external PES defined")
		exit(EXIT_FAILURE);
	}
#else
	extern double EXTERNAL_PES_NAME(const double *r_ab,
	                                const double *r_bc,
	                                const double *r_ac);
#endif

/******************************************************************************

 Wrapper pes(): is an interface for the external user defined PES named pot_(),
 as function of a set of Jacobi coordinates, x = (r, R, theta), which are trans
 lated to internuclear distances, y = (r_ab, r_bc, r_ac).

 NOTE: If the macro USE_NON_REACTIVE_PES is defined, the coordinate system is
 not converted and Jacobi is used all along, assuming

 pot_(r, R, theta)

 Thus, make sure the order of each input parameter of pot_() follows the inter
 face given above.

******************************************************************************/

double pes(const jacobi_coor *x)
{
	#if defined(USE_NON_REACTIVE_PES)
		return pot_(&x->r, &x->R, &x->theta);
	#else
		internuc_coor y;
		coor_jacobi_to_internuc(x, &y);
		return pot_(&y.r_ab, &y.r_bc, &y.r_ac);
	#endif
}

/******************************************************************************

 Function pec(): returns the diatomic potential, V(r) as function of the inter
 nuclear distance, r, asymptotically as R -> inf, for a given PES arrangement.

******************************************************************************/

double pec(const char arrang, const double r)
{
	const jacobi_coor x =
	{
		.r = r,
		.R = INF,
		.theta = 90.0,
		.arrang = arrang
	};

	return pes(&x);
}

/******************************************************************************

 Function pes_about(): print in a given output file the conditions in which the
 module was compiled.

******************************************************************************/

void pes_about(FILE *output)
{
	ASSERT(output != NULL)

	fprintf(output, "# build date   = %s\n", __DATE__);
	fprintf(output, "# source code  = %s\n", __FILE__);
	fprintf(output, "# external PES = %s\n", PRINT_MACRO(EXTERNAL_PES_NAME));

	#if defined(USE_NON_REACTIVE_PES)
		fprintf(output, "# coordinates  = Jacobi\n");
		fprintf(output, "# interface    = %s(r, R, theta)\n", PRINT_MACRO(EXTERNAL_PES_NAME));
	#else
		fprintf(output, "# coordinates  = internuclear distances\n");
		fprintf(output, "# interface    = %s(r_ab, r_bc, r_ac)\n", PRINT_MACRO(EXTERNAL_PES_NAME));
	#endif
}
