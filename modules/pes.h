#if !defined(PES_HEADER)
	#define PES_HEADER
	#include "globals.h"
	#include "matrix.h"

	#define MAX_JACOBI_VECTOR 5
	#define MAX_ATOM (MAX_JACOBI_VECTOR + 1)
	#define MAX_INTERNUC_DISTANCE (3*MAX_ATOM - 6)

	/******************************************************************************

	 Type jacobi_coor: represent one set of Jacobi coordinates for triatomic systems
	 at a given arrangement indexed either by 'a' (A + BC), 'b' (B + CA) or 'c'
	 (C + AB).

	******************************************************************************/

	struct jacobi_coor
	{
		char arrang;
		double r, R, theta;
	};

	typedef struct jacobi_coor jacobi_coor;

	/******************************************************************************

	 Type internuc_coor: represent one set of internuclear distances for triatomic
	 systems; where, r is the respective distance between atoms A, B and C.

	******************************************************************************/

	struct internuc_coor
	{
		double r_ab, r_bc, r_ac;
	};

	typedef struct internuc_coor internuc_coor;

	double pes(const jacobi_coor *x);

//	double pes_sf(const jacobi_sf *set);

	void pes_init();

	double pec(const char arrang, const double r);

//	double pes_asymptotic_min(const char arrang, const double scan_step);

	matrix *pes_olson_smith_model(const double x);

	matrix *pes_tully_model(const int n, const double x);

	void pes_about(FILE *output);
#endif
