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

	struct jacobi_6d
	{
		double r1, r2, r3, theta12, theta3, phi3;
	};

	typedef struct jacobi_6d jacobi_6d;

	/******************************************************************************

	 Type internuc_coor: represent one set of internuclear distances for triatomic
	 systems; where, r is the respective distance between atoms A, B and C.

	******************************************************************************/

	struct internuc_coor
	{
		double r_ab, r_bc, r_ac;
	};

	typedef struct internuc_coor internuc_coor;

	struct internuc
	{
		double ab, bc, ac;
	};

	typedef struct internuc internuc;

	void pes_init_mass(FILE *input, const char atom);

	double pes_mass(const char atom);

	double pes_abc(const char arrang,
	               const double r, const double R, const double theta);

	double pes_bc(const int j, const double r);

	double pes_ac(const int j, const double r);

	double pes_ab(const int j, const double r);

	double pes(const jacobi_coor *x);

	void pes_init();

	double pec(const char arrang, const double r);

	matrix *pes_olson_smith_model(const double x);

	matrix *pes_tully_model(const int n, const double x);

	void pes_about(FILE *output);
#endif
