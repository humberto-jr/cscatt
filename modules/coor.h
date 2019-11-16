#if !defined(COOR_HEADER)
	#define COOR_HEADER

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

	/******************************************************************************

	 Type xyz_coor: represent one set of Cartesian coordinates.

	******************************************************************************/

	struct xyz_coor
	{
		double x, y, z;
	};

	typedef struct xyz_coor xyz_coor;

	/******************************************************************************

	 Function distance(): resolves the distance between two points, a = (x, y, z)
	 and b = (x', y', z'), in Cartesian coordinates.

	******************************************************************************/

	inline static double coor_xyz_distance(const xyz_coor *a, const xyz_coor *b)
	{
		return sqrt(pow(a->x - b->x, 2) + pow(a->y - b->y, 2) + pow(a->z - b->z, 2));
	}

	void coor_jacobi_to_internuc(const jacobi_coor *from, internuc_coor *to);
#endif
