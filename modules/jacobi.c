#include "jacobi.h"

/******************************************************************************

 Function jacobi_to_cartesian(): returns the n cartesian coordinates for the
 center-of-masses (COM), and n + 1 cartesian coordinates for the respective
 atoms, from a set of n Jacobi vectors in the space-fixed (SF) reference frame.
 Where, n = set->size.

******************************************************************************/

void jacobi_to_cartesian(const jacobi_sf *set, cartesian com[], cartesian atom[])
{
	ASSERT(set->size >= 2)
	ASSERT(set->size <= MAX_JACOBI_VECTOR)

	atom[0].x = 0.0;
	atom[0].y = 0.0;
	atom[0].z = 0.0;

	cartesian_from_spherical(&set->r[0], &atom[1]);

	double m = set->mass[0] + set->mass[1];

	com[0].x = set->mass[1]*atom[1].x/m;
	com[0].y = set->mass[1]*atom[1].y/m;
	com[0].z = set->mass[1]*atom[1].z/m;

	for (int n = 1; n < set->size; ++n)
	{
		cartesian_from_spherical(&set->r[n], &atom[n + 1]);

		atom[n + 1].x -= com[n - 1].x;
		atom[n + 1].y -= com[n - 1].y;
		atom[n + 1].z -= com[n - 1].z;

		const double com_mass = m;
		m += set->mass[n + 1];

		com[n].x = (com_mass*com[n - 1].x + set->mass[n + 1]*atom[n + 1].x)/m;
		com[n].y = (com_mass*com[n - 1].y + set->mass[n + 1]*atom[n + 1].y)/m;
		com[n].z = (com_mass*com[n - 1].z + set->mass[n + 1]*atom[n + 1].z)/m;
	}
}

/******************************************************************************

 Function jacobi_to_internuc(): returns the 3*n - 6 internuclear distances, r,
 from a set of n - 1 Jacobi vectors in the space-fixed (SF) reference frame.
 Where, n = set->size + 1.

******************************************************************************/

void jacobi_to_internuc(const jacobi_sf *set, double r[])
{
	cartesian com[MAX_JACOBI_VECTOR], atom[MAX_JACOBI_VECTOR + 1];

	jacobi_to_cartesian(set, com, atom);
	const int n_max = set->size + 1;

	int pair = 0;
	for (int n = 0; n < n_max; ++n)
	{
		for (int m = (n + 1); m < n_max; ++m)
		{
			r[pair] = cartesian_distance(&atom[n], &atom[m]);
			++pair;
		}
	}

	ASSERT(pair == 3*n_max - 6)
}
