/******************************************************************************

 About
 -----

 This module is a collection of routines and types to handle various sets of
 physical coordinates.

******************************************************************************/

#include "mass.h"
#include "coor.h"

/******************************************************************************

 Function coor_jacobi_to_internuc(): converts from a set of Jacobi coordinates
 to a set of internuclear distances for triatomic systems. The following order
 of arrangements is used: arrangement index 'a' for A + BC, 'b' for B + CA and
 'c' is C + AB. Where, A, B and C represents three atoms.

******************************************************************************/

void coor_jacobi_to_internuc(const jacobi_coor *from, internuc_coor *to)
{
	xyz_coor a, b, c;

	switch (from->arrang)
	{
		case 'a':
		{
			c.x = 0.0;
			c.y = from->r/2.0;
			c.z = 0.0;

			b.x =  0.0;
			b.y = -c.y;
			b.z =  0.0;

			const double mass_c = mass(atom_c);
			const double mass_b = mass(atom_b);
			const double cb_com = (c.y*mass_c + b.y*mass_b)/(mass_c + mass_b);

			a.x = 0.0;
			a.y = cb_com + from->R*sin(from->theta*M_PI/180.0);
			a.z = from->R*cos(from->theta*M_PI/180.0);

			to->r_bc = from->r;
			to->r_ac = coor_xyz_distance(&a, &c);
			to->r_ab = coor_xyz_distance(&a, &b);
		}
		break;

		case 'b':
		{
			c.x = 0.0;
			c.y = from->r/2.0;
			c.z = 0.0;

			a.x =  0.0;
			a.y = -c.y;
			a.z =  0.0;

			const double mass_c = mass(atom_c);
			const double mass_a = mass(atom_a);
			const double ca_com = (c.y*mass_c + a.y*mass_a)/(mass_c + mass_a);

			b.x = 0.0;
			b.y = ca_com + from->R*sin(from->theta*M_PI/180.0);
			b.z = from->R*cos(from->theta*M_PI/180.0);

			to->r_ac = from->r;
			to->r_bc = coor_xyz_distance(&b, &c);
			to->r_ab = coor_xyz_distance(&a, &b);
		}
		break;

		case 'c':
		{
			a.x = 0.0;
			a.y = from->r/2.0;
			a.z = 0.0;

			b.x =  0.0;
			b.y = -a.y;
			b.z =  0.0;

			const double mass_a = mass(atom_a);
			const double mass_b = mass(atom_b);
			const double ab_com = (a.y*mass_a + b.y*mass_b)/(mass_a + mass_b);

			c.x = 0.0;
			c.y = ab_com + from->R*sin(from->theta*M_PI/180.0);
			c.z = from->R*cos(from->theta*M_PI/180.0);

			to->r_ab = from->r;
			to->r_bc = coor_xyz_distance(&b, &c);
			to->r_ac = coor_xyz_distance(&a, &c);
		}
		break;

		default:
			PRINT_ERROR("invalid arrangement %c\n", from->arrang)
			exit(EXIT_FAILURE);
	}
}
