/******************************************************************************

 About
 -----

 This module is an interface to an external user defined potential energy
 surface routine provided during compilation by the macro EXTERNAL_PES_NAME.


 References
 ----------

 [1] R. E. Olson et al. Phys. Rev. A. 3, 1607 (1971)
     doi: https://doi.org/10.1103/PhysRevA.3.1607

 [2] John C. Tully. J. Chem. Phys. 93, 1061 (1990)
     doi: http://dx.doi.org/10.1063/1.459170

******************************************************************************/

#include "cartesian.h"
#include "matrix.h"
#include "mass.h"
#include "pes.h"

/******************************************************************************

 Macro EXTERNAL_PES_NAME: is the name of an external user defined routine that
 evaluates the PES as function of a vector of coordinates, and defined during
 compilation. If not defined, a place holder called pes_noname is created for
 compatibility reasons.

 NOTE: calling an undefined EXTERNAL_PES_NAME() function, i.e. pes_noname(),
 is considered a runtime error not compilation error. 

******************************************************************************/

#if !defined(EXTERNAL_PES_NAME)
	#define EXTERNAL_PES_NAME pes_noname

	inline static double EXTERNAL_PES_NAME(const double r[])
	{
		PRINT_ERROR("%s\n", "no external PES defined")
		return r[0];
	}
#endif

extern double EXTERNAL_PES_NAME(const double r[]);

/******************************************************************************

 Function jacobi_to_internuc(): converts from a set of Jacobi coordinates to a
 set of internuclear distances for triatomic systems. The following order of
 arrangements is used: arrangement index 'a' for A + BC, 'b' for B + CA and 'c'
 is C + AB. Where, A, B and C represents three atoms.

******************************************************************************/

void jacobi_to_internuc(const jacobi_coor *from, internuc_coor *to)
{
	cartesian a, b, c;

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
			to->r_ac = cartesian_distance(&a, &c);
			to->r_ab = cartesian_distance(&a, &b);
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
			to->r_bc = cartesian_distance(&b, &c);
			to->r_ab = cartesian_distance(&a, &b);
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
			to->r_bc = cartesian_distance(&b, &c);
			to->r_ac = cartesian_distance(&a, &c);
		}
		break;

		default:
			PRINT_ERROR("invalid arrangement %c\n", from->arrang)
			exit(EXIT_FAILURE);
	}
}

/******************************************************************************

 Function jacobi_to_cartesian(): converts from a set of Jacobi coordinates to a
 set of internuclear distances for triatomic systems. The following order of
 arrangements is used: arrangement index 'a' for A + BC, 'b' for B + CA and 'c'
 is C + AB. Where, A, B and C represents three atoms.

******************************************************************************/

void jacobi_to_cartesian(const jacobi_coor *from,
                         cartesian *a, cartesian *b, cartesian *c)
{
	switch (from->arrang)
	{
		case 'a':
		{
			c->x = 0.0;
			c->y = from->r/2.0;
			c->z = 0.0;

			b->x =  0.0;
			b->y = -c->y;
			b->z =  0.0;

			const double mass_c = mass(atom_c);
			const double mass_b = mass(atom_b);
			const double cb_com = (c->y*mass_c + b->y*mass_b)/(mass_c + mass_b);

			a->x = 0.0;
			a->y = cb_com + from->R*sin(from->theta*M_PI/180.0);
			a->z = from->R*cos(from->theta*M_PI/180.0);
		}
		break;

		case 'b':
		{
			c->x = 0.0;
			c->y = from->r/2.0;
			c->z = 0.0;

			a->x =  0.0;
			a->y = -c->y;
			a->z =  0.0;

			const double mass_c = mass(atom_c);
			const double mass_a = mass(atom_a);
			const double ca_com = (c->y*mass_c + a->y*mass_a)/(mass_c + mass_a);

			b->x = 0.0;
			b->y = ca_com + from->R*sin(from->theta*M_PI/180.0);
			b->z = from->R*cos(from->theta*M_PI/180.0);
		}
		break;

		case 'c':
		{
			a->x = 0.0;
			a->y = from->r/2.0;
			a->z = 0.0;

			b->x =  0.0;
			b->y = -a->y;
			b->z =  0.0;

			const double mass_a = mass(atom_a);
			const double mass_b = mass(atom_b);
			const double ab_com = (a->y*mass_a + b->y*mass_b)/(mass_a + mass_b);

			c->x = 0.0;
			c->y = ab_com + from->R*sin(from->theta*M_PI/180.0);
			c->z = from->R*cos(from->theta*M_PI/180.0);
		}
		break;

		default:
			PRINT_ERROR("invalid arrangement %c\n", from->arrang)
			exit(EXIT_FAILURE);
	}
}

/******************************************************************************

 Wrapper pes(): is an interface for the external user defined PES as function
 of a set of Jacobi coordinates, x = (r, R, theta), which are translated to
 internuclear distances, y = (r_ab, r_bc, r_ac).

 NOTE: If the macro USE_NON_REACTIVE_PES is defined, the coordinate system is
 not converted and Jacobi is used all along, assuming

 EXTERNAL_PES_NAME(r, R, theta)

 Thus, make sure the order of each input parameter of EXTERNAL_PES_NAME follows
 the interface given above.

******************************************************************************/

double pes(const jacobi_coor *x)
{
	#if defined(USE_NON_REACTIVE_PES)
		double r[3];
		r[0] = x->r;
		r[1] = x->R;
		r[2] = x->theta;
	#elif defined(USE_CARTESIAN_COORDINATES)
		cartesian a, b, c;
		jacobi_to_cartesian(x, &a, &b, &c);

		double r[3*MAX_ATOM];
		r[0] = a.x;
		r[1] = a.y;
		r[2] = a.z;
		r[3] = b.x;
		r[4] = b.y;
		r[5] = b.z;
		r[6] = c.x;
		r[7] = c.y;
		r[8] = c.z;
	#else
		internuc_coor y;
		jacobi_to_internuc(x, &y);

		double r[MAX_INTERNUC_DISTANCE];
		r[0] = y.r_ab;
		r[1] = y.r_bc;
		r[2] = y.r_ac;
	#endif

	return EXTERNAL_PES_NAME(r);
}

/******************************************************************************

 Wrapper pes_sf(): is an interface for the user defined PES as function of a
 set of n Jacobi vectors in a space-fixed (SF) reference frame. Where, n =
 set->size. Jacobi coordinates are transformed in the user's internuclear
 distances.

******************************************************************************/
/*
double pes_sf(const jacobi_sf *set)
{
	#if defined(USE_CARTESIAN_COORDINATES)
	{
		cartesian com[MAX_JACOBI_VECTOR], atom[MAX_JACOBI_VECTOR + 1];

		jacobi_to_cartesian(set, com, atom);
		return _EXTERNAL_PES_NAME(atom);
	}
	#else
	{
		double r[MAX_INTERNUC_DISTANCE];

		jacobi_to_internuc(set, r);
		return _EXTERNAL_PES_NAME(r);
	}
	#endif
}*/

/******************************************************************************
******************************************************************************/

void pes_init()
{
	const jacobi_coor x =
	{
		.r = 1000.0,
		.R = 1000.0,
		.theta = 90.0,
		.arrang = 'a'
	};

	const bool is_nan = isnan(pes(&x));

	ASSERT(is_nan == false)
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
		.R = 1000.0,
		.theta = 90.0,
		.arrang = arrang
	};

	return pes(&x);
}

/******************************************************************************

 Function pes_asymptotic_min(): returns the min value of a diatomic potential,
 V(r) as function of the internuclear distance, r, asymptotically as R -> inf,
 for a given PES arrangement.

 NOTE: the PES is scanned in the interval r = [0, 50], R = inf and theta = 90.

******************************************************************************/

double pes_asymptotic_min(const char arrang, const double scan_step)
{
	ASSERT(scan_step > 0.0)

	double min = INF;
	const int grid_size = as_int(50.0/scan_step);

	for (int n = 0; n < grid_size; ++n)
	{
		const double r = as_double(n)*scan_step;
		const double energy = pec(arrang, r);

		if (energy < min) min = energy;
	}

	return min;
}

/******************************************************************************

 Function pes_olson_smith_model(): return the 2-by-2 model potential of Olson
 and Smith for the system Ne + He^+ at a given internuclear distance x. See
 Eq. (45) and Table I from Ref. [1] for details.

 NOTE: a handy benchmark for the algorithms.

******************************************************************************/

matrix *pes_olson_smith_model(const double x)
{
	matrix *m = matrix_alloc(2, 2, false);

	matrix_set(m, 0, 0, 21.1*exp(-x/0.678)/x);
	matrix_set(m, 0, 1, 0.170*exp(-x/0.667));
	matrix_set(m, 1, 0, 0.170*exp(-x/0.667));
	matrix_set(m, 1, 1, (21.1/x - 12.1)*exp(-x/0.678) + 16.8/27.2113961);

	return m;
}

/******************************************************************************

 Function pes_tully_model(): return one of the three (n = 1, 2, 3) 2-by-2 model
 potentials of J. C. Tully at a given x. See Eq. (21), (23) and (24) of
 Ref. [2].

******************************************************************************/

matrix *pes_tully_model(const int n, const double x)
{
	matrix *m = matrix_alloc(2, 2, false);

	switch (n)
	{
		case 1:
			if (x > 0.0)
			{
				matrix_set(m, 0, 0,  0.01*(1.0 - exp(-1.6*x)));
				matrix_set(m, 1, 1, -0.01*(1.0 - exp(-1.6*x)));
			}
			else
			{
				matrix_set(m, 0, 0, -0.01*(1.0 - exp(1.6*x)));
				matrix_set(m, 1, 1,  0.01*(1.0 - exp(1.6*x)));
			}

			matrix_set(m, 0, 1, 0.005*exp(-1.0*x*x));
			matrix_set(m, 1, 0, 0.005*exp(-1.0*x*x));
		break;

		case 2:
			matrix_set(m, 0, 0,  0.0);
			matrix_set(m, 0, 1,  0.015*exp(-0.06*x*x));
			matrix_set(m, 1, 0,  0.015*exp(-0.06*x*x));
			matrix_set(m, 1, 1, -0.10*exp(-0.28*x*x) + 0.05);
		break;

		case 3:
			matrix_set(m, 0, 0,  0.0006);
			matrix_set(m, 1, 1, -0.0006);

			matrix_set(m, 1, 0, 0.0);

			if (x < 0.0)
				matrix_set(m, 0, 1, 0.10*exp(0.90*x));

			else
				matrix_set(m, 0, 1, 0.10*(2.0 - exp(-0.90*x)));
		break;

		default:
			PRINT_ERROR("invalid index n = %d\n", n)
			exit(EXIT_FAILURE);
	}

	return m;
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
