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
#include "nist.h"
#include "math.h"
#include "file.h"
#include "pes.h"

static double mass_a = 0.0, mass_b = 0.0, mass_c = 0.0, mass_d = 0.0;

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

 Function pes_read_mass(): an auxiliary routine that scan a given input file
 searching for the first occurrence of a keyword of format "[key] = [value]".
 Where, [value] is expected as an atomic symbol, e.g. "1H", "2H", "4He", ...,
 "12C", etc. The actual atomic mass is returned if the keyword is found.

 NOTE: Lines starting by '#' are ignored.

******************************************************************************/

static double pes_read_mass(FILE *input, const char key[])
{
	ASSERT(input != NULL)

	char *line = file_find_string(input, key);
	char *token = strtok(line, "=");

	while (token != NULL)
	{
		if (strstr(line, key) != NULL)
		{
			token = strtok(NULL, "=");
			token = trim(token);

			const isotope s = nist_isotope(token);

			/* NOTE: masses are returned in atomic units. */
			return nist_atomic_mass(s)*1822.8873843;
		}

		token = strtok(NULL, "=");
	}

	PRINT_ERROR("no entry '%s' found\n", key)
	exit(EXIT_FAILURE);
}

/******************************************************************************

 Function pes_init_mass(): reads from a given input the atomic symbols for
 atoms A, B, C, etc using the following format:

 atom_a = [symbol]
 atom_b = [symbol]
 atom_c = [symbol]
 .
 .
 .

 Where, [symbol] is expected as "1H", "2H", "4He", ..., "12C", etc.

 The routine will store the actual masses for a later use by pes_mass().

 NOTE: Lines starting by '#' are ignored.

******************************************************************************/

void pes_init_mass(FILE *input, const char atom)
{
	ASSERT(input != NULL)

	switch (atom)
	{
		case 'a': return mass_a = read_atomic_mass(input, "atom_a");
		case 'b': return mass_b = read_atomic_mass(input, "atom_b");
		case 'c': return mass_c = read_atomic_mass(input, "atom_c");
		case 'd': return mass_d = read_atomic_mass(input, "atom_d");

		default:
			PRINT_ERROR("invalid atom %c\n", atom)
			exit(EXIT_FAILURE);
	}
}

/******************************************************************************

 Function pes_mass(): return the respective masses (in atomic units) of each
 atom a, b, c, ... etc previously initialized by pes_init_mass().

******************************************************************************/

double pes_mass(const char atom)
{
	switch (atom)
	{
		case 'a': return mass_a;
		case 'b': return mass_b;
		case 'c': return mass_c;
		case 'd': return mass_d;

		default:
			PRINT_ERROR("invalid atom %c\n", atom)
			exit(EXIT_FAILURE);
	}
}

/******************************************************************************

 Wrapper pes_abc(): an interface for the external user defined triatomic PES as
 function of a set of Jacobi coordinates (r, R, theta) for a given arrangement,
 which are translated to internuclear distances (ab, bc, ac).

 NOTE: If the macro USE_NON_REACTIVE_PES is defined, the coordinate system is
 not converted and Jacobi is used all along, assuming

 EXTERNAL_PES_NAME(x) and x[0] = r, x[1] = R, x[2] = theta

 Thus, make sure the order of each input parameter of EXTERNAL_PES_NAME follows
 the interface given above.

******************************************************************************/

double pes_abc(const char arrang,
               const double r, const double R, const double theta)
{
	#if defined(USE_JACOBI_COORDINATES)
	{
		const double jacobi[3] = {r, R, theta};		
		return EXTERNAL_PES_NAME(jacobi);
	}
	#endif

	/* NOTE: bc = 0, ac = 1, ab = 2. */
	double internuc[3] = {0.0, 0.0, 0.0};

	cartesian a, b, c;

	switch (arrang)
	{
		case 'a':
		{
			c.x = 0.0;
			c.y = r/2.0;
			c.z = 0.0;

			b.x =  0.0;
			b.y = -c.y;
			b.z =  0.0;

			const double cb_com = (c.y*mass_c + b.y*mass_b)/(mass_c + mass_b);

			a.x = 0.0;
			a.y = cb_com + R*sin(theta*M_PI/180.0);
			a.z = R*cos(theta*M_PI/180.0);

			internuc[0] = r;
			internuc[1] = cartesian_distance(&a, &c);
			internuc[2] = cartesian_distance(&a, &b);
		}
		break;

		case 'b':
		{
			c.x = 0.0;
			c.y = r/2.0;
			c.z = 0.0;

			a.x =  0.0;
			a.y = -c.y;
			a.z =  0.0;

			const double ca_com = (c.y*mass_c + a.y*mass_a)/(mass_c + mass_a);

			b.x = 0.0;
			b.y = ca_com + R*sin(theta*M_PI/180.0);
			b.z = R*cos(theta*M_PI/180.0);

			internuc[1] = r;
			internuc[0] = cartesian_distance(&b, &c);
			internuc[2] = cartesian_distance(&a, &b);
		}
		break;

		case 'c':
		{
			a.x = 0.0;
			a.y = r/2.0;
			a.z = 0.0;

			b.x =  0.0;
			b.y = -a.y;
			b.z =  0.0;

			const double ab_com = (a.y*mass_a + b.y*mass_b)/(mass_a + mass_b);

			c.x = 0.0;
			c.y = ab_com + R*sin(theta*M_PI/180.0);
			c.z = R*cos(theta*M_PI/180.0);

			internuc[2] = r;
			internuc[0] = cartesian_distance(&b, &c);
			internuc[1] = cartesian_distance(&a, &c);
		}
		break;

		default:
			PRINT_ERROR("invalid arrangement %c\n", arrang)
			exit(EXIT_FAILURE);
	}

	#if defined(USE_CARTESIAN_COORDINATES)
	{
		const double xyz[9] = {a.x, a.y, a.z, b.x, b.y, b.z, c.x, c.y, c.z};		
		return EXTERNAL_PES_NAME(xyz);
	}

	return EXTERNAL_PES_NAME(internuc);
}

/******************************************************************************

 Function pes_bc(): returns the asymptotic external user defined triatomic PES
 for arrangement a and angular momentum j, i.e. the BC diatomic interaction.

******************************************************************************/

double pes_bc(const int j, const double r)
{
	const double mass_bc = mass_b*mass_c/(mass_b + mass_c);
	ASSERT(mass_bc != 0.0)

	return pes_abc('a', r, 1000.0, 90.0) + as_double(j*(j + 1))/(2.0*mass_bc*r*r);
}

/******************************************************************************

 Function pes_ac(): returns the asymptotic external user defined triatomic PES
 for arrangement b and angular momentum j, i.e. the AC diatomic interaction.

******************************************************************************/

double pes_ac(const int j, const double r)
{
	const double mass_ac = mass_a*mass_c/(mass_a + mass_c);
	ASSERT(mass_ac != 0.0)

	return pes_abc('b', r, 1000.0, 90.0) + as_double(j*(j + 1))/(2.0*mass_ac*r*r);
}

/******************************************************************************

 Function pes_ab(): returns the asymptotic external user defined triatomic PES
 for arrangement c and angular momentum j, i.e. the AB diatomic interaction.

******************************************************************************/

double pes_ab(const int j, const double r)
{
	const double mass_ab = mass_a*mass_b/(mass_a + mass_b);
	ASSERT(mass_ab != 0.0)

	return pes_abc('b', r, 1000.0, 90.0) + as_double(j*(j + 1))/(2.0*mass_ab*r*r);
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

 Function legendre_integrand(): is an auxiliary routine which returns the atom-
 diatom Legendre multipole integrand for the inner integral in Eq. (22) of Ref.
 [1] at a given (r, R) Jacobi coordinate.

******************************************************************************/

struct legendre_params
{
	const char arrang;
	const int lambda;
	const double r;
	const double R;
};

static double pes_legendre_integrand(const double theta, void *params)
{
	struct legendre_params *p = (struct legendre_params *) params;

	const double pot_energy = pes_abc(p->arrang, p->r, p->R,   theta*180.0/M_PI)
	                        - pes_abc(p->arrang, p->r, 1000.0, theta*180.0/M_PI);

	/* Eq. (22) of Ref. [1], inner integrand */
	return pot_energy*math_legendre_poly(p->lambda, cos(theta))*sin(theta);
}

/******************************************************************************

 Function pes_legendre_multipole(): return the inner integral (in theta) of the
 triatomic PES at a given (r, R) Jacobi coordinate, as shown in Eq. (22) of
 Ref. [1].

 NOTE: Integration over theta in [0, pi] performed by the QAG algorithm.

******************************************************************************/

double pes_legendre_multipole(const char arrang,
                              const int lambda, const double r, const double R)
{
	struct legendre_params p =
	{
		.arrang = arrang,
		.lambda = lambda,
		.r = r,
		.R = R
	};

	double factor = 1.0, theta_max = M_PI;

	if (arrang == 'a' && mass_b == mass_c)
	{
		 factor = 2.0;
		 theta_max = M_PI/2.0;
	}

	if (arrang == 'b' && mass_c == mass_a)
	{
		 factor = 2.0;
		 theta_max = M_PI/2.0;
	}

	if (arrang == 'c' && mass_a == mass_b)
	{
		 factor = 2.0;
		 theta_max = M_PI/2.0;
	}

	const double result = math_qags(0.0, theta_max, &p, pes_legendre_integrand);

	return as_double(2*lambda + 1)*factor*result/2.0;
}

/******************************************************************************

 Function legendre_integrand(): is an auxiliary routine which returns the atom-
 diatom Legendre multipole integrand for the inner integral in Eq. (22) of Ref.
 [1] at a given (r, R) Jacobi coordinate.

******************************************************************************/

struct _params
{
	double phi;
	const int l, m;
	const double r1, r2, r3, theta12;
};

static double pes_theta_integrand(const double theta, void *params)
{
	struct _params *p = (struct _params *) params;

	const jacobi_6d x =
	{
		.r1 = p->r1,
		.r2 = p->r2,
		.r3 = p->r3,
		.theta12 = p->theta12,
		.theta3 = theta*180.0/M_PI,
		.phi3 = p->phi*180.0/M_PI
	};

	const jacobi_6d x_inf =
	{
		.r1 = p->r1,
		.r2 = p->r2,
		.r3 = 1000.0,
		.theta12 = p->theta12,
		.theta3 = theta*180.0/M_PI,
		.phi3 = p->phi*180.0/M_PI
	};

	const double pot_energy = pes(&x) - pes(&x_inf);

	return pot_energy*math_sphe_harmonics(p->l, p->m, theta, p->phi)*sin(theta);
}

/******************************************************************************

 Function legendre_integrand(): is an auxiliary routine which returns the atom-
 diatom Legendre multipole integrand for the inner integral in Eq. (22) of Ref.
 [1] at a given (r, R) Jacobi coordinate.

******************************************************************************/

static double pes_phi_integrand(const double phi, void *params)
{
	struct _params *p = (struct _params *) params;

	p->phi = phi;

	return math_qags(0.0, M_PI, &p, pes_theta_integrand);
}

/******************************************************************************

 Function pes_legendre_multipole(): return the inner integral (in theta) of the
 triatomic PES at a given (r, R) Jacobi coordinate, as shown in Eq. (22) of
 Ref. [1].

 NOTE: Integration over theta in [0, pi] performed by the QAG algorithm.

******************************************************************************/

double pes_harmonic_multipole(const int l,
                              const int m,
                              const double r1,
                              const double r2,
                              const double r3,
                              const double theta)
{
	struct _params p =
	{
		.r1 = r1,
		.r2 = r2,
		.r3 = r3,
		.theta12 = theta,
		.l = l,
		.m = m
	};

	return math_qags(0.0, 2.0*M_PI, &p, pes_phi_integrand);
}

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
