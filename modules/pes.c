/******************************************************************************

 About
 -----

 This module is an interface to an external user defined function that returns
 a potential energy surface (PES), provided during compilation by the macro
 EXTERNAL_PES_NAME, i.e. EXTERNAL_PES_NAME(x). Where, x is a vector with
 the coordinates of the system. Internuclear distances are assume as a default
 coordinate system. The module also contains all routines needed to read, save,
 and output the respective atomic masses of the problem.


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

static double mass_a = 0.0, mass_b = 0.0, mass_c = 0.0, mass_d = 0.0, inf = 1000.0;

/******************************************************************************

 Macro MULTIPOLE_FILE_FORMAT:

******************************************************************************/

#if !defined(MULTIPOLE_FILE_FORMAT)
	#define MULTIPOLE_FILE_FORMAT "multipole_arrang=%c_n=%d.%s"
#endif

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

 Function pes_set_inf(): to define an asymptotic limit used internally. Default
 value is 1000 Bohrs.

******************************************************************************/

void pes_set_inf(const double x_inf)
{
	inf = x_inf;
}

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

	char *line = file_find(input, key);
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
 atoms A, B, C, etc. using the following format:

 atom_a = [symbol]
 atom_b = [symbol]
 atom_c = [symbol]
 (...)

 Where, [symbol] is expected as "1H", "2H", "4He", ..., "12C", etc. The routine
 will store the actual masses for a later use by pes_mass(), pes_mass_bc() etc.

 NOTE: Lines starting by '#' are ignored.

 NOTE: On entry, atom is given in lowercase, i.e. atom = 'a', 'b' or 'c'.

******************************************************************************/

void pes_init_mass(FILE *input, const char atom)
{
	ASSERT(input != NULL)

	switch (atom)
	{
		case 'a':
			mass_a = pes_read_mass(input, "atom_a");
			ASSERT(mass_a != 0.0)
			return;

		case 'b':
			mass_b = pes_read_mass(input, "atom_b");
			ASSERT(mass_b != 0.0)
			return;

		case 'c':
			mass_c = pes_read_mass(input, "atom_c");
			ASSERT(mass_c != 0.0)
			return;

		case 'd':
			mass_d = pes_read_mass(input, "atom_d");
			ASSERT(mass_d != 0.0)
			return;

		default:
			PRINT_ERROR("invalid atom %c\n", atom)
			exit(EXIT_FAILURE);
	}
}

/******************************************************************************

 Function pes_init(): a dummy call to the user defined PES such that it can
 initialize its own internal data, if any. Should be called before any other
 function from this module but after pes_init_mass().

******************************************************************************/

void pes_init()
{
	ASSERT(mass_a != 0.0)
	ASSERT(mass_b != 0.0)
	ASSERT(mass_c != 0.0)

	bool is_nan = true;

	if (mass_d == 0.0)
		is_nan = isnan(pes_abc('a', inf, inf, 90.0));
	else
		is_nan = isnan(pes_abcd(inf, inf, inf, 90.0, 90.0, 0.0));

	ASSERT(is_nan == false)

	math_no_gsl_handler();
	math_set_error(1.0e-2);
}

/******************************************************************************

 Function pes_mass(): return the respective masses of each atom a, b, c, etc.,
 previously initialized by pes_init_mass().

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

 Function pes_mass_bc(): return the reduced mass of BC.

******************************************************************************/

double pes_mass_bc()
{
	return mass_b*mass_c/(mass_b + mass_c);
}

/******************************************************************************

 Function pes_mass_ac(): return the reduced mass of AC.

******************************************************************************/

double pes_mass_ac()
{
	return mass_a*mass_c/(mass_a + mass_c);
}

/******************************************************************************

 Function pes_mass_ab(): return the reduced mass of AB.

******************************************************************************/

double pes_mass_ab()
{
	return mass_a*mass_b/(mass_a + mass_b);
}

/******************************************************************************

 Function pes_mass_abc(): return the reduced mass of ABC for a given
 arrangement.

******************************************************************************/

double pes_mass_abc(const char arrang)
{
	switch (arrang)
	{
		case 'a': return mass_a*(mass_b + mass_c)/(mass_a + mass_b + mass_c);
		case 'b': return mass_b*(mass_a + mass_c)/(mass_a + mass_b + mass_c);
		case 'c': return mass_c*(mass_a + mass_b)/(mass_a + mass_b + mass_c);

		default:
			PRINT_ERROR("invalid arrangement %c\n", arrang)
			exit(EXIT_FAILURE);
	}
}

/******************************************************************************

 Function pes_mass_abc(): return the reduced mass of the A+BCD arrangement.

******************************************************************************/

double pes_mass_abcd()
{
	return mass_a*(mass_b + mass_c + mass_d)/(mass_a + mass_b + mass_c + mass_d);
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

			const double cb = (c.y*mass_c + b.y*mass_b)/(mass_c + mass_b);

			a.x = 0.0;
			a.y = cb + R*sin(theta*M_PI/180.0);
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

			const double ca = (c.y*mass_c + a.y*mass_a)/(mass_c + mass_a);

			b.x = 0.0;
			b.y = ca + R*sin(theta*M_PI/180.0);
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

			const double ab = (a.y*mass_a + b.y*mass_b)/(mass_a + mass_b);

			c.x = 0.0;
			c.y = ab + R*sin(theta*M_PI/180.0);
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
	#endif

	return EXTERNAL_PES_NAME(internuc);
}

/******************************************************************************

 Wrapper pes_abcd(): an interface for the external user defined tetratomic PES
 as function of a set of Jacobi coordinates (r1, r2, theta12, R, theta, phi)
 for arrangement 'a' only (A + BC-D). Coordinates are translated to
 internuclear distances (ab, ac, ad, bc, bd, cd).

 NOTE: If the macro USE_NON_REACTIVE_PES is defined, the coordinate system is
 not converted and Jacobi is used all along, assuming

 EXTERNAL_PES_NAME(x) and x[0] = r1, x[1] = r2, x[2] = theta23, x[3] = R, x[4]
 = theta, x[5] = phi.

 Thus, make sure the order of each input parameter of EXTERNAL_PES_NAME follows
 the interface given above.

******************************************************************************/

double pes_abcd(const double r_bc,
                const double r_bcd,
                const double r_abcd,
                const double theta_bc,
                const double theta_a,
                const double phi_a)
{
	#if defined(USE_JACOBI_COORDINATES)
	{
		const double jacobi[6] = {r_bc, r_bcd, theta_bc, r_abcd, theta_a, phi_a};
		return EXTERNAL_PES_NAME(jacobi);
	}
	#endif

	/* NOTE: ab = 0, ac = 1, ad = 2, bc = 3, bd = 4, cd = 5. */
	double internuc[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

	cartesian a, b, c, d, cb, cbd;

	c.x = 0.0;
	c.y = r_bc/2.0;
	c.z = 0.0;

	b.x =  0.0;
	b.y = -c.y;
	b.z =  0.0;

	cb.x = 0.0;
	cb.y = (c.y*mass_c + b.y*mass_b)/(mass_c + mass_b);
	cb.z = 0.0;

	d.x = cb.x;
	d.y = cb.y + r_bcd*sin(theta_bc*M_PI/180.0);
	d.z = cb.z + r_bcd*cos(theta_bc*M_PI/180.0);

	const double mass_cb = mass_c + mass_b;

	cbd.x = (d.x*mass_d + cb.x*mass_cb)/(mass_d + mass_cb);
	cbd.y = (d.y*mass_d + cb.y*mass_cb)/(mass_d + mass_cb);
	cbd.z = (d.z*mass_d + cb.z*mass_cb)/(mass_d + mass_cb);

	a.x = cbd.x + r_abcd*sin(theta_a*M_PI/180.0)*cos(phi_a*M_PI/180.0);
	a.y = cbd.y + r_abcd*sin(theta_a*M_PI/180.0)*sin(phi_a*M_PI/180.0);
	a.z = cbd.z + r_abcd*cos(theta_a*M_PI/180.0);

	#if defined(USE_CARTESIAN_COORDINATES)
	{
		const double xyz[12]
			= {a.x, a.y, a.z, b.x, b.y, b.z, c.x, c.y, c.z, d.x, d.y, d.z};

		return EXTERNAL_PES_NAME(xyz);
	}
	#endif

	internuc[0] = cartesian_distance(&a, &b);
	internuc[1] = cartesian_distance(&a, &c);
	internuc[2] = cartesian_distance(&a, &d);
	internuc[3] = cartesian_distance(&b, &c);
	internuc[4] = cartesian_distance(&b, &d);
	internuc[5] = cartesian_distance(&c, &d);

	return EXTERNAL_PES_NAME(internuc);
}

/******************************************************************************

 Function pes_bc(): returns the asymptotic external user defined triatomic PES
 for arrangement 'a' and angular momentum j, i.e. the BC diatomic interaction.

******************************************************************************/

double pes_bc(const int j, const double r)
{
	const double mass = pes_mass_bc();

	return pes_abc('a', r, inf, 90.0) + as_double(j*(j + 1))/(2.0*mass*r*r);
}

/******************************************************************************

 Function pes_ac(): returns the asymptotic external user defined triatomic PES
 for arrangement 'b' and angular momentum j, i.e. the AC diatomic interaction.

******************************************************************************/

double pes_ac(const int j, const double r)
{
	const double mass = pes_mass_ac();

	return pes_abc('b', r, inf, 90.0) + as_double(j*(j + 1))/(2.0*mass*r*r);
}

/******************************************************************************

 Function pes_ab(): returns the asymptotic external user defined triatomic PES
 for arrangement 'c' and angular momentum j, i.e. the AB diatomic interaction.

******************************************************************************/

double pes_ab(const int j, const double r)
{
	const double mass = pes_mass_ab();

	return pes_abc('c', r, inf, 90.0) + as_double(j*(j + 1))/(2.0*mass*r*r);
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

static double pes_legendre_integrand(const double theta, const void *params)
{
	struct legendre_params *p = (struct legendre_params *) params;

	const double v = pes_abc(p->arrang, p->r, p->R, theta*180.0/M_PI);

	const double v_inf = pes_abc(p->arrang, p->r, inf, theta*180.0/M_PI);

	/* Eq. (22) of Ref. [1], inner integrand */
	return (v - v_inf)*math_legendre_poly(p->lambda, cos(theta))*sin(theta);
}

/******************************************************************************

 Function pes_legendre_multipole(): return the inner integral (in theta) of the
 triatomic PES at a given (r, R) Jacobi coordinate, as shown in Eq. (22) of
 Ref. [1].

 NOTE: Integration over theta in [0, pi] performed by a Gauss-Legendre
 quadrature rule with 64 Gauss points.

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

	const double result
		= math_gauss_legendre(0.0, M_PI, 64, &p, pes_legendre_integrand);

	return as_double(2*lambda + 1)*result/2.0;
}

/******************************************************************************

 Function pes_legendre_multipole(): return the inner integral (in theta) of the
 triatomic PES at a given (r, R) Jacobi coordinate, as shown in Eq. (22) of
 Ref. [1].

 NOTE: Integration over theta in [0, pi] performed by the QAG algorithm.

******************************************************************************/

/*
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
}*/

/******************************************************************************

 Function legendre_integrand(): is an auxiliary routine which returns the atom-
 diatom Legendre multipole integrand for the inner integral in Eq. (22) of Ref.
 [1] at a given (r, R) Jacobi coordinate.

******************************************************************************/

struct harmonics_params
{
	double phi_a;
	const int eta, m_eta;
	const double r_bc, r_bcd, r_abcd, theta_bc;
};

static double pes_theta_integrand(const double theta, void *params)
{
	struct harmonics_params *p = (struct harmonics_params *) params;

	const double theta_a = theta*180.0/M_PI;

	const double v
		= pes_abcd(p->r_bc, p->r_bcd, p->r_abcd, p->theta_bc, theta_a, p->phi_a);

	const double v_inf
		= pes_abcd(p->r_bc, p->r_bcd, inf, p->theta_bc, theta_a, p->phi_a);

	const double y
		= math_sphe_harmonics(p->eta, p->m_eta, theta_a, p->phi_a);

	return (v - v_inf)*y*sin(theta);
}

/******************************************************************************

 Function legendre_integrand(): is an auxiliary routine which returns the atom-
 diatom Legendre multipole integrand for the inner integral in Eq. (22) of Ref.
 [1] at a given (r, R) Jacobi coordinate.

******************************************************************************/

static double pes_phi_integrand(const double phi, void *params)
{
	struct harmonics_params *p = (struct harmonics_params *) params;

	p->phi_a = phi*180.0/M_PI;

	return math_qags(0.0, M_PI, &p, pes_theta_integrand);
}

/******************************************************************************

 Function pes_legendre_multipole(): return the inner integral (in theta) of the
 triatomic PES at a given (r, R) Jacobi coordinate, as shown in Eq. (22) of
 Ref. [1].

 NOTE: Integration over theta in [0, pi] performed by the QAG algorithm.

******************************************************************************/

double pes_harmonics_multipole(const int eta,
                               const int m_eta,
                               const double r_bc,
                               const double r_bcd,
                               const double r_abcd,
                               const double theta_bc)
{
	struct harmonics_params p =
	{
		.r_bc = r_bc,
		.r_bcd = r_bcd,
		.r_abcd = r_abcd,
		.theta_bc = theta_bc,

		.eta = eta,
		.m_eta = abs(m_eta)
	};

	const double phase = (m_eta < 0? pow(-1.0, m_eta) : 1.0);

	return phase*math_qags(0.0, 2.0*M_PI, &p, pes_phi_integrand);
}

/******************************************************************************

 Function pes_multipole_file(): opens the multipole file for the n-th grid
 index, arrangement and total angular momentum, J. Where, mode is the file
 access mode of fopen() from the C library.

 NOTE: mode = "wb" for write + binary format and "rb" for read + binary format,
 extension used is .bin otherwise .dat is used assuming text mode ("w" or "r").

******************************************************************************/

FILE *pes_multipole_file(const char arrang,
                         const int n, const char mode[], const bool verbose)
{
	char filename[MAX_LINE_LENGTH], ext[4];

	if (strlen(mode) > 1 && mode[1] == 'b')
		sprintf(ext, "%s", "bin");
	else
		sprintf(ext, "%s", "dat");

	sprintf(filename, MULTIPOLE_FILE_FORMAT, arrang, n, ext);

	FILE *stream = fopen(filename, mode);

	if (verbose)
	{
		if (stream != NULL && mode[0] == 'r') printf("# Reading %s\n", filename);
		if (stream != NULL && mode[0] == 'w') printf("# Writing %s\n", filename);
	}

	return stream;
}

/******************************************************************************

 Function pes_multipole_init(): allocates a set of multipole coefficients for a
 given number of grid points and lambda terms. Where, on entry m->value is NULL.

******************************************************************************/

void pes_multipole_init(pes_multipole *m)
{
	ASSERT(m->value == NULL)
	ASSERT(m->grid_size > 0)
	ASSERT(m->lambda_min > -1)
	ASSERT(m->lambda_step > 0)
	ASSERT(m->lambda_max >= m->lambda_min)

	m->value = allocate(m->lambda_max + 1, sizeof(double *), true);

	for (int lambda = m->lambda_min; lambda <= m->lambda_max; lambda += m->lambda_step)
	{
		m->value[lambda] = allocate(m->grid_size, sizeof(double), false);
	}
}

/******************************************************************************

 Function pes_multipole_init_all(): the same as pes_multipole_init() but for an
 array of n sets.

******************************************************************************/

void pes_multipole_init_all(const int n_max, pes_multipole m[])
{
	ASSERT(n_max > 0)

	for (int n = 0; n < n_max; ++n)
	{
		pes_multipole_init(&m[n]);
	}
}

/******************************************************************************

 Function pes_multipole_write(): writes in the disk a set of multipole terms in
 binary format.

******************************************************************************/

void pes_multipole_write(const pes_multipole *m, FILE *output)
{
	ASSERT(m != NULL)
	ASSERT(m->value != NULL)

	file_write(&m->R, sizeof(double), 1, output);

	file_write(&m->r_min, sizeof(double), 1, output);

	file_write(&m->r_max, sizeof(double), 1, output);

	file_write(&m->r_step, sizeof(double), 1, output);

	file_write(&m->lambda_min, sizeof(int), 1, output);

	file_write(&m->lambda_max, sizeof(int), 1, output);

	file_write(&m->lambda_step, sizeof(int), 1, output);

	file_write(&m->grid_size, sizeof(int), 1, output);

	for (int lambda = m->lambda_min; lambda <= m->lambda_max; lambda += m->lambda_step)
		file_write(m->value[lambda], sizeof(double), m->grid_size, output);
}

/******************************************************************************

 Function pes_multipole_write_all(): the same as pes_multipole_write() but for
 an array of n sets. 

******************************************************************************/

void pes_multipole_write_all(const int n_max,
                             const pes_multipole m[], FILE *output)
{
	ASSERT(n_max > 0)

	for (int n = 0; n < n_max; ++n)
	{
		pes_multipole_write(&m[n], output);
	}
}

/******************************************************************************

 Function pes_multipole_read(): reads from the disk a set of multipole terms in
 binary format.

******************************************************************************/

void pes_multipole_read(pes_multipole *m, FILE *input)
{
	file_read(&m->R, sizeof(double), 1, input, 0);

	file_read(&m->r_min, sizeof(double), 1, input, 0);

	file_read(&m->r_max, sizeof(double), 1, input, 0);

	file_read(&m->r_step, sizeof(double), 1, input, 0);

	file_read(&m->lambda_min, sizeof(int), 1, input, 0);

	file_read(&m->lambda_max, sizeof(int), 1, input, 0);

	file_read(&m->lambda_step, sizeof(int), 1, input, 0);

	file_read(&m->grid_size, sizeof(int), 1, input, 0);

	pes_multipole_init(m);

	for (int lambda = m->lambda_min; lambda <= m->lambda_max; lambda += m->lambda_step)
	{
		file_read(m->value[lambda], sizeof(double), m->grid_size, input, 0);
	}
}

/******************************************************************************

 Function pes_multipole_read_all(): the same as pes_multipole_read() but for
 an array of n sets.

******************************************************************************/

void pes_multipole_read_all(const int n_max, pes_multipole m[], FILE *input)
{
	ASSERT(n_max > 0)

	for (int n = 0; n < n_max; ++n)
	{
		pes_multipole_read(&m[n], input);
	}
}

/******************************************************************************

 Function pes_multipole_free():

******************************************************************************/

void pes_multipole_free(pes_multipole *m)
{
	ASSERT(m != NULL)

	if (m->value == NULL) return;

	for (int lambda = m->lambda_min; lambda <= m->lambda_max; lambda += m->lambda_step)
		if (m->value[lambda] != NULL) free(m->value[lambda]);

	free(m->value);
	m->value = NULL;
}

/******************************************************************************

 Function pes_about(): prints in a given output file the conditions in which
 the module was compiled.

******************************************************************************/

void pes_about(FILE *output)
{
	ASSERT(output != NULL)

	fprintf(output, "# build date   = %s\n", __DATE__);
	fprintf(output, "# source code  = %s\n", __FILE__);
	fprintf(output, "# external PES = %s\n", PRINT_MACRO(EXTERNAL_PES_NAME));

	#if defined(USE_JACOBI_COORDINATES)
		fprintf(output, "# coordinates  = Jacobi\n");
	#endif

	#if defined(USE_CARTESIAN_COORDINATES)
		fprintf(output, "# coordinates  = Cartesian\n");
	#endif

	#if !defined(USE_JACOBI_COORDINATES) && !defined(USE_CARTESIAN_COORDINATES)
		fprintf(output, "# coordinates  = internuclear distances\n");
	#endif
}
