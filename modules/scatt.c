/******************************************************************************

 AUTHOR
 ======

 Written by Humberto Jr, 2018-2019
 Last modified: Jul 02, 2019


 ABOUT
 =====

 This module is a collection of routines designed to solve problems related to
 atomic and molecular scattering. It includes general-purpose functions to man
 ipulate files, pointers, linear algebra, several types of data structures and
 physics.

 All written in C (C99 standard) and the main dependency is the GNU Scientific
 Library (GSL).

 Atomic units are used through unless explicitly stated otherwise.


 HOW TO BUILD
 ============

 Suppose a C compiler CC (e.g. CC=gcc or CC=icc) and a GSL installation in the
 path given by GSLROOT (e.g. GSLROOT=/usr/local),

 $CC -W -Wall -std=c99 -pedantic -fopenmp -O3 -I$GSLROOT/include -c filename.c

 is a fairly general method to compile the module, and it is relatively simple
 to adapt for a given choice of CC. Specifying the C standard and pedantic opt
 ion is always recommended in order to avoid the use of C++ rules and others C
 specifications.

 Some optional macros are available to tune the compilation:

 USE_MKL              = switch on the use of Intel Math Kernel Library.
 USE_ESSL             = switch on the use of IBM Engineering and Scientific
                        Subroutine Library.
 USE_LAPACKE          = switch on the use of LAPACKE library.
 USE_DUMMY_PES        = switch off the use of an user defined PES routine.
 MAX_LINE_LENGTH      = maximum length of all strings used in the module.
 GSL_MAX_WORKSPACE    = maximum number of elements used in all GSL workspace.
 USE_NON_REACTIVE_PES = uses Jacobi coordinates instead of internuclear dista
                        nces when invoking the user defined PES routine.

 Often, the compiler option used to define macros is -D, e.g.

 -DUSE_MKL -DUSE_NON_REACTIVE_PES

 Thus, check out the documentation of your choice of CC.


 NAMING CONVENTION
 =================

 pes        = potential energy surface.
 pec        = potential energy curve.
 grid_size  = number of points in a grid mesh.
 grid_step  = step of a given variable in a grid.
 input      = a FILE object open with mode = 'r' (read only).
 output     = a FILE object open with mode = 'w' (write only).
 max_row    = maximum number of rows (matrices, data files etc).
 max_col    = maximum number of columns (matrices, data files etc).
 max_ch     = maximum number of (scattering) channels.
 ch         = (scattering) channel.
 wavef      = wavefunction.
 arrang     = arrangement; 'a' for A+BC, 'b' for B+CA, 'c' for C+AB etc.
 eigenvec   = eigenvector.
 eigenval   = eigenvalue.
 spin_mult  = spin multiplicity; 1 for singlet, 2 for doublet etc.
 pot_energy = potential energy (sometimes named v).
 kin_energy = kinetic energy.
 tot_energy = total energy.
 use_omp    = switch on/off the use of OpenMP.
 wavenum    = wavenumber; sqrt(2.0*mass*E), where E = E_tot - V(inf).


 DEV NOTES
 =========

 1) Variables are set const whenever possible, for the sake of correctness.

 2) Due to (1), on exit, all non-const input parameters are necessarily
    modified and/or output.

 3) On entry, all input representing sizes and pointers are always checked.
    The use of size_t types is not recommended.

 4) Maximum length for strings is given by the macro MAX_LINE_LENGTH.

 5) Matrices are always ordinary vectors with elements stored in a row-major
    scheme: vector[n*max_col + m], where n = [0, max_row) and m = [0, max_col).

 6) Priority for efficiency. However complicated algorithms are sometimes
    written in a non-efficient approach, provided the resulting code is
    more readble and error-free.

 7) The use of OpenMP in this module (often, use_omp = true) implies it is not
    used anywhere else, including the external linear algebra libraries.

 8) This module is thread-safe intended.


 REFERENCES
 ==========

 [1] M. Hernandez Vera et al. J. Chem. Phys. 146, 124310 (2017)
     doi: https://doi.org/10.1063/1.4978475

 [2] R. T. Pack J. Chem. Phys. 60, 633 (1974)
     doi: https://doi.org/10.1063/1.1681085

 [3] O. Dulieu et al. J. Chem. Phys. 103 (1) (1995)
     doi: 10.1063/1.469622

 [4] B. R. Johnson J. Chem. Phys. 69, 4678 (1978)
     doi: https://doi.org/10.1063/1.436421

 [5] V. Kokoouline et al. J. Chem. Phys. 110, 20 (1999)
     doi:

 [6] W. H. Miller. J. Chem. Phys. 1, 50 (1969)
     doi:

 [7] R. B. Bernstein et al. Proc. R. Soc. Lond. A 1963 274 , 427-442
     doi: https://doi.org/10.1098/rspa.1963.0142

 [8] A. E. DePristo et al. J. Phys. B: Atom. Molec. Phys., Vol. 9, No. 3 (1976)
     doi:

 [9] D. D. Lopez-Duran et al. Com. Phys. Comm. 179 (2008) 821-838
     doi: https://doi.org/10.1016/j.cpc.2008.07.017

 [10] B. R. Johnson J. Comp. Phys. 13, 445-449 (1973)
      doi: https://doi.org/10.1016/0021-9991(73)90049-1

 [11] B. R. Johnson J. Chem. Phys. 67, 4086 (1977)
      doi: https://doi.org/10.1063/1.435384

 [12] M. Monnerville et al. J. Chem. Phys. 101, 7580 (1994)
      doi: http://dx.doi.org/10.1063/1.468252

 [13] R. E. Olson et al. Phys. Rev. A. 3, 1607 (1971)
      doi: https://doi.org/10.1103/PhysRevA.3.1607

 [14] B. R. Johnson. J. Comp. Phys. 13, 445 (1973)
      doi: https://doi.org/10.1016/0021-9991(73)90049-1

 [15] John C. Tully. J. Chem. Phys. 93, 1061 (1990)
      doi: http://dx.doi.org/10.1063/1.459170

******************************************************************************/

#include <omp.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_coupling.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_permutation.h>

#include "globals.h"
#include "matrix.h"
#include "amo.h"

#if !defined(GSL_MAX_WORKSPACE)
	#define GSL_MAX_WORKSPACE 50000
#endif

#if !defined(HAVE_INLINE)
	/* NOTE: to force GSL inline functions. */
	#define HAVE_INLINE
#endif

#if defined(USE_DUMMY_PES)
	inline static double pot_(const double *r_ab, const double *r_bc, const double *r_ac)
	{
		PRINT_ERROR("dummy PES invoked for r = (%f, %f, %f)\n", *r_ab, *r_bc, *r_ac)
		exit(EXIT_FAILURE);
	}
#else
	extern double pot_(const double *r_ab, const double *r_bc, const double *r_ac);
#endif

/******************************************************************************

 Function save_ch(): writes an instance, ch, of type scatt_basis in a binary
 file named filename.

 TODO: rename to save_scatt_basis. Used in umatrix.

******************************************************************************/

void save_ch(const char filename[], const scatt_basis *ch)
{
	FILE *output = fopen(filename, "wb");

	if (output == NULL) {
		PRINT_ERROR("unable to open %s\n", filename)
		exit(EXIT_FAILURE);
	}

	int status = 0;
	status = fwrite(&ch->v, sizeof(int), 1, output);
	status = fwrite(&ch->j, sizeof(int), 1, output);
	status = fwrite(&ch->l, sizeof(int), 1, output);
	status = fwrite(&ch->J, sizeof(int), 1, output);
	status = fwrite(&ch->spin_mult, sizeof(int), 1, output);
	status = fwrite(&ch->grid_size, sizeof(int), 1, output);

	status = fwrite(&ch->arrang, sizeof(char), 1, output);

	status = fwrite(&ch->n, sizeof(double), 1, output);
	status = fwrite(&ch->r_min, sizeof(double), 1, output);
	status = fwrite(&ch->r_max, sizeof(double), 1, output);
	status = fwrite(&ch->energy, sizeof(double), 1, output);

	if (ch->grid_size > 0 && ch->wavef != NULL) {
		status = fwrite(ch->wavef, sizeof(double), ch->grid_size, output);

		if (status != ch->grid_size) {
			PRINT_ERROR("only %d elements from a total of %d were written\n", status, ch->grid_size)
			exit(EXIT_FAILURE);
		}
	}

	fclose(output);
}

/******************************************************************************

 Function load_ch(): reads an instance, ch, of type scatt_basis in a binary
 file named filename.

 TODO: rename to load_scatt_basis. Used in umatrix.

******************************************************************************/

void load_ch(const char filename[], scatt_basis *ch)
{
	FILE *input = fopen(filename, "rb");

	if (input == NULL) {
		PRINT_ERROR("unable to open %s\n", filename)
		exit(EXIT_FAILURE);
	}

	int status = 0;
	status = fread(&ch->v, sizeof(int), 1, input);
	status = fread(&ch->j, sizeof(int), 1, input);
	status = fread(&ch->l, sizeof(int), 1, input);
	status = fread(&ch->J, sizeof(int), 1, input);
	status = fread(&ch->spin_mult, sizeof(int), 1, input);
	status = fread(&ch->grid_size, sizeof(int), 1, input);

	status = fread(&ch->arrang, sizeof(char), 1, input);

	status = fread(&ch->n, sizeof(double), 1, input);
	status = fread(&ch->r_min, sizeof(double), 1, input);
	status = fread(&ch->r_max, sizeof(double), 1, input);
	status = fread(&ch->energy, sizeof(double), 1, input);

	if (ch->grid_size > 0 && ch->wavef == NULL) {
		ch->wavef = malloc(sizeof(double)*ch->grid_size);

		if (ch->wavef == NULL) {
			PRINT_ERROR("unable to allocate %d elements\n", ch->grid_size)
			exit(EXIT_FAILURE);
		}

		status = fread(ch->wavef, sizeof(double), ch->grid_size, input);

		if (status != ch->grid_size) {
			PRINT_ERROR("only %d elements from a total of %d were read\n", status, ch->grid_size)
			exit(EXIT_FAILURE);
		}
	}

	fclose(input);
}

/*****************************************************************************

 Function save_scatt_matrix(): saves a s_matrix object in the disk and its
 associated energy, if any. All in binary format.

******************************************************************************/

void save_scatt_matrix(const char filename[],
                       const double energy, const s_matrix *s)
{
	ASSERT(filename != NULL)
	FILE *output = fopen(filename, "wb");

	if (output == NULL)
    {
		PRINT_ERROR("unable to open %s\n", filename)
		exit(EXIT_FAILURE);
	}

	ASSERT(s != NULL)
	ASSERT(s->max_ch > 0)
	ASSERT(s->re_part != NULL)
	ASSERT(s->im_part != NULL)

	int status = fwrite(&energy, sizeof(double), 1, output);
	ASSERT(status == 1)

	status = fwrite(&s->max_ch, sizeof(int), 1, output);
	ASSERT(status == 1)

	status = fwrite(s->re_part, sizeof(double), s->max_ch*s->max_ch, output);
	ASSERT(status == s->max_ch*s->max_ch)

	status = fwrite(s->im_part, sizeof(double), s->max_ch*s->max_ch, output);
	ASSERT(status == s->max_ch*s->max_ch)

	fclose(output);
}

void debug_stamp(FILE *output)
{
	ASSERT(output != NULL)

	fprintf(output, "#\n");
	fprintf(output, "# Build: %s @%s\n", __DATE__, __TIME__);
	fprintf(output, "# LINEAR_ALGEBRA_LIB = %s\n", LINEAR_ALGEBRA_LIB);
	fprintf(output, "# GSL_MAX_WORKSPACE  = %d\n", GSL_MAX_WORKSPACE);
	fprintf(output, "# MAX_LINE_LENGTH    = %d\n", MAX_LINE_LENGTH);
}

/******************************************************************************

 Function mass(): works as an interface between an input file (plain text given
 in the C stdin) and the remaining routines of this module. It returns the many
 kind of masses needed, in atomic units. The user can define masses for atoms
 A, B and C in an external file using the following keywords:

 mass_a = [value]
 mass_b = [value]
 mass_c = [value]

 The routine will store these values and return them on demand, as well as any
 derived quantities defined by the enumeration mass_case. If mass_a, mass_b or
 mass_c are given not in atomic units, the respective conversion factor can be
 provided by the following keyword (1.0 by default):

 mass_scale_factor = [value]

******************************************************************************/

double mass(const enum mass_case m)
{
	static bool first_call = true;

	static double mass_a, mass_b, mass_c,
	              mass_ab, mass_bc, mass_ac, mass_abc, mass_bca, mass_cab;

	if (first_call)
	{
		first_call = false;

		mass_a = get_key(stdin, "mass_a", 0.0, INFINITY, 0.0);
		mass_b = get_key(stdin, "mass_b", 0.0, INFINITY, 0.0);
		mass_c = get_key(stdin, "mass_c", 0.0, INFINITY, 0.0);

		ASSERT(mass_a > 0.0 && mass_b > 0.0 && mass_c > 0.0)

		const double factor
			= get_key(stdin, "mass_scale_factor", -INFINITY, INFINITY, 1.0);

		mass_a = factor*mass_a;
		mass_b = factor*mass_b;
		mass_c = factor*mass_c;

		mass_ab = mass_a*mass_b/(mass_a + mass_b);
		mass_bc = mass_b*mass_c/(mass_b + mass_c);
		mass_ac = mass_a*mass_c/(mass_a + mass_c);

		mass_abc = mass_a*(mass_b + mass_c)/(mass_a + mass_b + mass_c);
		mass_bca = mass_b*(mass_c + mass_a)/(mass_a + mass_b + mass_c);
		mass_cab = mass_c*(mass_a + mass_b)/(mass_a + mass_b + mass_c);
	}

	switch (m)
	{
		case atom_a:   return mass_a;
		case atom_b:   return mass_b;
		case atom_c:   return mass_c;
		case pair_ab:  return mass_ab;
		case pair_bc:  return mass_bc;
		case pair_ac:  return mass_ac;
		case arrang_a: return mass_abc;
		case arrang_b: return mass_bca;
		case arrang_c: return mass_cab;
		case total:    return (mass_a + mass_b + mass_c);

		default:
			PRINT_ERROR("%s", "invalid choice of enum mass_case\n")
			exit(EXIT_FAILURE);
	}
}

/******************************************************************************

 Function jacobi_to_internuc(): converts from a set of Jacobi coordinates to a
 set of internuclear distances for triatomic systems. The following convention
 of arrangement notation is used: arrangement index 'a' for A + BC, 'b' for
 B + CA and 'c' is C + AB. Where, A, B and C represents three atoms.

******************************************************************************/

void jacobi_to_internuc(const jacobi_coor *from, internuc_coor *to)
{
	xyz_coor a, b, c;

	switch (from->arrang)
	{
		case 'a':
			c.x = 0.0;
			c.y = from->r/2.0;
			c.z = 0.0;

			b.x =  0.0;
			b.y = -c.y;
			b.z =  0.0;

			a.x = 0.0;

			a.y = (c.y*mass(atom_c) + b.y*mass(atom_b))/(mass(atom_c) + mass(atom_b))
			    + from->R*sin(from->theta*M_PI/180.0);

			a.z = from->R*cos(from->theta*M_PI/180.0);

			to->r_bc = from->r;
			to->r_ac = distance(&a, &c);
			to->r_ab = distance(&a, &b);
		break;

		case 'b':
			c.x = 0.0;
			c.y = from->r/2.0;
			c.z = 0.0;

			a.x =  0.0;
			a.y = -c.y;
			a.z =  0.0;

			b.x = 0.0;

			b.y = (c.y*mass(atom_c) + a.y*mass(atom_a))/(mass(atom_c) + mass(atom_a))
			    + from->R*sin(from->theta*M_PI/180.0);

			b.z = from->R*cos(from->theta*M_PI/180.0);

			to->r_ac = from->r;
			to->r_bc = distance(&b, &c);
			to->r_ab = distance(&a, &b);
		break;

		case 'c':
			a.x = 0.0;
			a.y = from->r/2.0;
			a.z = 0.0;

			b.x =  0.0;
			b.y = -a.y;
			b.z =  0.0;

			c.x = 0.0;

			c.y = (a.y*mass(atom_a) + b.y*mass(atom_b))/(mass(atom_a) + mass(atom_b))
			    + from->R*sin(from->theta*M_PI/180.0);

			c.z = from->R*cos(from->theta*M_PI/180.0);

			to->r_ab = from->r;
			to->r_bc = distance(&b, &c);
			to->r_ac = distance(&a, &c);
		break;

		default:
			PRINT_ERROR("invalid arrangement %c\n", from->arrang)
			exit(EXIT_FAILURE);
	}
}

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
		jacobi_to_internuc(x, &y);
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
		.R = 100000.0,
		.theta = 90.0,
		.arrang = arrang
	};

	return pes(&x);
}

/******************************************************************************

 Function riccati_bessel(): returns the Riccati-Bessel function J(l, x) if type
 = 'j', or N(l, x) if type = 'n', as defined in Eq. (14) and (15) of Ref. [10].
 These are made upon spherical Bessel functions, j(l, x) and y(l, x).

******************************************************************************/

double riccati_bessel(const char type,
                      const int l, const double wavenum, const double x)
{
	switch (type)
	{
		case 'j': return  (wavenum*x)*gsl_sf_bessel_jl(l, wavenum*x)/sqrt(wavenum);
		case 'n': return -(wavenum*x)*gsl_sf_bessel_yl(l, wavenum*x)/sqrt(wavenum);

		default: return 0.0;
	}
}

/******************************************************************************

 Function modif_spher_bessel(): returns the modified spherical Bessel functions
 J(l, x) if type = 'j', or N(l, x) if type = 'n', as defined in Eq. (16) and
 (17) of Ref. [10]. These are made upon regular and irregular modified
 spherical Bessel functions, i(l, x) and k(l, x).

******************************************************************************/

double modif_spher_bessel(const char type,
                          const int l, const double wavenum, const double x)
{
	switch (type)
	{
		case 'j': return sqrt(wavenum*x)*gsl_sf_bessel_il_scaled(l, wavenum*x);
		case 'n': return sqrt(wavenum*x)*gsl_sf_bessel_kl_scaled(l, wavenum*x);

		default: return 0.0;
	}
}

/******************************************************************************

 Function norm(): normalize to unity a radial wavefunction, g = g(r); where,
 the integral of the amplitude square over r is performed by the 1/3-Simpson
 quadrature rule.

******************************************************************************/

void norm(const int grid_size, const double grid_step, double wavef[])
{
	ASSERT(grid_size > 0)
	ASSERT(wavef != NULL)

	register const int n_max
		= (grid_size%2 == 0? grid_size : grid_size - 1);

	register double sum
		= wavef[0]*wavef[0] + wavef[n_max - 1]*wavef[n_max - 1];

	for (int n = 1; n < (n_max - 2); n += 2)
	{
		sum += 4.0*wavef[n]*wavef[n] + 2.0*wavef[n + 1]*wavef[n + 1];
	}

	sum = grid_step*sum/3.0;
	sum = sqrt(sum);

	for (int n = 0; n < grid_size; ++n)
	{
		wavef[n] = wavef[n]/sum;
	}
}

/******************************************************************************

 Function fgh_dvr(): return all eigenvalues/eigenvectors of a discrete variable
 representation (DVR) of the Fourier grid Hamiltonian (FGH), H = T + V, subject
 to a given single channel potential energy and reduced mass. As shown in Eq.
 (6a) and (6b) of Ref. [3].

 NOTE: The m-th eigenvector is the m-th column in a grid_size-by-grid_size row-
 major matrix: eigenvec[n*grid_size + m], where n = m = [0, grid_size).

******************************************************************************/

eigensystem *fgh_dvr(const int grid_size, const double grid_step,
                     const double mass, const double pot_energy[])
{
	ASSERT(pot_energy != NULL)

	double *h_matrix = allocate(grid_size*grid_size, false);

	register const double box_length
		= as_double(grid_size - 1)*grid_step;

	register const double factor
		= (M_PI*M_PI)/(mass*box_length*box_length);

	for (int n = 0; n < grid_size; ++n)
	{
		h_matrix[n*grid_size + n]
			= factor*as_double(grid_size*grid_size + 2)/6.0 + pot_energy[n];

		for (int m = (n + 1); m < grid_size; ++m)
		{
			register const double nm = as_double((n + 1) - (m + 1));
			register const double param = nm*M_PI/as_double(grid_size);

			h_matrix[n*grid_size + m]
				= pow(-1.0, nm)*factor/(sin(param)*sin(param));

			h_matrix[m*grid_size + n] = h_matrix[n*grid_size + m];
		}
	}

	return diag_symm_matrix(grid_size, h_matrix);
}

/******************************************************************************

 Function multich_fgh_dvr(): the same as fgh_dvr(), except that the potential
 energy is the one of a problem with max_ch channels within grid_size points:
 pot_energy[p].data[n*max_ch + m] where, p = [0, grid_size) and n = m = [0,
 max_ch). For details see Eq. (19a), Eq. (19b) and Eq. (19c) of Ref. [12].

 NOTE: in the multichannel problem there are (grid_size*max_ch) eigenvalues
 and the eigenvectors are stored in a (grid_size*max_ch)-by-(grid_size*max_ch)
 row-major matrix. Thus, each multichannel eigenvector will have one component
 in every channel of the problem following the order in which the channels (or
 blocks) are arranged in the original potential matrix representation.

******************************************************************************/

eigensystem *multich_fgh_dvr(const int max_ch,
                             const int grid_size,
                             const double grid_step,
                             const double mass,
                             const matrix pot_energy[])
{
	ASSERT(pot_energy != NULL)

	double *h_matrix = allocate(grid_size*max_ch*grid_size*max_ch, false);

	const double box_length = as_double(grid_size - 1)*grid_step;
	const double factor = (M_PI*M_PI)/(mass*box_length*box_length);

	int row = 0;
	for (int p = 0; p < max_ch; ++p)
	{
		for (int n = 0; n < grid_size; ++n)
		{
			int col = 0;
			for (int q = 0; q < max_ch; ++q)
			{
				for (int m = 0; m < grid_size; ++m)
				{
					if (n == m && p == q)
					{
						register const double param
							= factor*as_double(grid_size*grid_size + 2)/6.0;

						h_matrix[row*grid_size*max_ch + col]
							= param + pot_energy[n].data[p*max_ch + p];
					}

					if (n != m && p == q)
					{
						register const double nm
							= as_double((n + 1) - (m + 1));

						register const double param
							= nm*M_PI/as_double(grid_size);

						h_matrix[row*grid_size*max_ch + col]
							= pow(-1.0, nm)*factor/(sin(param)*sin(param));
					}

					if (n == m && p != q)
					{
						h_matrix[row*grid_size*max_ch + col]
							= pot_energy[n].data[p*max_ch + q];
					}

					if (n != m && p != q)
					{
						h_matrix[row*grid_size*max_ch + col] = 0.0;
					}

					++col;
				}
			}

			++row;
		}
	}

	ASSERT(row == grid_size*max_ch)
	return diag_symm_matrix(grid_size*max_ch, h_matrix);
}

/******************************************************************************

 Function fgh_dvr_wavef(): interpolate one eigenvector computed by fgh_dvr()
 using Eq. (4.1) of Ref. [5] and return its amplitude value at a given new r.

******************************************************************************/

double fgh_dvr_wavef(const int grid_size, const double eigenvec[],
                     const double r_min, const double r_max, const double r_new)
{
	register double result = 0.0;
	register const double r_inc = (r_max - r_min)/as_double(grid_size);

	for (int n = 0; n < grid_size; ++n)
	{
		register const double r_old
			= r_min + as_double(n)*r_inc;

		register const double param
			= M_PI*(r_new - r_old)/r_inc;

		result += (fabs(param) > 1.0e-7? eigenvec[n]*sinc(param) : eigenvec[n]);
	}

	return result;
}

/******************************************************************************

 Function fgh_dvr_product(): interpolate two eigenvectors computed by fgh_dvr()
 using Eq. (4.1) of Ref. [5] and return the product of their amplitude value at
 a given new r.

******************************************************************************/

double fgh_dvr_product(const int grid_size,
                       const double eigenvec_a[],
                       const double eigenvec_b[],
                       const double r_min,
                       const double r_max,
                       const double r_new)
{
	ASSERT(grid_size > 0)
	ASSERT(eigenvec_a != NULL)
	ASSERT(eigenvec_b != NULL)

	register double a = 0.0, b = 0.0;

	register const double r_inc
		= (r_max - r_min)/as_double(grid_size);

	for (int n = 0; n < grid_size; ++n)
	{
		register const double r_old
			= r_min + as_double(n)*r_inc;

		register const double param
			= M_PI*(r_new - r_old)/r_inc;

		if (fabs(param) > 1.0e-7)
		{
			a += eigenvec_a[n]*sinc(param);
			b += eigenvec_b[n]*sinc(param);
		}
		else
		{
			a += eigenvec_a[n];
			b += eigenvec_b[n];
		}
	}

	return a*b;
}

/******************************************************************************

 Function rot_integrand(): is an auxiliary routine which return the integrand
 for the multipolar expansion coefficients computed by qag_multipolar_coeff()
 at a given (r, R) Jacobi coordinate. See the inner integral in Eq. (22) of
 Ref. [6].

******************************************************************************/

struct multipolar_params
{
	const int lambda;
	jacobi_coor x;
};

static inline double rot_integrand(const double theta, void *params)
{
	struct multipolar_params *p = (struct multipolar_params *) params;

	p->x.theta = theta*180.0/M_PI;
	register const double pot_energy = pes(&p->x) - pec(p->x.arrang, p->x.r);

	/* Eq. (22) of Ref. [6], inner integrand */
	return pot_energy*gsl_sf_legendre_Pl(p->lambda, cos(theta))*sin(theta);
}

/******************************************************************************

 Function qag_multipolar_coeff(): return the multipolar expansion coefficient,
 V(lambda, r, R), of the PES at a given (r, R) Jacobi coordinate. See Eq. (4)
 of Ref. [1] and Eq. (22) of Ref. [6]. Integration over theta, in [0, pi],
 performed by the QAG algorithm.

******************************************************************************/

double qag_multipolar_coeff(const char arrang,
                            const int lambda, const double r, const double R)
{
	gsl_integration_workspace *work
		= gsl_integration_workspace_alloc(GSL_MAX_WORKSPACE);

	ASSERT(work != NULL)

	struct multipolar_params p =
	{
		.lambda = lambda,
		.x.arrang = arrang,
		.x.r = r,
		.x.R = R
	};

	gsl_function f =
	{
		.params = &p,
		.function = &rot_integrand
	};

	double error = 0.0, result = 0.0;

	gsl_integration_qag(&f, 0.0, M_PI, 0.1, 1.0e-12, GSL_MAX_WORKSPACE,
	                         GSL_INTEG_GAUSS61, work, &result, &error);

	gsl_integration_workspace_free(work);
	return as_double(2*lambda + 1)*result/2.0;
}

/******************************************************************************

 Function vib_integrand(): is an auxiliary routine which return the integrand
 for the vibrational integral computed by qag_vib_quadr() at a given R Jacobi
 coordinate. See the outer integral in Eq. (22) of Ref. [6].

 NOTE: any factor r^2 of Eq. (22) is not included as it is absorbed in the
 definition of vibrational wavefunctions.

******************************************************************************/

struct vib_quadr_params
{
	const char arrang;
	const int lambda, grid_size;
	const double r_min, r_max, *wavef_a, *wavef_b, R;
};

static inline double vib_integrand(const double r, void *params)
{
	struct vib_quadr_params *p = (struct vib_quadr_params *) params;

	register const double wavef_ab
		= fgh_dvr_product(p->grid_size, p->wavef_a, p->wavef_b, p->r_min, p->r_max, r);

	register const double rot_integral
		= qag_multipolar_coeff(p->arrang, p->lambda, r, p->R);

	/* Eq. (22) of Ref. [6], outer integrand */
	return wavef_ab*rot_integral;
}

/******************************************************************************

 Function qag_vib_quadr(): return the lambda-dependent vibrational coefficient,
 V(lambda, R), as function of the R Jacobi coordinate, for the expansion of the
 atom-diatom coupling potential, as shown in Eq. (15) of Ref. [2]. Integration
 over r, in [r_min, r_max], performed by the QAG algorithm.

 NOTE: the subscript n of Eq. (15) is written lambda here and, likewise,
 variables (r, R) are interchanged with respect to Ref. [2].

 NOTE: r_min and r_max are expected the same for both vib. wavefunctions.

******************************************************************************/

double qag_vib_quadr(const char arrang,
                     const int lambda,
                     const int grid_size,
                     const double wavef_a[],
                     const double wavef_b[],
                     const double r_min,
                     const double r_max,
                     const double R)
{
	gsl_integration_workspace *work
		= gsl_integration_workspace_alloc(GSL_MAX_WORKSPACE);

	ASSERT(work != NULL)

	struct vib_quadr_params p =
	{
		.R = R,
		.r_min = r_min,
		.r_max = r_max,
		.lambda = lambda,
		.arrang = arrang,
		.wavef_a = wavef_a,
		.wavef_b = wavef_b,
		.grid_size = grid_size
	};

	gsl_function f =
	{
		.params = &p,
		.function = &vib_integrand
	};

	double error = 0.0, result = 0.0;

	gsl_integration_qag(&f, r_min, r_max, 0.1, 1.0e-12, GSL_MAX_WORKSPACE,
	                            GSL_INTEG_GAUSS61, work, &result, &error);

	gsl_integration_workspace_free(work);
	return result;
}

/******************************************************************************

 Function wigner_3j():

******************************************************************************/

double wigner_3j(const int a, const int b, const int c,
                 const int d, const int e, const int f)
{
	return gsl_sf_coupling_3j(2*a, 2*b, 2*c, 2*d, 2*e, 2*f);
}

/******************************************************************************

 Function wigner_6j():

******************************************************************************/

double wigner_6j(const int a, const int b, const int c,
                 const int d, const int e, const int f)
{
	return gsl_sf_coupling_6j(2*a, 2*b, 2*c, 2*d, 2*e, 2*f);
}

/******************************************************************************

 Function percival_seaton(): return the so-called Percival & Seaton term often
 used in the definition of the atom-diatom collisional coupling matrix. Where,

 j1     = diatomic rotational angular momentum quantum number for channel 1
 j2     = diatomic rotational angular momentum quantum number for channel 2
 l1     = atom-diatom orbital angular momentum quantum number for channel 1
 l2     = atom-diatom orbital angular momentum quantum number for channel 2
 J      = total angular momentum of the problem
 lambda = an integer parameter (0, 1, 2, ...)

 and also depends parametrically on the electronic spin multiplicity (1 for
 singlet, 2 for doublet and 3 for triplet).

 See Ref. [7] for more, in particular, Eq. (A1).

******************************************************************************/

double percival_seaton(const int spin_mult, const int j1, const int j2,
                       const int l1, const int l2, const int lambda, const int J)
{
	register double result = 0.0;

	switch (spin_mult)
	{
		case 1: /* Eq. (3) of Ref. [8] */
			result  = pow(-1.0, j1 + j2 - J);
			result *= wigner_3j(l1, l2, lambda, 0, 0, 0);
			result *= wigner_3j(j1, j2, lambda, 0, 0, 0);
			result *= wigner_6j(j1, j2, lambda, l2, l1, J);
			result *= sqrt(as_double(2*j1 + 1)*as_double(2*j2 + 1)*as_double(2*l1 + 1)*as_double(2*l2 + 1));
		break;

		case 2: /* Eq. (C.21) of Ref. [9] */
		break;

		case 3: /* Eq. (C.22) of Ref. [9] */
		break;

		default:
			PRINT_ERROR("invalid spin multiplicity = %d\n", spin_mult)
	}

	return result;
}

/******************************************************************************
******************************************************************************/

double asympt_scatt_wavef(const char mode, const double mass,
                          const double kin_energy, const double R)
{
	const double wavenum = sqrt(2.0*mass*kin_energy);

	switch (mode)
	{
		case 'i':
			if (kin_energy > 0.0)
				return sin(wavenum*R)/sqrt(wavenum/mass);
			else
				return exp(fabs(wavenum)*R)/sqrt(2.0*wavenum/mass);

		case 'o':
			if (kin_energy > 0.0)
				return cos(wavenum*R)/sqrt(wavenum/mass);
			else
				return exp(-fabs(wavenum)*R)/sqrt(2.0*wavenum/mass);

		default:
			PRINT_ERROR("invalid mode = %c\n", mode)
			return 0.0;
	}
}

/******************************************************************************

 Function renorm_numerov(): use the method of B. R. Johnson, Ref. [11], to cons
 truct the wavefunction at a given trial energy subject to a single channel pot
 ential energy and reduced mass. The return pointer contains the unormalized am
 plitude in a grid_size-long vector.

 NOTE: on exit, a negative number of nodes implies that no matching point has
 been found (see Ref. [11]).

******************************************************************************/

double *renorm_numerov(const int grid_size,
                       const double grid_step,
                       const double pot_energy[],
                       const double trial_energy,
                       const double mass,
                       double *error,
                       int *nodes)
{
	ASSERT(grid_size > 0)
	ASSERT(pot_energy != NULL)

	double *T = allocate(grid_size, false);
	double *F = allocate(grid_size, false);

	double *inward_R  = allocate(grid_size, false);
	double *outward_R = allocate(grid_size, false);

/*
 *	Solve Eq. (32), (35) and (46) for each T, U and R (inward) coefficient:
 */

	inward_R[grid_size - 1] = 1.0E20;

	T[grid_size - 1] =
	-pow(grid_step, 2)*2.0*mass*(trial_energy - pot_energy[grid_size - 1])/12.0;

	int M = 0;
	for (int n = (grid_size - 2); n > -1; --n)
	{
		T[n] = -pow(grid_step, 2)*2.0*mass*(trial_energy - pot_energy[n])/12.0;
		inward_R[n] = (2.0 + 10.0*T[n])/(1.0 - T[n]) - 1.0/inward_R[n + 1];

/*
 *		Check if the matching point, M, is found:
 */

		if (inward_R[n] <= -1.0)
		{
			M = n;
			break;
		}
	}

	if (M == 0)
	{
		free(T);
		free(F);
		free(inward_R);
		free(outward_R);

		*error =  1.0E20;
		*nodes = -1;
		return NULL;
	}

/*
 *	Solve Eq. (32), (35) and (37) for each T, U and R (outward) coefficient:
 */

	outward_R[0] = 1.0E20;

	T[0] = -pow(grid_step, 2)*2.0*mass*(trial_energy - pot_energy[0])/12.0;

	*nodes = 0;
	for (int n = 1; n < M; ++n)
	{
		T[n] = -pow(grid_step, 2)*2.0*mass*(trial_energy - pot_energy[n])/12.0;
		outward_R[n] = (2.0 + 10.0*T[n])/(1.0 - T[n]) - 1.0/outward_R[n - 1];

/*
 *		Count the number of nodes in the wavefunction:
 */

		if (outward_R[n] < 0.0) ++(*nodes);
	}

/*
 *	Solve Eq. (36) and (45) for each F coefficient:
 */

	F[M] = 1.0;
	for (int n = (M - 1); n > -1; --n)
	{
		F[n] = F[n + 1]/outward_R[n];
	}

	for (int n = (M + 1); n < grid_size; ++n)
	{
		F[n] = F[n - 1]/inward_R[n];
	}

/*
 *	Solve Eq. (33) for the wavefunction (unnormalized):
 */

	double *wavef = allocate(grid_size, false);

	for (int n = 0; n < grid_size; ++n)
	{
		wavef[n] = F[n]/(1.0 - T[n]);
	}

/*
 *	Solve Eq. (48) for the difference, D (named as error here), between the
 *	left and right integrations:
 */

	*error = (0.5 - T[M + 1])*(1.0/inward_R[M + 1] - outward_R[M])/(1.0 - T[M + 1]) -
		     (0.5 - T[M - 1])*(inward_R[M] - 1.0/outward_R[M - 1])/(1.0 - T[M - 1]);

	*error *= (1.0 - T[M]);

	free(T);
	free(F);
	free(inward_R);
	free(outward_R);

	return wavef;
}

/******************************************************************************

 Function multich_renorm_numerov(): use the multichannel renormalized Numerov
 algorithm of B. R. Johnson, Ref. [4], in order to propagate the multichannel
 wavefunction ratio matrix, R, from the radial grid point (n - 1) to n driven
 by the interaction potential V and subjected to a given total energy. See Eq.
 (22), (23) and (24) of Ref. [4]. The potential matrix is used as an auxiliar
 chunk of memory and, therefore, its content is destroyed on exit.

******************************************************************************/

void multich_renorm_numerov(const int max_ch, const double mass, const double grid_step,
                            const double tot_energy, double v[], double r[], const bool use_omp)
{
	ASSERT(v != NULL)
	ASSERT(r != NULL)
	ASSERT(max_ch > 0)

	double *w = v;
	#pragma omp parallel sections num_threads(2) if(use_omp)
	{
/*		#pragma omp section
		{
			if (!is_null(max_ch*max_ch, r)) inverse_symm_matrix(max_ch, r);
		}*/

		#pragma omp section
		{
/*
 *			NOTE: From Eq. (2) and (17) of Ref. [4] the following numerical
 *			factor is defined in atomic units:
 */

			register const double factor
				= -grid_step*grid_step*2.0*mass/12.0;

/*
 *			Resolve Eq. (23) of Ref. [4] with Eq. (2) and (17) plugged in:
 */

			for (int n = 0; n < max_ch; ++n)
			{
				w[n*max_ch + n]
					= 1.0 - factor*(tot_energy - v[n*max_ch + n]);

				for (int m = (n + 1); m < max_ch; ++m)
				{
					w[n*max_ch + m] = factor*v[n*max_ch + m];
					w[m*max_ch + n] = factor*v[m*max_ch + n];
				}
			}

			inverse_symm_matrix(max_ch, w);
		}
	}

/*
 *	Solve Eq. (22) and (24) of Ref. [4]:
 */

	for (int n = 0; n < max_ch; ++n)
	{
		r[n*max_ch + n] = (12.0*w[n*max_ch + n] - 10.0) - r[n*max_ch + n];

		for (int m = (n + 1); m < max_ch; ++m)
		{
			r[n*max_ch + m] = 12.0*w[n*max_ch + m] - r[n*max_ch + m];
			r[m*max_ch + n] = 12.0*w[m*max_ch + n] - r[m*max_ch + n];
		}
	}

	w = NULL;
}

/******************************************************************************

 Function johnson_numerov(): use the multichannel renormalized Numerov
 algorithm of B. R. Johnson, Ref. [4], in order to propagate the multichannel
 ratio matrix, R, from the radial grid point (n - 1) to n driven by the inter
 action potential V and subjected to a given total energy. See Eq. (22), (23)
 and (24) of Ref. [4]. The potential matrix is used as workspace and its
 content is destroyed on exit.

******************************************************************************/

void johnson_numerov(const double grid_step, const double mass,
                     const double tot_energy, matrix *v, matrix *r)
{
	ASSERT(v != NULL)
	ASSERT(r != NULL)

	if (!matrix_null(r))
	{
		matrix_inv(r);
	}

/*
 *	NOTE: From Eq. (2) and (17) of Ref. [4] the following numerical
 *	factor is defined in atomic units:
 */

	register const double factor
		= -grid_step*grid_step*2.0*mass/12.0;

/*
 *	Resolve Eq. (23) of Ref. [4] with Eq. (2) and (17) plugged in:
 */

	matrix *w = v;

	for (int n = 0; n < matrix_row(v); ++n)
	{
		register double t = factor*(tot_energy - matrix_get(v, n, n));
		matrix_set(w, n, n, 1.0 - t);

		for (int m = (n + 1); m < matrix_col(v); ++m)
		{
			t = factor*(0.0 - matrix_get(v, n, m));
			matrix_set(w, n, m, 0.0 - t);

			t = factor*(0.0 - matrix_get(v, m, n));
			matrix_set(w, m, n, 0.0 - t);
		}
	}

/*
 *	Solve Eq. (22) and (24) of Ref. [4]:
 */

	matrix_inv(w);

	for (int n = 0; n < matrix_row(v); ++n)
	{
		register double u = 12.0*matrix_get(w, n, n) - 10.0;
		matrix_set(r, n, n, u - matrix_get(r, n, n));

		for (int m = (n + 1); m < matrix_col(v); ++m)
		{
			u = 12.0*matrix_get(w, n, m);
			matrix_set(r, n, m, u - matrix_get(r, n, m));

			u = 12.0*matrix_get(w, m, n);
			matrix_set(r, m, n, u - matrix_get(r, m, n));
		}
	}

	w = NULL;
}

/******************************************************************************

 Function johnson_logd(): use the algorithm of B. R. Johnson, Ref. [14], in
 order to propagate the multichannel log derivative matrix, Y, from the radial
 grid point (n - 1) to n driven by the interaction potential V.

 NOTE: matrix V is defined as Eq. (2) of Ref. [14].

******************************************************************************/

void johnson_logd(const int n, const int grid_size,
                  const double grid_step, const matrix *v, matrix *y)
{
	ASSERT(n > 0)
	ASSERT(grid_size > 0)
	ASSERT(grid_size >= n)
	ASSERT(GSL_IS_ODD(grid_size) == 1)

	ASSERT(v != NULL)
	ASSERT(y != NULL)

	/* NOTE: z[1] = left term and u[1] = right term of Eq. (6). */
	matrix *z = matrix_alloc(2, matrix_row(v), matrix_col(v), false);
	matrix *u = matrix_alloc(2, matrix_row(v), matrix_col(v), false);

/*
 *	Solve Eq. (7) of Ref. [14]:
 */

	if (GSL_IS_EVEN(n) == 1)
	{
		register const double w = (n == 0? 1.0 : 2.0);

		for (int p = 0; p < matrix_row(v); ++p)
		{
			matrix_set(&z[0], p, p, 1.0 + matrix_get(y, p, p)*grid_step);
			matrix_set(&u[1], p, p, matrix_get(v, p, p)*grid_step*w/3.0);

			for (int q = (p + 1); q < matrix_col(v); ++q)
			{
				matrix_set(&z[0], p, q, 0.0 + matrix_get(y, p, q)*grid_step);
				matrix_set(&u[1], p, q, matrix_get(v, p, q)*grid_step*w/3.0);

				matrix_set(&z[0], q, p, 0.0 + matrix_get(y, q, p)*grid_step);
				matrix_set(&u[1], q, p, matrix_get(v, q, p)*grid_step*w/3.0);
			}
		}
	}
	else
	{
		register const double w = (n == grid_size? 1.0 : 4.0);

		for (int p = 0; p < matrix_row(v); ++p)
		{
			matrix_set(&z[0], p, p, 1.0 + matrix_get(y, p, p)*grid_step);
			matrix_set(&u[0], p, p, 1.0 + matrix_get(v, p, p)*grid_step*grid_step/6.0);

			for (int q = (p + 1); q < matrix_col(v); ++q)
			{
				matrix_set(&z[0], p, q, 0.0 + matrix_get(y, p, q)*grid_step);
				matrix_set(&u[0], p, q, 0.0 + matrix_get(v, p, q)*grid_step*grid_step/6.0);

				matrix_set(&z[0], q, p, 0.0 + matrix_get(y, q, p)*grid_step);
				matrix_set(&u[0], q, p, 0.0 + matrix_get(v, q, p)*grid_step*grid_step/6.0);
			}
		}

		matrix_inv(&u[0]);
		matrix_mul(grid_step*w/3.0, &u[0], v, 0.0, &u[1]);
	}

	matrix_inv(&z[0]);
	matrix_mul(1.0, &z[0], y, 0.0, &z[1]);

/*
 *	Solve Eq. (6) of Ref. [14]:
 */

	matrix_sub(1.0, &z[1], 1.0, &u[1], y, matrix_row(v) > 100);

	matrix_free_all(2, z);
	matrix_free_all(2, u);
}

/******************************************************************************

 Function react_matrix(): build the reactant matrix K as described by Eq. (A19)
 of Ref. [4], from the ratio matrix of the multichannel scattering wavefunction
 for a given total energy, reduced mass and asymptotic value of R.

******************************************************************************/

double *react_matrix(const int max_ch,
                     const int l[],
                     const double mass,
                     const double grid_step,
                     const double tot_energy,
                     const double level[],
                     const double r[],
                     const double R)
{
	ASSERT(l != NULL)
	ASSERT(r != NULL)
	ASSERT(max_ch > 0)
	ASSERT(level != NULL)

/*
 *	NOTE: from Eq. (2) and (17) of Ref. [4] the following numerical factor is
 *	defined in atomic units:
 */

	register const double factor
		= -grid_step*grid_step*2.0*mass/12.0;

/*
 *	Step 1: build the j(R) and n(R) diagonal matrices defined by Eq. (A16) and
 *	(A17) of Ref. [4]; where, Eq. (23) is also used:
 */

	double *n = allocate(max_ch*max_ch, true);
	double *j = allocate(max_ch*max_ch, true);

	for (int m = 0; m < max_ch; ++m)
	{
		if (level[m] < tot_energy)
		{
			register const double wavenum
				= sqrt(2.0*mass*(tot_energy - level[m]));

			register const double w
				= 1.0 - factor*(tot_energy - centr_term(l[m], mass, R));

			n[m*max_ch + m] = w*riccati_bessel('n', l[m], wavenum, R);
			j[m*max_ch + m] = w*riccati_bessel('j', l[m], wavenum, R);
		}
		else
		{
			n[m*max_ch + m] = 1.0;
			j[m*max_ch + m] = 1.0;
		}
	}

/*
 *	Step 2: build the product rn and rj at the grid point R, as shown in Eq.
 *	(A19) of Ref. [4]:
 */

	double *rn = allocate(max_ch*max_ch, false);
	double *rj = allocate(max_ch*max_ch, false);

	call_dgemm('n', 'n', max_ch, max_ch, max_ch,
	           1.0, r, max_ch, n, max_ch, 0.0, rn, max_ch);

	call_dgemm('n', 'n', max_ch, max_ch, max_ch,
	           1.0, r, max_ch, j, max_ch, 0.0, rj, max_ch);

/*
 *	Step 3: subtract j(R + grid_step) and n(R + grid_step) diagonal matrices
 *	from the rn(R) and rj(R) products, following Eq. (A19) of Ref. [4]:
 */

	for (int m = 0; m < max_ch; ++m)
	{
		register const double wavenum
			= sqrt(2.0*mass*(tot_energy - level[m]));

		register const double w_prime
			= 1.0 - factor*(tot_energy - centr_term(l[m], mass, R + grid_step));

		register const double n_prime
			= w_prime*riccati_bessel('n', l[m], wavenum, R + grid_step);

		register const double j_prime
			= w_prime*riccati_bessel('j', l[m], wavenum, R + grid_step);

		if (level[m] < tot_energy)
		{
			rn[m*max_ch + m] -= n_prime;
			rj[m*max_ch + m] -= j_prime;
		}
		else
		{
			register const double w
				= 1.0 - factor*(tot_energy - centr_term(l[m], mass, R));

			rn[m*max_ch + m] -= n_prime/(w*modif_spher_bessel('n', l[m], wavenum, R));
			rj[m*max_ch + m] -= j_prime/(w*modif_spher_bessel('j', l[m], wavenum, R));
		}
	}

	inverse_symm_matrix(max_ch, rn);

/*
 *	Step 4: resolve Eq. (A19) of Ref. [4] for the K-matrix:
 */

	double *k = allocate(max_ch*max_ch, false);

	call_dgemm('n', 'n', max_ch, max_ch, max_ch,
	           1.0, rn, max_ch, rj, max_ch, 0.0, k, max_ch);

	for (int m = 0; m < max_ch*max_ch; ++m)
	{
		k[m] = -k[m];
	}

	free(j);
	free(n);
	free(rn);
	free(rj);
	return k;
}

/******************************************************************************

 Function scatt_matrix(): return both imaginary and real parts of a scattering
 matrix, S, defined as

 S := (I + iK)inv(I - iK)
    = (I - KK)inv(I + KK) + i(2K)inv(I + KK)

 Where, I is the unit matrix and K is a reactant matrix for max_ch channels,
 while open channels are defined with respect each asymptotic energy (level)
 and a given total energy.

 For details see Eq. (20) of Ref. [10].

******************************************************************************/

s_matrix *scatt_matrix(const int max_ch, const double levels[],
                       const double k[], const double tot_energy)
{
	ASSERT(max_ch > 0)
	ASSERT(k != NULL)
	ASSERT(levels != NULL)

/*
 *	Step 1: count the number of open channels:
 */

	register int open_ch = 0;
	for (int n = 0; n < max_ch; ++n)
	{
		if (levels[n] < tot_energy) ++open_ch;
	}

	if (open_ch < 1)
	{
		PRINT_ERROR("no open channels for total energy = %f\n", tot_energy)
		return NULL;
	}

/*
 *	Step 2: build the open-open block of K, K_OO:
 */

	double *k_oo = allocate(open_ch*open_ch, false);

	int p = 0;
	for (int n = 0; n < max_ch; ++n)
	{
		if (levels[n] > tot_energy) continue;

		int q = 0;
		for (int m = 0; m < max_ch; ++m)
		{
			if (levels[m] > tot_energy) continue;

			k_oo[p*open_ch + q] = k[n*max_ch + m];
			++q;
		}

		++p;
	}

	ASSERT(open_ch == p)

/*
 *	Step 3: compute A = KK, B = inv(I + KK) and C = (I - KK):
 */

	double *a = allocate(open_ch*open_ch, false);
	double *b = allocate(open_ch*open_ch, false);
	double *c = allocate(open_ch*open_ch, false);

	call_dgemm('n', 'n', open_ch, open_ch, open_ch,
	           1.0, k_oo, open_ch, k_oo, open_ch, 0.0, a, open_ch);

	for (int n = 0; n < open_ch; ++n)
	{
		b[n*open_ch + n] = 1.0 + a[n*open_ch + n];
		c[n*open_ch + n] = 1.0 - a[n*open_ch + n];

		for (int m = (n + 1); m < open_ch; ++m)
		{
			b[n*open_ch + m] =  a[n*open_ch + m];
			c[n*open_ch + m] = -a[n*open_ch + m];

			b[m*open_ch + n] =  a[m*open_ch + n];
			c[m*open_ch + n] = -a[m*open_ch + n];
		}
	}

	inverse_symm_matrix(open_ch, b);

/*
 *	Step 4: build the S matrix, Re(S) = C*inv(B) and Im(S) = 2*K*inv(B):
 */

	s_matrix *s = calloc(1, sizeof(s_matrix));
	ASSERT(s != NULL)

	s->re_part = allocate(open_ch*open_ch, false);
	s->im_part = allocate(open_ch*open_ch, false);

	call_dgemm('n', 'n', open_ch, open_ch, open_ch, 1.0, c,
	           open_ch, b, open_ch, 0.0, s->re_part, open_ch);

	call_dgemm('n', 'n', open_ch, open_ch, open_ch, 2.0, k_oo,
	           open_ch, b, open_ch, 0.0, s->im_part, open_ch);

	s->max_ch = open_ch;

	free(a);
	free(b);
	free(c);
	free(k_oo);
	return s;
}

/******************************************************************************

 Function r_to_k(): build the reactant matrix K as described by Eq. (A19) of
 Ref. [4], from the ratio matrix of the multichannel scattering wavefunction
 for a given total energy, reduced mass and asymptotic value of R.

******************************************************************************/

matrix *r_to_k(const scatt_param *p,
               const matrix *r, const double grid_step, const double R)
{
	ASSERT(p != NULL)
	ASSERT(r != NULL)
	ASSERT(p->max_ch == matrix_row(r))
	ASSERT(p->max_ch == matrix_col(r))

/*
 *	NOTE: from Eq. (2) and (17) of Ref. [4] the following numerical factor is
 *	defined in atomic units:
 */

	register const double factor
		= -grid_step*grid_step*2.0*p->mass/12.0;

/*
 *	Step 1: build the j(R) and n(R) diagonal matrices defined by Eq. (A16) and
 *	(A17) of Ref. [4]; where, Eq. (23) is also used:
 */

	matrix *n = matrix_alloc(1, p->max_ch, p->max_ch, true);
	matrix *j = matrix_alloc(1, p->max_ch, p->max_ch, true);

	for (int m = 0; m < p->max_ch; ++m)
	{
		if (p->kin_energy[m] > 0.0)
		{
			register const double w
				= 1.0 - factor*(p->tot_energy - centr_term(p->l[m], p->mass, R));

			matrix_set(n, m, m, w*riccati_bessel('n', p->l[m], p->wavenum[m], R));
			matrix_set(j, m, m, w*riccati_bessel('j', p->l[m], p->wavenum[m], R));
		}
		else
		{
			matrix_set(n, m, m, 1.0);
			matrix_set(j, m, m, 1.0);
		}
	}

/*
 *	Step 2: build the product rn and rj at the grid point R, as shown in Eq.
 *	(A19) of Ref. [4]:
 */

	matrix *rn = matrix_alloc(1, p->max_ch, p->max_ch, false);
	matrix *rj = matrix_alloc(1, p->max_ch, p->max_ch, false);

	matrix_mul(1.0, r, n, 0.0, rn);
	matrix_mul(1.0, r, j, 0.0, rj);

/*
 *	Step 3: subtract j(R + grid_step) and n(R + grid_step) diagonal matrices
 *	from the rn(R) and rj(R) products, following Eq. (A19) of Ref. [4]:
 */

	for (int m = 0; m < p->max_ch; ++m)
	{
		register const double w_prime
			= 1.0 - factor*(p->tot_energy - centr_term(p->l[m], p->mass, R + grid_step));

		register const double n_prime
			= w_prime*riccati_bessel('n', p->l[m], p->wavenum[m], R + grid_step);

		register const double j_prime
			= w_prime*riccati_bessel('j', p->l[m], p->wavenum[m], R + grid_step);

		if (p->kin_energy[m] > 0.0)
		{
			matrix_decr(rn, m, m, n_prime);
			matrix_decr(rj, m, m, j_prime);
		}
		else
		{
			register const double w
				= 1.0 - factor*(p->tot_energy - centr_term(p->l[m], p->mass, R));

			matrix_decr(rn, m, m, n_prime/(w*modif_spher_bessel('n', p->l[m], p->wavenum[m], R)));
			matrix_decr(rj, m, m, j_prime/(w*modif_spher_bessel('j', p->l[m], p->wavenum[m], R)));
		}
	}

/*
 *	Step 4: resolve Eq. (A19) of Ref. [4] for the K-matrix:
 */

	matrix_inv(rn);

	matrix *k = matrix_alloc(1, p->max_ch, p->max_ch, false);
	matrix_mul(-1.0, rn, rj, 0.0, k);

	matrix_free(j);
	matrix_free(n);
	matrix_free(rn);
	matrix_free(rj);
	return k;
}

/******************************************************************************

 Function k_to_s(): return both real, s[0], and imaginary s[1], parts of a
 scattering matrix, S, defined as

 S := (I + iK)inv(I - iK)
    = (I - KK)inv(I + KK) + i(2K)inv(I + KK)

 Where, I is the unit matrix and K is the open-open block of a reactant matrix.

 For details see Eq. (20) of Ref. [10].

******************************************************************************/

matrix *k_to_s(const matrix *k)
{
	matrix *a = matrix_alloc(1, matrix_row(k), matrix_col(k), false);
	matrix *b = matrix_alloc(1, matrix_row(k), matrix_col(k), false);
	matrix *c = matrix_alloc(1, matrix_row(k), matrix_col(k), false);

/*
 *	Resolve A = KK, B = (I + A) and C = (I - A):
 */

	matrix_mul(1.0, k, k, 0.0, a);

	for (int n = 0; n < matrix_row(k); ++n)
	{
		matrix_set(b, n, n, 1.0 + matrix_get(a, n, n));
		matrix_set(c, n, n, 1.0 - matrix_get(a, n, n));

		for (int m = (n + 1); m < matrix_col(k); ++m)
		{
			matrix_set(b, n, m, 0.0 + matrix_get(a, n, m));
			matrix_set(b, m, n, 0.0 + matrix_get(a, m, n));

			matrix_set(c, n, m, 0.0 - matrix_get(a, n, m));
			matrix_set(c, m, n, 0.0 - matrix_get(a, m, n));
		}
	}

/*
 *	Build the S matrix, Re(S) = C*inv(B) and Im(S) = 2*K*inv(B):
 */

	matrix_inv(b);

	matrix *s = matrix_alloc(2, matrix_row(k), matrix_col(k), false);

	matrix_mul(1.0, c, b, 0.0, &s[0]);
	matrix_mul(2.0, k, b, 0.0, &s[1]);

	matrix_free(a);
	matrix_free(b);
	matrix_free(c);
	return s;
}

/******************************************************************************

 Function olson_smith_model(): return the 2-by-2 model potential of Olson and
 Smith for the system Ne + He^+ at a given internuclear distance x. See Eq.
 (45) and Table I from Ref. [13] for details.

 NOTE: a handy benchmark for the algorithms.

******************************************************************************/

matrix *olson_smith_model(const double x)
{
	matrix *m = matrix_alloc(1, 2, 2, false);

	matrix_set(m, 0, 0, 21.1*exp(-x/0.678)/x);
	matrix_set(m, 0, 1, 0.170*exp(-x/0.667));
	matrix_set(m, 1, 0, 0.170*exp(-x/0.667));
	matrix_set(m, 1, 1, (21.1/x - 12.1)*exp(-x/0.678) + 16.8/27.2113961);

	return m;
}

/******************************************************************************

 Function tully_model(): return one of the three (n = 1, 2, 3) 2-by-2 model
 potentials of J. C. Tully at a given x. See Eq. (21), (23) and (24) of
 Ref. [15].

******************************************************************************/

matrix *tully_model(const int n, const double x)
{
	matrix *m = matrix_alloc(1, 2, 2, false);

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

			if (x < 0.0)
			{
				matrix_set(m, 0, 1, 0.10*exp(0.90*x));
			}
			else
			{
				matrix_set(m, 0, 1, 0.10*(2.0 - exp(-0.90*x)));
			}

			matrix_set(m, 1, 0, 0.0);
		break;

		default:
			PRINT_ERROR("invalid index n = %d\n", n)
			exit(EXIT_FAILURE);
	}

	return m;
}
