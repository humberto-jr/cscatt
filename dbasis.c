#include "modules/pes.h"
#include "modules/mass.h"
#include "modules/file.h"
#include "modules/matrix.h"
#include "modules/globals.h"

#include "mass_config.h"
#include "basis_config.h"

#define FORMAT "# %5d   %5d   %5d   %5d   %+5d     % -8e  % -8e  % -8e\n"

/******************************************************************************

 Function parity(): return the parity p = (-1)^l for a given l; either p = +1
 or p = -1.

******************************************************************************/

inline static int parity(const int l)
{
	return (int) pow(-1.0, l);
}

/******************************************************************************

 Function centr_term(): returns the actual centrifugal potential term at x for
 a given angular momentum l, provided the mass under consideration.

******************************************************************************/

inline static double centr_term(const int l,
                                const double mass, const double x)
{
	return as_double(l*(l + 1))/(2.0*mass*x*x);
}

/******************************************************************************

 Function fgh_matrix(): return the discrete variable representation (DVR) of a
 Fourier grid Hamiltonian (FGH) subjected to a given single channel potential
 energy and reduced mass. As shown in Eq. (6a) and (6b) of Ref. [1].

 NOTE: The m-th eigenvector is the m-th column in a grid_size-by-grid_size row-
 major matrix: eigenvec[n*grid_size + m], where n = m = [0, grid_size).

******************************************************************************/

matrix *fgh_matrix(const int grid_size,
                   const double grid_step,
                   const double pot_energy[],
                   const double mass)
{
	ASSERT(pot_energy != NULL)

	matrix *result
		= matrix_alloc(grid_size, grid_size, false);

	const double box_length
		= as_double(grid_size - 1)*grid_step;

	const double factor
		= (M_PI*M_PI)/(mass*box_length*box_length);

	const double nn_term
		= factor*as_double(grid_size*grid_size + 2)/6.0;

	for (int n = 0; n < grid_size; ++n)
	{
		matrix_diag_set(result, n, nn_term + pot_energy[n]);

		for (int m = (n + 1); m < grid_size; ++m)
		{
			const double nm
				= as_double((n + 1) - (m + 1));

			const double nm_term
				= sin(nm*M_PI/as_double(grid_size));

			matrix_symm_set(result, n, m, pow(-1.0, nm)*factor/pow(nm_term, 2));
		}
	}

	return result;
}

/******************************************************************************

 Function norm(): normalize to unity a radial eigenvector of grid_size points,
 using a 1/3-Simpson quadrature rule.

******************************************************************************/

void norm(const int grid_size,
          const double grid_step, const bool use_omp, double eigenvec[])
{
	ASSERT(grid_size > 0)
	ASSERT(eigenvec != NULL)

	const int n_max
		= (grid_size%2 == 0? grid_size : grid_size - 1);

	double sum
		= eigenvec[0]*eigenvec[0] + eigenvec[n_max - 1]*eigenvec[n_max - 1];

	#pragma omp parallel for default(none) shared(eigenvec) reduction(+:sum) if(use_omp)
	for (int n = 1; n < (n_max - 2); n += 2)
	{
		sum += 4.0*eigenvec[n]*eigenvec[n];
		sum += 2.0*eigenvec[n + 1]*eigenvec[n + 1];
	}

	sum = grid_step*sum/3.0;
	sum = sqrt(sum);

	#pragma omp parallel for default(none) shared(eigenvec, sum) if(use_omp)
	for (int n = 0; n < grid_size; ++n)
	{
		eigenvec[n] = eigenvec[n]/sum;
	}
}

int main(int argc, char *argv[])
{
	ASSERT(argc > 1)
	file_init_stdin(argv[1]);

/*
 *	Total angular momentum, J:
 */

	const int J
		= (int) file_get_key(stdin, "J", 0.0, INF, 0.0);

	const int J_parity
		= (int) file_get_key(stdin, "parity", -1.0, 1.0, 0.0);

/*
 *	Vibrational quantum numbers, v:
 */

	const int v_min
		= (int) file_get_key(stdin, "v_min", 0.0, INF, 0.0);

	const int v_max
		= (int) file_get_key(stdin, "v_max", as_double(v_min), INF, as_double(v_min));

	const int v_step
		= (int) file_get_key(stdin, "v_step", 1.0, INF, 1.0);

/*
 *	Rotational quantum numbers, j:
 */

	const int j_min
		= (int) file_get_key(stdin, "j_min", 0.0, INF, 0.0);

	const int j_max
		= (int) file_get_key(stdin, "j_max", as_double(j_min), INF, as_double(j_min));

	const int j_step
		= (int) file_get_key(stdin, "j_step", 1.0, INF, 1.0);

/*
 *	Vibrational grid:
 */

	const int grid_size
		= (int) file_get_key(stdin, "rovib_grid_size", as_double(v_max + 1), INF, 500.0);

	const double r_min
		= file_get_key(stdin, "r_min", 0.0, INF, 0.5);

	const double r_max
		= file_get_key(stdin, "r_max", r_min, INF, r_min + 30.0);

	const double grid_step
		= (r_max - r_min)/as_double(grid_size);

/*
 *	Arrangement (1 == a, 2 == b, 3 == c) and atomic masses:
 */

	const char arrang
		= 96 + (int) file_get_key(stdin, "arrang", 1.0, 3.0, 1.0);

	const enum mass_case m
		= init_atomic_masses(stdin, arrang, 'p', 0);

	pes_init();

/*
 *	Resolve the diatomic eigenvalue for each j-case and sort results as scatt. channels:
 */

	printf("#\n");
	printf("# J = %d, J-parity = %d\n", J, J_parity);
	printf("#   Ch.       v       j       l       p        E (a.u.)       E (cm-1)        E (eV)   \n");
	printf("# -------------------------------------------------------------------------------------\n");

/*
 *	Step 1: loop over rotational states j of the diatom and solve the diatomic eigenvalue problem
 *	using the Fourier grid Hamiltonian discrete variable representation (FGH-DVR) method.
 */

	int max_ch = 0;
	const int max_state = -666;

	for (int j = j_min; j <= j_max; j += j_step)
	{
		double *pot_energy = allocate(grid_size, sizeof(double), false);

		for (int n = 0; n < grid_size; ++n)
		{
			const double r = r_min + as_double(n)*grid_step;
			pot_energy[n] = pec(arrang, r) + centr_term(j, mass(m), r);
		}

		matrix *eigenvec = fgh_matrix(grid_size, grid_step, pot_energy, mass(m));

		free(pot_energy);

		double *eigenval = matrix_symm_eigen(eigenvec, 'v');

/*
 *		Step 2: loop over the vibrational states v of the diatom, solutions of step 2, and select
 *		only those of interest.
 */

		for (int v = v_min; v <= v_max; v += v_step)
		{
			double *wavef = matrix_raw_col(eigenvec, v, false);

			norm(grid_size, grid_step, (grid_size >= 3000), wavef);

/*
 *			Step 3: loop over all partial waves l of the atom around the diatom given by the
 *			respective J and j.
 */

			for (int l = abs(J - j); l <= (J + j); ++l)
			{
				if (parity(j + l) != J_parity && J_parity != 0) continue;

				printf(FORMAT, max_ch, v, j, l, parity(j + l),
				       eigenval[v], eigenval[v]*219474.63137054, eigenval[v]*27.211385);

/*
 *				Step 4: save each basis function |vjl> in the disk and increment the counter of
 *				channels.
 */

				FILE *output = open_basis_file("wb", arrang, max_ch, J);

				fwrite(&v, sizeof(int), 1, output);
				fwrite(&j, sizeof(int), 1, output);
				fwrite(&l, sizeof(int), 1, output);
				fwrite(&max_state, sizeof(int), 1, output);

				fwrite(&r_min, sizeof(double), 1, output);
				fwrite(&r_max, sizeof(double), 1, output);
				fwrite(&grid_step, sizeof(double), 1, output);

				fwrite(&eigenval[v], sizeof(double), 1, output);

				fwrite(&grid_size, sizeof(int), 1, output);
				fwrite(wavef, sizeof(double), grid_size, output);

				fclose(output);
				++max_ch;
			}

			free(wavef);
		}

		matrix_free(eigenvec);
		free(eigenval);
	}

	printf("\n#\n# A total of %d basis functions are computed with %d grid "
	       "points in r = [%f, %f)\n", max_ch, grid_size, r_min, r_max);

	return EXIT_SUCCESS;
}
