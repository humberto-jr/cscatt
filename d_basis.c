#include "modules/pes.h"
#include "modules/file.h"
#include "modules/matrix.h"
#include "modules/globals.h"

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

	matrix *result = matrix_alloc(grid_size, grid_size, false);

	const double box_length = as_double(grid_size - 1)*grid_step;

	const double factor = (M_PI*M_PI)/(mass*box_length*box_length);

	const double nn_term = factor*as_double(grid_size*grid_size + 2)/6.0;

	for (int n = 0; n < grid_size; ++n)
	{
		matrix_diag_set(result, n, nn_term + pot_energy[n]);

		for (int m = (n + 1); m < grid_size; ++m)
		{
			const double nm = as_double((n + 1) - (m + 1));

			const double nm_term = sin(nm*M_PI/as_double(grid_size));

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

	const int J = (int) file_keyword(stdin, "J", 0.0, INF, 0.0);

	const int J_parity = (int) file_keyword(stdin, "parity", -1.0, 1.0, 0.0);

/*
 *	Vibrational quantum numbers, v:
 */

	const int v_min = (int) file_keyword(stdin, "v_min", 0.0, INF, 0.0);

	const int v_max = (int) file_keyword(stdin, "v_max", 0.0, INF, 0.0);

	const int v_step = (int) file_keyword(stdin, "v_step", 1.0, INF, 1.0);

	ASSERT(v_max >= v_min)

/*
 *	Rotational quantum numbers, j:
 */

	const int j_min = (int) file_keyword(stdin, "j_min", 0.0, INF, 0.0);

	const int j_max = (int) file_keyword(stdin, "j_max", 0.0, INF, 0.0);

	const int j_step = (int) file_keyword(stdin, "j_step", 1.0, INF, 1.0);

	ASSERT(j_max >= j_min)

/*
 *	Vibrational grid:
 */

	const int n_max = (int) file_keyword(stdin, "rovib_grid_size", 1.0, INF, 500.0);

	const double r_min = file_keyword(stdin, "r_min", 0.0, INF, 0.5);

	const double r_max = file_keyword(stdin, "r_max", r_min, INF, r_min + 30.0);

	const double r_step = (r_max - r_min)/as_double(n_max);

	ASSERT(n_max >= v_max + 1)

/*
 *	Arrangement (a = 1, b = 2, c = 3), atomic masses and PES:
 */

	const char arrang = 96 + (int) file_keyword(stdin, "arrang", 1.0, 3.0, 1.0);

	pes_init_mass(stdin, 'a');
	pes_init_mass(stdin, 'b');
	pes_init_mass(stdin, 'c');

	pes_init();

	double (*pec)(const int, const double), mass = 0.0;

	switch (arrang)
	{
		case 'a':
			pec = pes_bc;
			mass = pes_mass_bc();
		break;

		case 'b':
			pec = pes_ac;
			mass = pes_mass_ac();
		break;

		case 'c':
			pec = pes_ab;
			mass = pes_mass_ab();
		break;

		default:
			PRINT_ERROR("invalid arrangement %c\n", arrang)
			exit(EXIT_FAILURE);
	}

	ASSERT(mass != 0.0)

/*
 *	Resolve the diatomic eigenvalue for each j-case and sort results as scatt. channels:
 */

	printf("#\n");
	printf("# J = %d, J-parity = %d, reduced mass = %f a.u.\n", J, J_parity, mass);
	printf("#   Ch.       v       j       l       p        E (a.u.)       E (cm-1)        E (eV)   \n");
	printf("# -------------------------------------------------------------------------------------\n");

/*
 *	Step 1: loop over rotational states j of the diatom and solve the diatomic eigenvalue problem
 *	using the Fourier grid Hamiltonian discrete variable representation (FGH-DVR) method.
 */

	int max_ch = 0;
	const int i = 0;

	for (int j = j_min; j <= j_max; j += j_step)
	{
		double *pot_energy = allocate(n_max, sizeof(double), false);

		for (int n = 0; n < n_max; ++n)
		{
			pot_energy[n] = pec(j, r_min + as_double(n)*r_step);
		}

		matrix *eigenvec = fgh_matrix(n_max, r_step, pot_energy, mass);

		free(pot_energy);

		double *eigenval = matrix_symm_eigen(eigenvec, 'v');

/*
 *		Step 2: loop over the vibrational states v of the diatom, solutions of step 2, and select
 *		only those of interest.
 */

		for (int v = v_min; v <= v_max; v += v_step)
		{
			double *wavef = matrix_raw_col(eigenvec, v, false);

			norm(n_max, r_step, (n_max >= 3000), wavef);

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
				fwrite(&i, sizeof(int), 1, output);

				fwrite(&r_min, sizeof(double), 1, output);
				fwrite(&r_max, sizeof(double), 1, output);
				fwrite(&r_step, sizeof(double), 1, output);

				fwrite(&eigenval[v], sizeof(double), 1, output);

				fwrite(&n_max, sizeof(int), 1, output);
				fwrite(wavef, sizeof(double), n_max, output);

				fclose(output);
				++max_ch;
			}

			free(wavef);
		}

		matrix_free(eigenvec);
		free(eigenval);
	}

	return EXIT_SUCCESS;
}
