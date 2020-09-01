#include "modules/pes.h"
#include "modules/file.h"
#include "modules/matrix.h"
#include "modules/mpi_lib.h"
#include "modules/globals.h"

#include "utils.h"

#define FORMAT "  %3d       %06f    %f\n"

/******************************************************************************

 Function johnson_jcp78_numerov(): uses the method of B. R. Johnson, Ref. [1],
 to propagate the ratio matrix of a multichannel wavefunction from the radial
 grid point (n - 1) to n at a given total energy.

 NOTE: the potential matrix is used as workspace and its content is destroyed
 if workspace = NULL.

******************************************************************************/

void johnson_jcp78_numerov(const double grid_step,
                           const double mass, const double tot_energy,
                           matrix *pot_energy, matrix *ratio, matrix *workspace)
{
	ASSERT(ratio != NULL)
	ASSERT(pot_energy != NULL)

	if (!matrix_null(ratio)) matrix_inverse(ratio);

/*
 *	NOTE: From Eq. (2) and (17) of Ref. [1] the following numerical
 *	factor is defined in atomic units:
 */

	const double factor = -grid_step*grid_step*2.0*mass/12.0;

/*
 *	Resolve Eq. (23) of Ref. [1] with Eq. (2) and (17) plugged in:
 */

	matrix *w = NULL;

	if (workspace == NULL)
		w = pot_energy;
	else
		w = workspace;

	const int n_max = matrix_row(pot_energy);

	for (int n = 0; n < n_max; ++n)
	{
		const double t_nn = factor*(tot_energy - matrix_get(pot_energy, n, n));

		matrix_diag_set(w, n, 1.0 - t_nn);

		for (int m = (n + 1); m < n_max; ++m)
		{
			const double t_nm = factor*(0.0 - matrix_get(pot_energy, n, m));

			matrix_symm_set(w, n, m, 0.0 - t_nm);
		}
	}

/*
 *	Solve Eq. (22) and (24) of Ref. [1]:
 */

	matrix_inverse(w);

	for (int n = 0; n < n_max; ++n)
	{
		const double u_nn = 12.0*matrix_get(w, n, n) - 10.0;

		matrix_diag_set(ratio, n, u_nn - matrix_get(ratio, n, n));

		for (int m = (n + 1); m < n_max; ++m)
		{
			const double u_nm = 12.0*matrix_get(w, n, m);

			matrix_symm_set(ratio, n, m, u_nm - matrix_get(ratio, n, m));
		}
	}

	w = NULL;
}

int main(int argc, char *argv[])
{
	mpi_init(argc, argv);

	file_init_stdin(argv[1]);

/*
 *	Arrangement (a = 1, b = 2, c = 3) and atomic masses:
 */

	const char arrang = 96 + (int) file_keyword(stdin, "arrang", 1.0, 3.0, 1.0);

	pes_init_mass(stdin, 'a');
	pes_init_mass(stdin, 'b');
	pes_init_mass(stdin, 'c');

	const double mass = pes_mass_abc(arrang);

/*
 *	Total angular momentum, J:
 */

	const int J = (int) file_keyword(stdin, "J", 0.0, INF, 0.0);

/*
 *	Collision energy grid:
 */

	const int coll_grid_size = (int) file_keyword(stdin, "coll_grid_size", 1.0, INF, 100.0);

	const double E_min = file_keyword(stdin, "E_min", -INF, INF, 0.0);

	const double E_max = file_keyword(stdin, "E_max", E_min, INF, E_min);

	const double E_step = (E_max - E_min)/as_double(coll_grid_size);

/*
 *	Scattering grid:
 */

	const int scatt_grid_size = (int) file_keyword(stdin, "scatt_grid_size", 1.0, INF, 500.0);

	const double R_min = file_keyword(stdin, "R_min", 0.0, INF, 0.5);

	const double R_max = file_keyword(stdin, "R_max", R_min, INF, R_min + 30.0);

	const double R_step = (R_max - R_min)/as_double(scatt_grid_size);

/*
 *	MPI:
 */

	mpi_set_tasks(coll_grid_size);

/*
 *	Resolve all tasks:
 */

	if (mpi_rank() == 0)
	{
		printf("# MPI CPUs = %d, OpenMP threads = %d, num. of energies = %d\n", mpi_comm_size(), max_threads(), coll_grid_size);
		printf("# CPU    R (a.u.)    wall time (s)\n");
		printf("# --------------------------------\n");
	}

	matrix *pot_energy = NULL, *ratio = NULL, *workspace = NULL;

	for (int n = 0; n < scatt_grid_size; ++n)
	{
		pot_energy = coupling_read(arrang, n, J, false);

		if (n == 0) workspace = matrix_alloc_as(pot_energy, false);

		const double start_time = wall_time();

		for (int m = mpi_first_task(); m <= mpi_last_task(); ++m)
		{
			extra_step:
			ratio = ratio_read(arrang, m, J, false);

			const double tot_energy = E_min + as_double(m)*E_step;

			johnson_jcp78_numerov(R_step, mass,
			                      tot_energy, pot_energy, ratio, workspace);

			ratio_write(arrang, m, J, false, ratio);
			matrix_free(ratio);

			if (m == mpi_last_task() && mpi_extra_task() > 0)
			{
				m = mpi_extra_task();
				goto extra_step;
			}
		}

		const double end_time = wall_time(), R = R_min + as_double(n)*R_step;

		if (mpi_rank() == 0) printf(FORMAT, mpi_rank(), R, end_time - start_time);

		matrix_free(pot_energy);
	}

	matrix_free(workspace);

	mpi_end();
	return EXIT_SUCCESS;
}
