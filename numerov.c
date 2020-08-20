#include "config.h"

/******************************************************************************

 Function johnson_jcp78_numerov(): use the method of B. R. Johnson, Ref. [1],
 to propagate the ratio matrix of a multichannel wavefunction from the radial
 grid point (n - 1) to n at a given total energy.

 NOTE: the potential matrix is used as workspace and its content is destroyed.

******************************************************************************/

void johnson_jcp78_numerov(const double grid_step, matrix *pot_energy,
                           const double tot_energy, const double mass, matrix *ratio)
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

	matrix *w = pot_energy;

	for (int n = 0; n < matrix_row(pot_energy); ++n)
	{
		double t = factor*(tot_energy - matrix_get(pot_energy, n, n));
		matrix_set(w, n, n, 1.0 - t);

		for (int m = (n + 1); m < matrix_col(pot_energy); ++m)
		{
			t = factor*(0.0 - matrix_get(pot_energy, n, m));
			matrix_set(w, n, m, 0.0 - t);

			t = factor*(0.0 - matrix_get(pot_energy, m, n));
			matrix_set(w, m, n, 0.0 - t);
		}
	}

/*
 *	Solve Eq. (22) and (24) of Ref. [1]:
 */

	matrix_inverse(w);

	for (int n = 0; n < matrix_row(pot_energy); ++n)
	{
		double u = 12.0*matrix_get(w, n, n) - 10.0;
		matrix_set(ratio, n, n, u - matrix_get(ratio, n, n));

		for (int m = (n + 1); m < matrix_col(pot_energy); ++m)
		{
			u = 12.0*matrix_get(w, n, m);
			matrix_set(ratio, n, m, u - matrix_get(ratio, n, m));

			u = 12.0*matrix_get(w, m, n);
			matrix_set(ratio, m, n, u - matrix_get(ratio, m, n));
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
 *	OpenMP:
 */

	const bool use_omp = (bool) file_keyword(stdin, "use_omp", 0.0, 1.0, 1.0);

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

	matrix *workspace = NULL;

	for (int n = 0; n < scatt_grid_size; ++n)
	{
		matrix *pot_energy = coupling_read(arrang, n, J, false);

		if (n == 0)
		{
			const int max_channel = matrix_row(pot_energy);
			workspace = matrix_alloc(max_channel, max_channel, false);
		} 

		const double start_time = wall_time();

		for (int m = mpi_first_task(); m <= mpi_last_task(); ++m)
		{
			extra_step: /* This is an empty statement */;

			const double tot_energy = E_min + as_double(m)*E_step;

			matrix *ratio = rmatrix_load(arrang, n, J, false);

			johnson_jcp78_numerov(R_step,
			                      pot_energy, tot_energy, mass, ratio, workspace);

			rmatrix_save(arrang, n, J, false, ratio);

			if (m == mpi_last_task() && mpi_extra_task() > 0)
			{
				m = mpi_extra_task();
				goto extra_step;
			}
		}

		const double end_time = wall_time();

		matrix_free(pot_energy);
	}

	mpi_end();
	return EXIT_SUCCESS;
}
