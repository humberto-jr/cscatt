#include "modules/pes.h"
#include "modules/file.h"
#include "modules/mpi_lib.h"
#include "modules/globals.h"

#define FORMAT "  %3d       %3d       %3d    %06f         %f\n"

/******************************************************************************

 Function driver(): performs the calculation of all multipole coefficients for
 all lambda- and r-values.

******************************************************************************/

void driver(const char arrang, const bool use_omp, pes_multipole *m)
{
	ASSERT(m != NULL)
	ASSERT(m->value != NULL)

	#pragma omp parallel for default(none) shared(m) schedule(static) if(use_omp)
	for (int lambda = m->lambda_min; lambda <= m->lambda_max; lambda += m->lambda_step)
	{
		const double start_time = wall_time();

		for (int n = 0; n < m->grid_size; ++n)
		{
			const double r = m->r_min + as_double(n)*m->r_step;

			m->value[lambda][n] = pes_legendre_multipole(arrang, lambda, r, m->R);
		}

		const double end_time = wall_time();

		#pragma omp critical
		printf(FORMAT, mpi_rank(), thread_id(), lambda, m->R, end_time - start_time);
	}
}

/******************************************************************************

 Function save_multipole(): writes in the disk the multipole for a given R grid
 index, n, and arrangement. Binary format is used.

******************************************************************************/

void save_multipole(const char arrang, const int n, const pes_multipole *m)
{
	FILE *output = pes_multipole_file(arrang, n, "wb", false);

	pes_multipole_write(m, output);

	fclose(output);
}

/******************************************************************************
******************************************************************************/

int main(int argc, char *argv[])
{
	mpi_init(argc, argv);

	file_init_stdin(argv[1]);

/*
 *	Arrangement (a = 1, b = 2, c = 3), atomic masses and PES:
 */

	const char arrang = 96 + (int) file_keyword(stdin, "arrang", 1.0, 3.0, 1.0);

	pes_init_mass(stdin, 'a');
	pes_init_mass(stdin, 'b');
	pes_init_mass(stdin, 'c');

	pes_init();

/*
 *	Vibrational grid:
 */

	const int rovib_grid_size = (int) file_keyword(stdin, "rovib_grid_size", 1.0, INF, 500.0);

	const double r_min = file_keyword(stdin, "r_min", 0.0, INF, 0.5);

	const double r_max = file_keyword(stdin, "r_max", r_min, INF, r_min + 30.0);

	const double r_step = (r_max - r_min)/as_double(rovib_grid_size);

/*
 *	Scattering grid:
 */

	const int scatt_grid_size = (int) file_keyword(stdin, "scatt_grid_size", 1.0, INF, 500.0);

	const double R_min = file_keyword(stdin, "R_min", 0.0, INF, 0.5);

	const double R_max = file_keyword(stdin, "R_max", R_min, INF, R_min + 30.0);

	const double R_step = (R_max - R_min)/as_double(scatt_grid_size);

/*
 *	Multipoles:
 */

	const int lambda_min = (int) file_keyword(stdin, "lambda_min", 0.0, INF, 0.0);

	const int lambda_max = (int) file_keyword(stdin, "lambda_max", 0.0, INF, 20.0);

	const int lambda_step = (int) file_keyword(stdin, "lambda_step", 1.0, INF, 2.0);

	ASSERT(lambda_max >= lambda_min)

/*
 *	OpenMP:
 */

	const bool use_omp = (bool) file_keyword(stdin, "use_omp", 0.0, 1.0, 1.0);

/*
 *	MPI:
 */

	mpi_set_tasks(scatt_grid_size);

/*
 *	Resolve all tasks:
 */

	if (mpi_rank() == 0)
	{
		printf("# MPI CPUs = %d, OpenMP threads = %d\n", mpi_comm_size(), max_threads());
		printf("# CPU    thread    lambda    R (a.u.)    wall time (s)\n");
		printf("# ----------------------------------------------------\n");
	}

	pes_multipole m =
	{
		.R = 0.0,
		.value = NULL,
		.r_min = r_min,
		.r_max = r_max,
		.r_step = r_step,
		.lambda_min = lambda_min,
		.lambda_max = lambda_max,
		.lambda_step = lambda_step,
		.grid_size = rovib_grid_size
	};

	pes_multipole_init(&m);

	for (int n = mpi_first_task(); n <= mpi_last_task(); ++n)
	{
		extra_step:
		m.R = R_min + as_double(n)*R_step;

		driver(arrang, use_omp, &m);
		save_multipole(arrang, n, &m);

		if (n == mpi_last_task() && mpi_extra_task() > 0)
		{
			n = mpi_extra_task();
			goto extra_step;
		}
	}

	pes_multipole_free(&m);

	mpi_end();
	return EXIT_SUCCESS;
}
