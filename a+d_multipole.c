#include "modules/pes.h"
#include "modules/file.h"
#include "modules/mpi_lib.h"
#include "modules/globals.h"

struct tasks
{
	double r;
	size_t lambda, index;
};

/******************************************************************************

 Function driver(): performs the calculation of multipole coefficients, m, for
 all max_task (lambda, r)-values stored in the array job at a fixed R. The
 wall time in seconds is returned.

******************************************************************************/

double driver(const char arrang, const size_t max_task,
              const struct tasks job[], pes_multipole *m, const bool use_omp)
{
	ASSERT(m != NULL)
	ASSERT(job != NULL)
	ASSERT(m->value != NULL)

	const double start_time = wall_time();

	#pragma omp parallel for default(none) shared(job, m) schedule(static) if(use_omp)
	for (size_t task = 0; task < max_task; ++task)
	{
		const size_t n = job[task].index;
		const size_t lambda = job[task].lambda;
		const double r = job[task].r, R = m->R;

		m->value[lambda][n] = pes_legendre_multipole(arrang, lambda, r, R);
	}

	const double end_time = wall_time();

	return (end_time - start_time);
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

	const char arrang = 96 + read_int_keyword(stdin, "arrang", 1, 3, 1);

	pes_init_mass(stdin, 'a');
	pes_init_mass(stdin, 'b');
	pes_init_mass(stdin, 'c');

	pes_init();

/*
 *	Vibrational grid:
 */

	const size_t rovib_grid_size = read_int_keyword(stdin, "rovib_grid_size", 1, 1000000, 1000);

	const double r_min = read_dbl_keyword(stdin, "r_min", 0.0, INF, 0.5);

	const double r_max = read_dbl_keyword(stdin, "r_max", r_min, INF, r_min + 30.0);

	const double r_step = (r_max - r_min)/as_double(rovib_grid_size);

/*
 *	Scattering grid:
 */

	const size_t scatt_grid_size = read_int_keyword(stdin, "scatt_grid_size", 1, 1000000, 500);

	const double R_min = read_dbl_keyword(stdin, "R_min", 0.0, INF, 0.5);

	const double R_max = read_dbl_keyword(stdin, "R_max", R_min, INF, R_min + 30.0);

	const double R_step = (R_max - R_min)/as_double(scatt_grid_size);

/*
 *	Multipoles:
 */

	const size_t lambda_min = read_int_keyword(stdin, "lambda_min", 0, 10000, 0);

	const size_t lambda_max = read_int_keyword(stdin, "lambda_max", 0, 10000, 20);

	const size_t lambda_step = read_int_keyword(stdin, "lambda_step", 1, 10000, 2);

	ASSERT(lambda_max >= lambda_min)

/*
 *	OpenMP: for a given R-value, each (lambda, r)-dependent multipole coefficient
 *	is computed by an OpenMP thread, if any. In order to impose a good workload
 *	of tasks per thread with minimal sync, an array of struct tasks is utilized
 *	to store precomputed values of (lambda, r). Thus, a single for-loop is used.
 */

	const bool use_omp = read_int_keyword(stdin, "use_omp", 0, 1, 0);

	struct tasks *list = NULL;

	size_t counter = 0;
	for (size_t lambda = lambda_min; lambda <= lambda_max; lambda += lambda_step)
	{
		for (size_t n = 0; n < rovib_grid_size; ++n)
		{
			++counter;
			list = realloc(list, sizeof(struct tasks)*counter);

			ASSERT(list != NULL)

			list[counter - 1].r = r_min + as_double(n)*r_step;
			list[counter - 1].lambda = lambda;
			list[counter - 1].index = n;
		}
	}

/*
 *	MPI: the main loop over chunks of R-values is handled by MPI processes, if any.
 */

	mpi_set_tasks(scatt_grid_size);

/*
 *	Directory to store all multipole coefficients:
 */

	char *dir = read_str_keyword(stdin, "multipole_dir", ".");

/*
 *	Resolve all tasks:
 */

	if (mpi_rank() == 0)
	{
		printf("# MPI CPUs = %zu, OpenMP threads = %d, num. of tasks = %zu, PES name = %s\n", mpi_comm_size(), max_threads(), counter*scatt_grid_size, pes_name());
		printf("#  CPU       R (a.u.)    wall time (s)\n");
		printf("# ------------------------------------\n");
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

	for (size_t n = mpi_first_task(); n <= mpi_last_task(); ++n)
	{
		extra_step:
		m.R = R_min + as_double(n)*R_step;

		const double wtime = driver(arrang, counter, list, &m, use_omp);

		pes_multipole_save(&m, dir, arrang, n);

		printf("  %4zu       %06f         %f\n", mpi_rank(), m.R, wtime);

		if (n == mpi_last_task() && mpi_extra_task() > 0)
		{
			n = mpi_extra_task();
			goto extra_step;
		}
	}

	pes_multipole_free(&m);
	free(list);
	free(dir);

	mpi_end();
	return EXIT_SUCCESS;
}
