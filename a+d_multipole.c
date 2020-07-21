#include "modules/pes.h"
#include "modules/file.h"
#include "modules/mpi_lib.h"
#include "modules/globals.h"

#define FORMAT "  %3d       %3d       %3d    %06f    %06f       % -8e         %f\n"

struct tasks
{
	double r, R, result;
	int lambda, grid_index;
};

/******************************************************************************

 Function send_results(): send all results from a given CPU to that same
 interval in the list of tasks of CPU 0.

******************************************************************************/

void send_results(const int max_task, struct tasks *list)
{
	if (mpi_comm_size() == 1) return;

	ASSERT(max_task > 0)
	ASSERT(list != NULL)

	mpi_barrier();

	if (mpi_rank() == 0)
	{
		for (int rank = 1; rank < mpi_comm_size(); ++rank)
		{
			do
			{
				int n = max_task - 1;
				mpi_receive(rank, 1, type_int, &n);
				mpi_receive(rank, 1, type_double, &list[n].result);
			}
			while (mpi_check(rank));
		}
	}
	else
	{
		for (int n = mpi_first_task(); n <= mpi_last_task(); ++n)
		{
			extra_step:
			mpi_send(0, 1, type_int, &n);
			mpi_send(0, 1, type_double, &list[n].result);

			if (n == mpi_last_task() && mpi_extra_task() > 0)
			{
				n = mpi_extra_task();
				goto extra_step;
			}
		}
	}
}

/******************************************************************************

 Function driver(): performs the calculation of a single task using one thread
 from one process.

******************************************************************************/

void driver(const char arrang, struct tasks *job)
{
	const double start_time = wall_time();

	job->result
		= pes_legendre_multipole(arrang, job->lambda, job->r, job->R);

	const double end_time = wall_time();

	#pragma omp critical
	{
		printf(FORMAT, mpi_rank(), thread_id(),
		       job->lambda, job->r, job->R, job->result, end_time - start_time);
	}
}

/******************************************************************************

 Function multipole_file(): opens the file for a given arrangement and grid
 index. Where, mode is the file access mode of fopen() from the C library.

 NOTE: mode = "wb" for write + binary format and "rb" for read + binary format.

******************************************************************************/

FILE *multipole_file(const char arrang,
                     const int grid_index, const char mode[])
{
	char filename[MAX_LINE_LENGTH];

	sprintf(filename, "multipole_arrang=%c_n=%d.bin", arrang, grid_index);

	return fopen(filename, mode);
}

/******************************************************************************

 Function sort_results(): writes the respective multipoles as function of r per
 lambda using binary format and filename labeled by the grid index, n, of R and
 arrangement.

******************************************************************************/

void sort_results(const char arrang,
                  const int max_task,
                  const int lambda_max,
                  const int grid_size,
                  const double r_min,
                  const double r_max,
                  const double r_step,
                  const struct tasks list[])
{
	FILE *output = NULL;
	int last_index = -1, last_lambda = -1;

	for (int n = 0; n < max_task; ++n)
	{
		if (list[n].grid_index != last_index)
		{
			if (output != NULL) fclose(output);
			output = multipole_file(arrang, list[n].grid_index, "wb");

			ASSERT(output != NULL)

			fwrite(&list[n].R, sizeof(double), 1, output);

			fwrite(&r_min, sizeof(double), 1, output);
			fwrite(&r_max, sizeof(double), 1, output);
			fwrite(&r_step, sizeof(double), 1, output);

			fwrite(&lambda_max, sizeof(int), 1, output);
			fwrite(&grid_size, sizeof(int), 1, output);

			last_lambda = -1;
		}

		if (list[n].lambda != last_lambda)
		{
			fwrite(&list[n].lambda, sizeof(int), 1, output);
		}

		fwrite(&list[n].result, sizeof(double), 1, output);

		last_index = list[n].grid_index;
		last_lambda = list[n].lambda;
	}

	if (output != NULL) fclose(output);
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

	const int rovib_grid_size
		= (int) file_keyword(stdin, "rovib_grid_size", 1.0, INF, 500.0);

	const double r_min
		= file_keyword(stdin, "r_min", 0.0, INF, 0.5);

	const double r_max
		= file_keyword(stdin, "r_max", r_min, INF, r_min + 30.0);

	const double r_step
		= (r_max - r_min)/as_double(rovib_grid_size);

/*
 *	Scattering grid:
 */

	const int scatt_grid_size
		= (int) file_keyword(stdin, "scatt_grid_size", 1.0, INF, 500.0);

	const double R_min
		= file_keyword(stdin, "R_min", 0.0, INF, 0.5);

	const double R_max
		= file_keyword(stdin, "R_max", R_min, INF, R_min + 30.0);

	const double R_step
		= (R_max - R_min)/as_double(scatt_grid_size);

/*
 *	Multipoles:
 */

	const int lambda_min
		= (int) file_keyword(stdin, "lambda_min", 0.0, INF, 0.0);

	const int lambda_max
		= (int) file_keyword(stdin, "lambda_max", as_double(lambda_min), INF, 20.0);

	const int lambda_step
		= (int) file_keyword(stdin, "lambda_step", 1.0, INF, 2.0);

	const int lambda_grid_size
		= (lambda_max - lambda_min)/lambda_step + 1;

/*
 *	OpenMP:
 */

	const bool use_omp = (bool) file_keyword(stdin, "use_omp", 0.0, 1.0, 1.0);

/*
 *	Sort each task (r, R, lambda) in an ordered list:
 */

	const int max_task = scatt_grid_size*rovib_grid_size*lambda_grid_size;

	struct tasks *list = allocate(max_task, sizeof(struct tasks), true);

	int counter = 0;
	for (int n = 0; n < scatt_grid_size; ++n)
	{
		for (int lambda = lambda_min; lambda <= lambda_max; lambda += lambda_step)
		{
			for (int m = 0; m < rovib_grid_size; ++m)
			{
				list[counter].r = r_min + as_double(m)*r_step;
				list[counter].R = R_min + as_double(n)*R_step;
				list[counter].lambda = lambda;
				list[counter].grid_index = n;
				++counter;
			}
		}
	}

	ASSERT(counter == max_task);

/*
 *	MPI:
 */

	mpi_set_tasks(max_task);

/*
 *	Resolve all tasks:
 */

	if (mpi_rank() == 0)
	{
		printf("# MPI CPUs = %d, OpenMP threads = %d, num. of tasks = %d\n", mpi_comm_size(), max_threads(), max_task);
		printf("# CPU    thread    lambda    r (a.u.)    R (a.u.)    multipole (a.u.)    wall time (s)\n");
		printf("# ------------------------------------------------------------------------------------\n");
	}

	#pragma omp parallel for default(none) shared(list) schedule(static) if(use_omp)
	for (int n = mpi_first_task(); n <= mpi_last_task(); ++n)
	{
		extra_step:
		driver(arrang, &list[n]);

		if (n == mpi_last_task() && mpi_extra_task() > 0)
		{
			n = mpi_extra_task();
			goto extra_step;
		}
	}

	send_results(max_task, list);

	if (mpi_rank() == 0)
	{
		sort_results(arrang, max_task, lambda_max,
	                rovib_grid_size, r_min, r_max, r_step, list);
	}

	free(list);

	mpi_end();
	return EXIT_SUCCESS;
}
