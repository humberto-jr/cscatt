#include "modules/pes.h"
#include "modules/file.h"
#include "modules/globals.h"

#if defined(USE_MPI)
	#include <mpi.h>
#endif

#define FORMAT "  %3d       %3d       %3d    %06f    %06f       % -8e         %f\n"

struct tasks
{
	double r, R, result;
	int lambda, grid_index;
};

/******************************************************************************

 Function mpi_print(): writes a short log from a given CPU in the C stdout of
 CPU 0. Notice, if OpenMP is used, only the thread 0 from each CPU sends data
 whereas only the thread 0 from CPU 0 receives.

******************************************************************************/

void mpi_print(const int cpu_id,
               const int max_cpu,
               const int max_task,
               const struct tasks *job,
               const double wall_time)
{
	ASSERT(max_cpu > 0)
	ASSERT(max_task > 0)

	#if defined(USE_MPI)
		static int counter = 0;
	#endif

	if (cpu_id == 0)
	{
		printf(FORMAT, cpu_id, thread_id(),
		       job->lambda, job->r, job->R, job->result, wall_time);

		#if defined(USE_MPI)
		{
			#pragma omp master
			{
				++counter;

				for (int id = 1; id < max_cpu; ++id)
				{
					if (counter == max_task) break;

					char log[MAX_LINE_LENGTH];

					MPI_Recv(&log, MAX_LINE_LENGTH, MPI_BYTE,
					         id, 666, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

					printf("%s", log);

					++counter;
				}
			}
		}
		#endif
	}
	else
	{
		#if defined(USE_MPI)
		{
			#pragma omp master
			{
				char log[MAX_LINE_LENGTH];

				sprintf(log, FORMAT, cpu_id, thread_id(),
				        job->lambda, job->r, job->R, job->result, wall_time);

				MPI_Request info;

				MPI_Isend(&log, MAX_LINE_LENGTH,
				          MPI_BYTE, 0, 666, MPI_COMM_WORLD, &info);
			}
		}
		#endif
	}
}

/******************************************************************************

 Function send_results(): send all results in [n_min, n_max] from a given CPU
 to the same interval in the list of tasks of CPU 0.

******************************************************************************/

void send_results(const int cpu_id,
                  const int max_cpu,
                  const int max_task,
                  const int n_min,
                  const int n_max,
                  struct tasks *list)
{
	ASSERT(cpu_id > -1)
	ASSERT(max_cpu > 0)
	ASSERT(max_task > 0)
	ASSERT(list != NULL)
	ASSERT(n_max >= n_min)

	#if defined(USE_MPI)
	{
		if (cpu_id == 0)
		{
			ASSERT(n_min == 0)

			int counter = n_max;

			for (int id = 1; id < max_cpu; ++id)
			{
				if (counter == max_task) break;

				int m_min = 0, m_max = 0;

				MPI_Recv(&m_min, 1,
				         MPI_INT, id, 666, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

				MPI_Recv(&m_max, 1,
				         MPI_INT, id, 667, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

				ASSERT(m_max >= m_min)

				for (int m = m_min; m <= m_max; ++m)
				{
					MPI_Recv(&list[m].result, 1,
					         MPI_DOUBLE, id, m, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					++counter;
				}
			}
		}
		else
		{
			MPI_Request info;

			MPI_Isend(&n_min, 1,
			          MPI_INT, 0, 666, MPI_COMM_WORLD, &info);

			MPI_Isend(&n_max, 1,
			          MPI_INT, 0, 667, MPI_COMM_WORLD, &info);

			for (int n = n_min; n <= n_max; ++n)
			{
				MPI_Isend(&list[n].result, 1,
				          MPI_DOUBLE, 0, n, MPI_COMM_WORLD, &info);
			}
		}
	}
	#endif

	return;
}

/******************************************************************************

 Function driver(): performs the calculation of a single task and prints a log
 in the C stdout of CPU 0.

******************************************************************************/

void driver(const char arrang,
            const int cpu_id,
            const int max_cpu,
            const int max_task,
            struct tasks *job)
{
	const double start_time = wall_time();

	job->result
		= pes_legendre_multipole(arrang, job->lambda, job->r, job->R);

	const double end_time = wall_time();

	mpi_print(cpu_id, max_cpu, max_task, job, end_time - start_time);
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
 arrangement. The data format used is as follows:

 [R]
 [r_min]
 [r_max]
 [r_step]
 [max. lambda]
 [rovib_grid_size]
 [1st lambda]
 [1st multipole]
 [2nd multipole]
 [3rd multipole]
 (...)
 [2nd lambda]
 [1st multipole]
 [2nd multipole]
 [3rd multipole]
 (...)

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
	FILE *output = multipole_file(arrang, list[0].grid_index, "wb");

	ASSERT(output != NULL)

	fwrite(&list[0].R, sizeof(double), 1, output);
	fwrite(&r_min, sizeof(double), 1, output);
	fwrite(&r_max, sizeof(double), 1, output);
	fwrite(&r_step, sizeof(double), 1, output);
	fwrite(&lambda_max, sizeof(int), 1, output);
	fwrite(&grid_size, sizeof(int), 1, output);

	for (int n = 0; n < max_task; ++n)
	{
		if (n > 0 && list[n].grid_index != list[n - 1].grid_index)
		{
			fclose(output);
			output = multipole_file(arrang, list[n].grid_index, "wb");

			ASSERT(output != NULL)

			fwrite(&list[n].R, sizeof(double), 1, output);
			fwrite(&r_min, sizeof(double), 1, output);
			fwrite(&r_max, sizeof(double), 1, output);
			fwrite(&r_step, sizeof(double), 1, output);
			fwrite(&lambda_max, sizeof(int), 1, output);
			fwrite(&grid_size, sizeof(int), 1, output);
		}

		if (n > 0 && list[n].lambda != list[n - 1].lambda)
		{
			fwrite(&list[n].lambda, sizeof(int), 1, output);
		}

		fwrite(&list[n].result, sizeof(double), 1, output);
	}

	fclose(output);
}

/******************************************************************************

******************************************************************************/

int main(int argc, char *argv[])
{
	ASSERT(argc > 1)

	int id = 0, mpi_size = 1;

	#if defined(USE_MPI)
		int thread_level = 0;
		MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &thread_level);
		MPI_Comm_rank(MPI_COMM_WORLD, &id);
		MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
	#endif

	const int cpu_id = id, max_cpu = mpi_size;

	file_init_stdin(argv[1]);

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
 *	Arrangement (a = 1, b = 2, c = 3), atomic masses and PES:
 */

	const char arrang
		= 96 + (int) file_keyword(stdin, "arrang", 1.0, 3.0, 1.0);

	pes_init_mass(stdin, 'a');
	pes_init_mass(stdin, 'b');
	pes_init_mass(stdin, 'c');
	pes_init();

/*
 *	OpenMP:
 */

	const bool use_omp
		= (bool) file_keyword(stdin, "use_omp", 0.0, 1.0, 1.0);

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

	const int mpi_chunk = max_task/max_cpu;

	const int n_min = cpu_id*mpi_chunk;

	const int n_max = n_min + (mpi_chunk - 1);

	const int mpi_last_task = (max_cpu - 1)*mpi_chunk + (mpi_chunk - 1);

	const int remainder = (max_task - 1) - mpi_last_task;

/*
 *	Resolve all tasks:
 */

	if (cpu_id == 0)
	{
		printf("# MPI CPUs = %d, OpenMP threads = %d, num. of tasks = %d\n", max_cpu, max_threads(), max_task);
		printf("# CPU    thread    lambda    r (a.u.)    R (a.u.)    multipole (a.u.)    wall time (s)\n");
		printf("#-------------------------------------------------------------------------------------\n");
	}

	#pragma omp parallel for default(none) shared(list) schedule(static) if(use_omp)
	for (int n = n_min; n <= n_max; ++n)
	{
		extra_step:
		driver(arrang, cpu_id, max_cpu, max_task, &list[n]);

		if (n == n_max && remainder > 0)
		{
			n = mpi_last_task + cpu_id + 1;
			if (n < max_task) goto extra_step;
		}
	}

	send_results(cpu_id, max_cpu, max_task, n_min, n_max, list);

	if (cpu_id == 0)
	{
		sort_results(arrang, max_task, lambda_max,
	                rovib_grid_size, r_min, r_max, r_step, list);
	}

	free(list);

	#if defined(USE_MPI)
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
	#endif

	return EXIT_SUCCESS;
}
