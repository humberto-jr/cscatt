#include "modules/pes.h"
#include "modules/mass.h"
#include "modules/file.h"
#include "modules/globals.h"

#include "mass_config.h"

#if defined(USE_MPI)
	#include <mpi.h>
#endif

#define FORMAT "%d    %d    %d    %f    %f\n"

static int cpu_id, max_cpu, max_task;

struct task
{
	char arrang;
	double r, R, result;
	int lambda, grid_index;
};

void mpi_print(const int lambda, const double R, const double wall_time)
{
	static int counter = 0;

	if (cpu_id == 0)
	{
		printf(FORMAT, cpu_id, thread_id(), lambda, R, wall_time);
		++counter;

		#if defined(USE_MPI)
		{
			for (int id = 1; id < max_cpu; ++id)
			{
				if (counter == max_task) break;

				char log[MAX_LINE_LENGTH];

				MPI_Recv(&log, MAX_LINE_LENGTH, MPI_BYTE,
				         id, 666, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

				printf(log);
				++counter;
			}
		}
		#endif
	}
	else
	{
		#if defined(USE_MPI)
		{
			char log[MAX_LINE_LENGTH];
			sprintf(log, FORMAT, cpu_id, thread_id(), lambda, R, wall_time);

			MPI_Request info;
			MPI_Isend(&log, MAX_LINE_LENGTH, MPI_BYTE, 0, 666, MPI_COMM_WORLD, &info);
		}
		#endif
	}
}

void driver(struct task *job)
{
	const double start_time = wall_time();

	job->result
		= pes_legendre_multipole(job->arrang, job->lambda, job->r, job->R);

	const double end_time = wall_time();

	mpi_print(job->lambda, job->R, end_time - start_time);
}

FILE *multipole_file(const char arrang, const int grid_index, const char mode[])
{
	char filename[MAX_LINE_LENGTH];
	sprintf(filename, "multipole_arrang=%c_n=%d.bin", arrang, grid_index);

	return fopen(filename, mode);
}

void sort_results(const int grid_size,
                  const double r_min,
                  const double r_max,
                  const double r_step,
                  const struct task list[])
{
	FILE *output = NULL;
	int prev_index = -1, prev_lambda = -1;

	for (int n = 0; n < max_task; ++n)
	{
		if (list[n].grid_index != prev_index)
		{
			if (output != NULL) fclose(output);

			output = multipole_file(list[n].arrang, list[n].grid_index, "wb");

			fwrite(&list[n].R, sizeof(double), 1, output);
			fwrite(&r_min, sizeof(double), 1, output);
			fwrite(&r_max, sizeof(double), 1, output);
			fwrite(&r_step, sizeof(double), 1, output);
			fwrite(&grid_size, sizeof(int), 1, output);

			prev_index = list[n].grid_index;
			prev_lambda = -1;
		}

		if (list[n].lambda != prev_lambda)
		{
			fwrite(&list[n].lambda, sizeof(int), 1, output);
			prev_lambda = list[n].lambda;
		}

		if (output != NULL)
		{
			fwrite(&list[n].result, sizeof(double), 1, output);
		}
	}
}

int main(int argc, char *argv[])
{
	cpu_id = 0;
	max_cpu = 1;

	#if defined(USE_MPI)
		int thread_level = 0;
		MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &thread_level);
		MPI_Comm_rank(MPI_COMM_WORLD, &cpu_id);
		MPI_Comm_size(MPI_COMM_WORLD, &max_cpu);
	#endif

	ASSERT(argc > 1)
	file_init_stdin(argv[1]);

	pes_init();

	const int rovib_grid_size
		= (int) file_get_key(stdin, "rovib_grid_size", 1.0, INF, 500.0);

	const double r_min
		= file_get_key(stdin, "r_min", 0.0, INF, 0.5);

	const double r_max
		= file_get_key(stdin, "r_max", r_min, INF, r_min + 30.0);

	const double r_step
		= (r_max - r_min)/as_double(rovib_grid_size);

	const int scatt_grid_size
		= (int) file_get_key(stdin, "scatt_grid_size", 1.0, INF, 500.0);

	const double R_min
		= file_get_key(stdin, "R_min", 0.0, INF, 0.5);

	const double R_max
		= file_get_key(stdin, "R_max", R_min, INF, R_min + 30.0);

	const double R_step
		= (R_max - R_min)/as_double(scatt_grid_size);

	const int lambda_min
		= (int) file_get_key(stdin, "lambda_min", 0.0, INF, 0.0);

	const int lambda_max
		= (int) file_get_key(stdin, "lambda_max", as_double(lambda_min), INF, 20.0);

	const int lambda_step
		= (int) file_get_key(stdin, "lambda_step", 1.0, INF, 2.0);

	const char arrang
		= 96 + (int) file_get_key(stdin, "arrang", 1.0, 3.0, 1.0);

	const bool use_omp
		= (bool) file_get_key(stdin, "use_omp", 0.0, 1.0, 1.0);

	const mass_case a
		= init_atomic_masses(stdin, arrang, 'a', cpu_id);

	max_task
		= scatt_grid_size*rovib_grid_size*((lambda_max - lambda_min)/lambda_step + 1);

	struct task *list = allocate(max_task, sizeof(struct task), true);

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
				list[counter].arrang = arrang;
				list[counter].grid_index = n;
				++counter;
			}
		}
	}

	ASSERT(counter == max_task);

	const int mpi_chunk = max_task/max_cpu;

	const int n_min = cpu_id*mpi_chunk;

	const int n_max = n_min + (mpi_chunk - 1);

	const int mpi_last_task = (max_cpu - 1)*mpi_chunk + (mpi_chunk - 1);

	const int remainder = (max_task - 1) - mpi_last_task;

	const int omp_chunk = max_threads()/(n_max - n_min + 1);

	if (cpu_id == 0)
	{
		printf("\n");
		printf("# CPU    thread    lambda    R (a.u.)    walltime (min)\n");
		printf("#------------------------------------------------------\n");
	}

	#pragma omp parallel for default(none) shared(list) schedule(dynamic, omp_chunk) if(use_omp)
	for (int n = n_min; n <= n_max; ++n)
	{
		extra_step:
		driver(&list[n]);

		if (n == n_max && remainder > 0)
		{
			n = mpi_last_task + cpu_id + 1;
			if (n < max_task) goto extra_step;
		}
	}

	sort_results(rovib_grid_size, r_min, r_max, r_step, list);
	free(list);

	#if defined(USE_MPI)
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
	#endif

	return EXIT_SUCCESS;
}
