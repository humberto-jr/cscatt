#include "modules/phys.h"
#include "modules/mass.h"
#include "modules/file.h"
#include "modules/miller.h"
#include "modules/matrix.h"
#include "modules/globals.h"

#include "mpi_config.h"
#include "mass_config.h"
#include "basis_config.h"
#include "coupl_config.h"

struct task
{
	int ch_a, ch_b;
	scatt_basis *basis_a, *basis_b;
};

typedef struct task task;

double integral(const int J,
                const char arrang,
                const int lambda_max,
                const int lambda_step,
                const scatt_basis *a,
                const scatt_basis *b,
                const double R)
{
	double result = 0.0;
	for (int lambda = 0; lambda <= lambda_max; lambda += lambda_step)
	{
		const double f
		= phys_percival_seaton(a->spin_mult, a->j, b->j, a->l, b->l, lambda, J);

		if (f == 0.0) continue;

		const double v
		= miller_jcp69_vib_integral(lambda, arrang, a->wavef, b->wavef, a->r_min, a->r_max, R);

		result += v*f;
	}

	return result;
}

int main(int argc, char *argv[])
{
	int cpu_id = 0, max_cpu = 1;

	#if defined(USE_MPI)
		int thread_level = 0;
		MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &thread_level);
		MPI_Comm_rank(MPI_COMM_WORLD, &cpu_id);
		MPI_Comm_size(MPI_COMM_WORLD, &max_cpu);
	#endif

	ASSERT(argc > 1)
	file_init_stdin(argv[1]);

	/* total angular momentum, J */

	const int J
		= (int) file_get_key(stdin, "J", 0.0, INF, 0.0);

	/* scattering grid, R */

	const int scatt_grid_size
		= (int) file_get_key(stdin, "scatt_grid_size", 1.0, INF, 500.0);

	const double R_min
		= file_get_key(stdin, "R_min", 0.0, INF, 0.5);

	const double R_max
		= file_get_key(stdin, "R_max", R_min, INF, R_min + 100.0);

	const double R_step
		= (R_max - R_min)/as_double(scatt_grid_size);

	/* arrangement (1 == a, 2 == b, 3 == c) and scatt. channels */

	const char arrang
		= 96 + (int) file_get_key(stdin, "arrang", 1.0, 3.0, 1.0);

	const int max_ch
		= count_basis(arrang, J);

	/* OpenMP (0 = false = not use, 1 = true = use) */

	const bool use_omp
		= (bool) file_get_key(stdin, "use_omp", 0.0, 1.0, 1.0);

	const int max_thread
		= (use_omp? omp_get_max_threads() : 1);

	const int omp_last_task
		= max_ch*(max_ch + 1)/2;

	const int omp_task_chunk
		= omp_last_task/max_thread;

	/* MPI */

	const int mpi_task_chunk
		= scatt_grid_size/max_cpu;

	const int n_min
		= cpu_id*mpi_task_chunk;

	const int n_max
		= n_min + (mpi_task_chunk - 1);

	const int mpi_last_task
		= (max_cpu - 1)*mpi_task_chunk + (mpi_task_chunk - 1);

	const int remainder
		= (scatt_grid_size - 1) - mpi_last_task;

	/* atomic masses and arrang. reduced mass */

	const mass_case a
		= init_atomic_masses(stdin, arrang, 'a', cpu_id);

	/* multipolar expansion coefficients, lambda */

	const int lambda_max
		= (int) file_get_key(stdin, "lambda_max", 0.0, INF, 20.0);

	const int lambda_step
		= init_lambda_step(arrang);

	if (cpu_id == 0)
	{
		printf("#\n");
		printf("# COUPLING MATRIX FOR J = %d\n", J);
		printf("# Scatt. grid  = %d\n", scatt_grid_size);
		printf("# MPI CPUs     = %d\n", max_cpu);
		printf("# Grid/CPU     = %d\n", mpi_task_chunk);
		printf("# Remainder    = %d\n", remainder);
		printf("# Max. channel = %d\n", max_ch);
		printf("# Tasks/grid   = %d\n", omp_last_task);
		printf("# OMP threads  = %d\n", max_thread);
		printf("# Tasks/thread = %d\n", omp_task_chunk);
		printf("# Max. lambda  = %d\n", lambda_max);
		printf("# Lambda step  = %d\n", lambda_step);
	}

	scatt_basis *basis
		= allocate(max_ch, sizeof(scatt_basis), true);

	for (int n = 0; n < max_ch; ++n)
	{
		load_basis(arrang, n, J, &basis[n]);
	}

/*
 *	NOTE: As the coupling is a symmetric matrix of shape max_ch-by-max_ch the
 *	usual way of using two for-loops to iterate over its triangular upper (or
 *	lower) part gives rise to a fairly unbalanced workload of tasks for each
 *	thread. In order to overcome this problem we will unfold the upper part in
 *	a long vector of size omp_last_task = max_ch*(max_ch + 1)/2. Thus, the team
 *	of threads can now divide the matrix elements into equally spaced chunks of
 *	tasks (as much as possible) being more efficient in terms of computing time.
 */

	task *list
		= allocate(omp_last_task, sizeof(task), true);

	int counter = 0;
	for (int n = 0; n < max_ch; ++n)
	{
		for (int m = n; m < max_ch; ++m)
		{
			list[counter].ch_a = n;
			list[counter].basis_a = &basis[n];

			list[counter].ch_b = m;
			list[counter].basis_b = &basis[m];

			++counter;
		}
	}

	ASSERT(counter == omp_last_task);
	counter = 0;

	for (int n = n_min; n <= n_max; ++n)
	{
		extra_step: /* This is an empty statement */;

		matrix *c = matrix_alloc(max_ch, max_ch, false);

		const double R = R_min + as_double(n)*R_step;
		const double start_time = wall_time();

		#pragma omp parallel for default(none) shared(c, list) schedule(dynamic, omp_task_chunk) if(use_omp)
		for (int m = 0; m < omp_last_task; ++m)
		{
			double result = integral(J, arrang, lambda_max, lambda_step,
			                         list[m].basis_a, list[m].basis_b, R);

			if (list[m].ch_a == list[m].ch_b)
			{
				result += list[m].basis_a->energy
				        + phys_centr_term(list[m].basis_a->l, mass(a), R);

				matrix_diag_set(c, list[m].ch_a, result);
			}
			else
			{
				matrix_symm_set(c, list[m].ch_a, list[m].ch_b, result);
			}
		}

		const double end_time = wall_time();

		char filename[MAX_LINE_LENGTH];
		sprintf(filename, CMATRIX_BUFFER_FORMAT, arrang, n, J);

		matrix_save(c, filename);
		matrix_free(c);

		if (cpu_id == 0)
		{
			printf("#\n");
			printf("# CPU         = %d\n", cpu_id);
			printf("# R value     = %f a.u.\n", R);
			printf("# Wall time   = %f min\n", (end_time - start_time)/60.0);
			printf("# Disk buffer = %s\n", filename);
			++counter;

			#if defined(USE_MPI)
				for (int their_id = 1; their_id < max_cpu; ++their_id)
				{
					if (counter == scatt_grid_size) break;

					char their_log[MAX_LINE_LENGTH];

					MPI_Recv(&their_log, MAX_LINE_LENGTH, MPI_BYTE,
					         their_id, 666, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

					printf(their_log);
					++counter;
				}
			#endif
		}
		else
		{
			#if defined(USE_MPI)
				MPI_Request info;
				char log[MAX_LINE_LENGTH];

				sprintf(log, "#\n");
				sprintf(log + strlen(log), "# CPU         = %d\n", cpu_id);
				sprintf(log + strlen(log), "# R value     = %f a.u.\n", R);
				sprintf(log + strlen(log), "# Wall time   = %f min\n", (end_time - start_time)/60.0);
				sprintf(log + strlen(log), "# Disk buffer = %s\n", filename);

				MPI_Isend(&log, MAX_LINE_LENGTH,
				          MPI_BYTE, 0, 666, MPI_COMM_WORLD, &info);
			#endif
		}

/*
 *		NOTE: If remainder > 0, each CPU shall compute one extra grid point after
 *		mpi_last_task. Where, CPU 0 handles the (mpi_last_task + 1)-th point and,
 *		likewise, CPU 1 goes for the (mpi_last_task + 2)-th, ..., until CPU N for
 *		the (mpi_last_task + N + 1). Thus, they have to redefine the grid index n
 *		and jump back to their new R value in order to loop over.
 */

		if (n == n_max && remainder > 0)
		{
			n = mpi_last_task + cpu_id + 1;
			if (n < scatt_grid_size) goto extra_step;
		}
	}

	for (int n = 0; n < max_ch; ++n)
	{
		matrix_free(basis[n].wavef);
	}

	free(basis);
	free(list);

	#if defined(USE_MPI)
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
	#endif

	return EXIT_SUCCESS;
}
