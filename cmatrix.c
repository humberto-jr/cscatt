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

/*
 *	Total angular momentum, J:
 */

	const int J
		= (int) file_get_key(stdin, "J", 0.0, INF, 0.0);

/*
 *	Scattering grid, R:
 */

	const int scatt_grid_size
		= (int) file_get_key(stdin, "scatt_grid_size", 1.0, INF, 100.0);

	const double R_min
		= file_get_key(stdin, "R_min", 0.0, INF, 0.5);

	const double R_max
		= file_get_key(stdin, "R_max", R_min, INF, R_min + 100.0);

	const double R_step
		= (R_max - R_min)/as_double(scatt_grid_size);

/*
 * Arrangement (1 == a, 2 == b, 3 == c) and channels:
 */

	const char arrang
		= 96 + (int) file_get_key(stdin, "arrang", 1.0, 3.0, 1.0);

	const int max_ch
		= count_basis(arrang, J);

/*
 *	OpenMP (0 = false, 1 = true):
 */

	const bool use_omp
		= (bool) file_get_key(stdin, "use_omp", 0.0, 1.0, 1.0);

	const int max_thread
		= (use_omp? omp_get_max_threads() : 1);

	const int omp_last_task
		= max_ch*(max_ch + 1)/2;

	const int omp_task_chunk
		= omp_last_task/max_thread;

/*
 *	MPI:
 *
 *	NOTE: On regarding MPI, chunks of R-values are computed independently by each
 *	CPU, where the grid points n = [n_min, n_max] are determined by their IDs. If
 *	the total number of tasks are not equally divided by the total number of CPUs
 *	(common case), a remainder shall be computed after, in a second round of
 *	calculations. If sequential running are used, the chunk of tasks is of
 *	size scatt_grid_size and, thus, remainder = 0, n_min = 0 and n_max =
 *	(scatt_grid_size - 1).
 *
 *	NOTE: For output, when MPI running are considered, each CPU shall open its
 *	own output file. Filenames are built using the CPU ID, the respective total
 *	angular momentum of the problem, J, and the arrangement index; thus, files
 *	belonging to different calculations and processes would not overlap. The
 *	filename template is defined in the macro STDOUT_BUFFER_FORMAT. If a
 *	sequential running is considered, then the output is the usual
 *	C stdout.
 */

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

	#if defined(USE_MPI)
		open_mpi_stdout(cpu_id, arrang, J, false);
	#endif

/*
 *	Atomic masses:
 */

	const enum mass_case m
		= init_atomic_masses(stdin, arrang, 'a');

/*
 * Multipolar coefficients, lambda:
 */

	const int lambda_max
		= (int) file_get_key(stdin, "lambda_max", 0.0, INF, 20.0);

	const int lambda_step
		= init_lambda_step(arrang);

/*
 *	Sum up the main entries before starting actual computations:
 */

	printf("#\n");
	printf("# COUPLING MATRIX FOR J = %d\n", J);
	printf("# Scatt. grid  = %d\n", scatt_grid_size);
	printf("# MPI CPUs     = %d\n", max_cpu);
	printf("# Grid/CPU     = %d\n", mpi_task_chunk);
	printf("# Grid used    = [%d, %d]\n", n_min, n_max);
	printf("# Remainder    = %d\n", remainder);
	printf("# Max. channel = %d\n", max_ch);
	printf("# Tasks/grid   = %d\n", omp_last_task);
	printf("# OMP threads  = %d\n", max_thread);
	printf("# Tasks/thread = %d\n", omp_task_chunk);
	printf("# Max. lambda  = %d\n", lambda_max);
	printf("# Lambda step  = %d\n", lambda_step);

	scatt_basis *basis
		= allocate(max_ch, sizeof(scatt_basis), true);

	for (int n = 0; n < max_ch; ++n)
	{
		load_basis(arrang, n, J, &basis[n]);
	}

/*
 *	NOTE: Once the basis set is built we shall compute the coupling matrix U as
 *	function of the scatt. coordinate R. If a MPI running is considered, then
 *	each CPU will resolve the matrices for given chunks of R values. If OMP
 *	is also used, threads will handle the ro-vibrational integrals of each
 *	matrix element, U(R) = <v'j'l'| V(r, R, theta) |vjl>.
 */

/*
 *	NOTE: As the coupling, U, is a symmetric matrix of shape max_ch-by-max_ch,
 *	the usual way of using two for-loops to iterate over its triangular upper
 *	(or lower) part would give rise a fairly unbalanced workload of tasks for
 *	each thread (when OMP is used). In order to overcome this problem we shall
 *	unfold the upper triangular part into a vector max_ch*(max_ch + 1)/2 long
 *	(including diagonal elements). Thus, the team of threads would divide the
 *	matrix elements into equally spaced chunks of tasks (as much as possible),
 *	being more efficient in terms of computing time.
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

/*
 *	Build the coupling matrix, U(Rn) for n = [n_min, n_max]:
 */

	for (int n = n_min; n <= n_max; ++n)
	{
		extra_step: /* This is an empty statement */;

		matrix *u = matrix_alloc(max_ch, max_ch, false);
		const double R = R_min + as_double(n)*R_step;
		const double start_time = wall_time();

		#pragma omp parallel for default(none) shared(u, list) schedule(static) if(use_omp)
		for (int m = 0; m < omp_last_task; ++m)
		{
			double result
				= integral(J, arrang, lambda_max, lambda_step, list[m].basis_a, list[m].basis_b, R);

			if (list[m].ch_a == list[m].ch_b)
			{
				result += phys_centr_term(list[m].basis_a->l, mass(m), R) + list[m].basis_a->energy;
				matrix_diag_set(u, list[m].ch_a, result);
			}
			else
			{
				matrix_symm_set(u, list[m].ch_a, list[m].ch_b, result);
			}
		}

		const double end_time = wall_time();

		char filename[MAX_LINE_LENGTH];
		sprintf(filename, CMATRIX_BUFFER_FORMAT, arrang, n, J);

		matrix_save(u, filename);
		matrix_free(u);

		printf("#\n");
		printf("# R           = %f\n", R);
		printf("# Wall time   = %f s\n", end_time - start_time);
		printf("# Disk buffer = %s\n", filename);

/*
 *		NOTE: If remainder > 0, each CPU shall compute one extra grid point
 *		after mpi_last_task. Where, CPU 0 handles the (mpi_last_task + 1)
 *		point and, likewise, CPU 1 for the (mpi_last_task + 2), ..., until
 *		CPU N for the (mpi_last_task + N + 1). Thus, they have to redefine
 *		the grid index n and jump back to their new R value in order to
 *		loop over.
 */

		if (n == n_max && remainder > 0)
		{
			n = mpi_last_task + cpu_id + 1;
			if (n < scatt_grid_size) goto extra_step;
		}
	}

	free(list);
	free(basis);

	#if defined(USE_MPI)
		fclose(output);
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
	#endif

	return EXIT_SUCCESS;
}
