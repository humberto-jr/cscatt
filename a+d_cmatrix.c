#include "modules/pes.h"
#include "modules/fgh.h"
#include "modules/file.h"
#include "modules/math.h"
#include "modules/matrix.h"
#include "modules/mpi_lib.h"
#include "modules/globals.h"

#define FORMAT "  %3d       %3d       %5d       %5d       %06f    %f\n"

#if !defined(COUPLING_MATRIX_FILE_FORMAT)
	#define COUPLING_MATRIX_FILE_FORMAT "cmatrix_arrang=%c_n=%d_J=%d.bin"
#endif

struct tasks
{
	int a, b;
	fgh_basis *basis_a, *basis_b;
};

typedef struct tasks tasks;

/******************************************************************************

 Function centr_term(): returns the centrifugal potential term at x for a given
 angular momentum l and mass.

******************************************************************************/

inline static double centr_term(const int l, const double mass, const double x)
{
	return as_double(l*(l + 1))/(2.0*mass*x*x);
}

/******************************************************************************

 Function simpson(): performs a 3/8-Simpson quadrature rule in the product of
 wavefunctions |a>, |b> and the respective interaction potential V, <a|V|b>.

******************************************************************************/

double simpson(const int grid_size,
               const double grid_step,
               const double potential[],
               const double wavef_a[],
               const double wavef_b[])
{
	int n_max = grid_size - 1;

	while (n_max%3 != 0) --n_max;

	ASSERT(n_max > 3)

	double sum = potential[0]*wavef_a[0]*wavef_b[0]
	           + potential[n_max]*wavef_a[n_max]*wavef_b[n_max];

	for (int n = 1; n < (n_max - 3); n += 3)
	{
		sum += 3.0*potential[n]*wavef_a[n]*wavef_b[n];

		sum += 3.0*potential[n + 1]*wavef_a[n + 1]*wavef_b[n + 1];

		sum += 2.0*potential[n + 2]*wavef_a[n + 2]*wavef_b[n + 2];
	}

	sum += 3.0*potential[n_max - 1]*wavef_a[n_max - 1]*wavef_b[n_max - 1];

	return 3.0*grid_step*sum/8.0;
}

/******************************************************************************

 Function ab_integral(): performs the integral on the product of wavefunctions
 |a>, |b> and the respective interaction potential V, i.e. <a|V|b>, for all
 lambda values upon which the potential depends.

******************************************************************************/

double ab_integral(const int J,
                   const pes_multipole *m, const fgh_basis *a, const fgh_basis *b)
{
	ASSERT(a->r_step == b->r_step)
	ASSERT(b->r_step == m->r_step)

	ASSERT(a->grid_size == b->grid_size)
	ASSERT(b->grid_size == m->grid_size)

	double result = 0.0;
	for (int lambda = m->lambda_min; lambda <= m->lambda_max; lambda += m->lambda_step)
	{
		ASSERT(m->value != NULL)

		const double f = math_percival_seaton(J, a->j, b->j, a->l, b->l, lambda);

		if (f == 0.0) continue;

		const double v = simpson(m->grid_size, m->r_step,
		                         m->value[lambda], a->eigenvec, b->eigenvec);
		result += v*f;
	}

	return result;
}

/******************************************************************************

 Function driver(): compute one task corresponding to one element of the coupling

******************************************************************************/

void driver(const int J,
            const double mass, const pes_multipole *m, const tasks *job, matrix *c)
{
	const double start_time = wall_time();

	double result = ab_integral(J, m, job->basis_a, job->basis_b);

	const double end_time = wall_time();

	if (job->a == job->b)
	{
		result += job->basis_a->eigenval + centr_term(job->basis_a->l, mass, m->R);
		matrix_diag_set(c, job->a, result);
	}
	else
	{
		matrix_symm_set(c, job->a, job->b, result);
	}

	#pragma omp critical
	printf(FORMAT, mpi_rank(), thread_id(), job->a, job->b, m->R, end_time - start_time);
}

/******************************************************************************

 Function save_coupling(): saves in the disk a coupling matrix, c, for a given
 grid point index, n, arrangement and total angular momentum J.

******************************************************************************/

void save_coupling(const char arrang, const int n, const int J, const matrix *c)
{
	char filename[MAX_LINE_LENGTH];

	sprintf(filename, COUPLING_MATRIX_FILE_FORMAT, arrang, n, J);

	matrix_save(c, filename);
}

/******************************************************************************
******************************************************************************/

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
 *	Scattering grid:
 */

	const int grid_size = pes_multipole_count(arrang);

/*
 *	OpenMP. Read all n basis functions from the disk and sort them in an ordered
 *	list of n*(n + 1)/2 tasks, i.e. the upper triangular part of the symmetric
 *	coupling matrix. Thus, a better workload of tasks per thread, if any.
 */

	const bool use_omp = (bool) file_keyword(stdin, "use_omp", 0.0, 1.0, 1.0);

	const int max_channel = fgh_basis_count(arrang, J);

	const int max_task = max_channel*(max_channel + 1)/2;

	fgh_basis *basis = allocate(max_channel, sizeof(fgh_basis), true);

	for (int n = 0; n < max_channel; ++n)
	{
		FILE *input = fgh_basis_file(arrang, n, J, "rb", false);

		fgh_basis_read(&basis[n], input);

		fclose(input);
	}

	tasks *list = allocate(max_task, sizeof(tasks), true);

	int counter = 0;
	for (int n = 0; n < max_channel; ++n)
	{
		for (int m = n; m < max_channel; ++m)
		{
			list[counter].a = n;
			list[counter].basis_a = &basis[n];

			list[counter].b = m;
			list[counter].basis_b = &basis[m];

			++counter;
		}
	}

	ASSERT(counter == max_task);

/*
 *	MPI. Whereas each OpenMP thread handle the integration of each matrix element,
 *	MPI processes are used to handle each scattering grid point, from R_min to
 *	R_max.
 */

	mpi_set_tasks(grid_size);

/*
 *	Resolve all tasks:
 */

	if (mpi_rank() == 0)
	{
		printf("# MPI CPUs = %d", mpi_comm_size());
		printf(", OMP threads = %d", max_threads());
		printf(", num. of channels = %d", max_channel);
		printf(", num. of grid points = %d\n", grid_size);

		printf("# CPU    thread       Ch. a       Ch. b       R (a.u.)    wall time (s)\n");
		printf("# ---------------------------------------------------------------------\n");
	}

	for (int n = mpi_first_task(); n <= mpi_last_task(); ++n)
	{
		extra_step: /* This is an empty statement */;

		FILE *input = pes_multipole_file(arrang, n, "rb", false);

		pes_multipole m = {.value = NULL};
		pes_multipole_read(&m, input);

		fclose(input);

		matrix *c = matrix_alloc(max_channel, max_channel, true);

		#pragma omp parallel for default(none) shared(list, m, c) schedule(static) if(use_omp)
		for (int task = 0; task < max_task; ++task)
			driver(J, mass, &m, &list[task], c);

		save_coupling(arrang, n, J, c);

		matrix_free(c);
		pes_multipole_free(&m);

		if (n == mpi_last_task() && mpi_extra_task() > 0)
		{
			n = mpi_extra_task();
			goto extra_step;
		}
	}

	for (int n = 0; n < max_channel; ++n)
		if (basis[n].eigenvec != NULL) free(basis[n].eigenvec);

	free(basis);
	free(list);

	mpi_end();
	return EXIT_SUCCESS;
}
