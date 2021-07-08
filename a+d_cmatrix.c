#include "modules/pes.h"
#include "modules/fgh.h"
#include "modules/file.h"
#include "modules/math.h"
#include "modules/matrix.h"
#include "modules/mpi_lib.h"
#include "modules/globals.h"

#if !defined(COUPLING_MATRIX_FILE_FORMAT)
	#define COUPLING_MATRIX_FILE_FORMAT "cmatrix_arrang=%c_n=%zu_J=%zu.bin"
#endif

struct tasks
{
	size_t a, b;
	fgh_basis *basis_a, *basis_b;
};

/******************************************************************************

 Function centr_term(): returns the centrifugal potential term at x for a given
 angular momentum l and mass.

******************************************************************************/

inline static double centr_term(const size_t l,
                                const double mass, const double x)
{
	return as_double(l*(l + 1))/(2.0*mass*x*x);
}

/******************************************************************************

 Function simpson(): performs a 3/8-Simpson quadrature rule in the product of
 wavefunctions |a>, |b> and the respective interaction potential V, <a|V|b>.

******************************************************************************/

double simpson(const size_t grid_size,
               const double grid_step,
               const double potential[],
               const double wavef_a[],
               const double wavef_b[])
{
	size_t n_max = grid_size - 1;

	while (n_max%3 != 0)
		--n_max;

	ASSERT(n_max > 10)

	double sum = potential[0]*wavef_a[0]*wavef_b[0]
	           + potential[n_max]*wavef_a[n_max]*wavef_b[n_max];

	for (size_t n = 1; n < (n_max - 3); n += 3)
	{
		sum += 3.0*potential[n]*wavef_a[n]*wavef_b[n];

		sum += 3.0*potential[n + 1]*wavef_a[n + 1]*wavef_b[n + 1];

		sum += 2.0*potential[n + 2]*wavef_a[n + 2]*wavef_b[n + 2];
	}

	sum += 3.0*potential[n_max - 1]*wavef_a[n_max - 1]*wavef_b[n_max - 1];

	return 3.0*grid_step*sum/8.0;
}

/******************************************************************************

 Function integral(): computes the matrix element <a|m|b> for all lambda values
 and a given total angular momentum, J, using a 3/8-Simpson quadrature rule.

******************************************************************************/

double integral(const size_t J,
                const pes_multipole *m, const fgh_basis *a, const fgh_basis *b)
{
	ASSERT(a->r_step == b->r_step)
	ASSERT(b->r_step == m->r_step)

	ASSERT(a->grid_size == b->grid_size)
	ASSERT(b->grid_size == m->grid_size)

	double result = 0.0;
	for (size_t lambda = m->lambda_min; lambda <= m->lambda_max; lambda += m->lambda_step)
	{
		ASSERT(m->value[lambda] != NULL)

		const double f
			= math_percival_seaton(J, a->j, b->j, a->l, b->l, lambda);

		if (f == 0.0) continue;

		const double v
			= simpson(m->grid_size, m->r_step, m->value[lambda], a->eigenvec, b->eigenvec);

		result += v*f;
	}

	return result;
}

/******************************************************************************

 Function driver(): compute one task corresponding to one element of the coupling

******************************************************************************/

double driver(const double mass, const size_t J, const size_t max_task,
              const struct tasks job[], const pes_multipole *m, matrix *c, const bool use_omp)
{
	const double start_time = wall_time();

	#pragma omp parallel for default(none) shared(job, m, c) schedule(static) if(use_omp)
	for (size_t task = 0; task < max_task; ++task)
	{
		double result = integral(J, m, job[task].basis_a, job[task].basis_b);

		if (job[task].a == job[task].b)
		{
			result += job[task].basis_a->eigenval;
			result += centr_term(job[task].basis_a->l, mass, m->R);

			matrix_set_diag(c, job[task].a, result);
		}
		else
		{
			matrix_set_symm(c, job[task].a, job[task].b, result);
		}
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
 *	Arrangement (a = 1, b = 2, c = 3) and atomic masses:
 */

	const char arrang = 96 + read_int_keyword(stdin, "arrang", 1, 3, 1);

	pes_init_mass(stdin, 'a');
	pes_init_mass(stdin, 'b');
	pes_init_mass(stdin, 'c');

	const double mass = pes_mass_abc(arrang);

/*
 *	Total angular momentum, J:
 */

	const size_t J_min = read_int_keyword(stdin, "J_min", 0, 10000, 0);

	const size_t J_max = read_int_keyword(stdin, "J_max", J_min, 10000, J_min);

	const size_t J_step = read_int_keyword(stdin, "J_step", 1, 10000, 1);

/*
 * Directory to read in all basis functions from:
 */

	char *dir = read_str_keyword(stdin, "basis_dir", ".");

/*
 *	MPI: whereas each OpenMP thread handle the integration of each matrix element,
 *	MPI processes are used to handle each scattering grid point, from R_min to R_max.
 */

	const size_t scatt_grid_size = pes_multipole_count(arrang);

	ASSERT(scatt_grid_size > 0)

	mpi_set_tasks(scatt_grid_size);

/*
 *	Load in all J-dependent FGH basis functions and compute the coupling matrix
 *	for each:
 */

	for (size_t J = J_min; J <= J_max; J += J_step)
	{
		const size_t max_channel = fgh_basis_count(dir, arrang, J);

		ASSERT(max_channel > 0)

		fgh_basis *basis = allocate(max_channel, sizeof(fgh_basis), true);

		for (size_t n = 0; n < max_channel; ++n)
			fgh_basis_load(&basis[n], dir, arrang, n, J);

/*
 *		OpenMP: after reading all n basis functions from the disk we sort them in
 *		an ordered list of n*(n + 1)/2 tasks, i.e. the upper triangular part of a
 *		symmetric matrix. Thus, a better workload of tasks per thread is made.
 */

		const bool use_omp = read_int_keyword(stdin, "use_omp", 0, 1, 0);

		const size_t max_task = max_channel*(max_channel + 1)/2;

		struct tasks *list = allocate(max_task, sizeof(struct tasks), true);

		size_t counter = 0;
		for (size_t n = 0; n < max_channel; ++n)
		{
			for (size_t m = n; m < max_channel; ++m)
			{
				list[counter].a = n;
				list[counter].basis_a = &basis[n];

				list[counter].b = m;
				list[counter].basis_b = &basis[m];

				++counter;
			}
		}

		ASSERT(counter == max_task)

/*
 *		Resolve all tasks:
 */

		if (mpi_rank() == 0 && J == J_min)
		{
			printf("# MPI CPUs = %zu, OMP threads = %d, num. of grid points = %zu, mass = %f\n",
			       mpi_comm_size(), max_threads(), scatt_grid_size, mass);

			printf("#  CPU      J     ch.     R (a.u.)      time (s)\n");
			printf("# ----------------------------------------------\n");
		}

		pes_multipole m;
		matrix *c = matrix_alloc(max_channel, max_channel, false);

		for (size_t n = mpi_first_task(); n <= mpi_last_task(); ++n)
		{
			extra_step:
			matrix_set_zero(c);
			pes_multipole_load(&m, arrang, n);

			const double wtime = driver(mass, J, max_task, list, &m, c, use_omp);

			char filename[MAX_LINE_LENGTH];
			sprintf(filename, COUPLING_MATRIX_FILE_FORMAT, arrang, n, J);

			matrix_save(c, filename);

			printf("  %4zu   %4zu   %4zu      %06f      %f\n", mpi_rank(), J, max_channel, m.R, wtime);

			pes_multipole_free(&m);

			if (n == mpi_last_task() && mpi_extra_task() > 0)
			{
				n = mpi_extra_task();
				goto extra_step;
			}
		}

		matrix_free(c);
		free(list);

		for (size_t n = 0; n < max_channel; ++n)
			if (basis[n].eigenvec != NULL) free(basis[n].eigenvec);

		free(basis);
	}

	mpi_end();
	return EXIT_SUCCESS;
}
