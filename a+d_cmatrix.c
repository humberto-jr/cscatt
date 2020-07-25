#include "modules/pes.h"
#include "modules/file.h"
#include "modules/matrix.h"
#include "modules/mpi_lib.h"
#include "modules/globals.h"

#include "utils.h"

#define FORMAT "  %3d       %3d       %3d    %06f    %06f       % -8e         %f\n"

struct tasks
{
	int a, b;
	basis *basis_a, *basis_b;
};

typedef struct tasks tasks;

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

 Function percival_seaton(): return the so-called Percival & Seaton term often
 used in the definition of the atom-diatom collisional coupling matrix. Where,

 j1     = diatomic rotational angular momentum quantum number for channel 1
 j2     = diatomic rotational angular momentum quantum number for channel 2
 l1     = atom-diatom orbital angular momentum quantum number for channel 1
 l2     = atom-diatom orbital angular momentum quantum number for channel 2
 J      = total angular momentum of the problem
 lambda = an integer parameter (0, 1, 2, ...)

 and also depends parametrically on the electronic spin multiplicity (1 for
 singlet, 2 for doublet and 3 for triplet).

 See Ref. [1] for more, in particular, Eq. (A1).

******************************************************************************/

double percival_seaton(const int j1, const int j2,
                       const int l1, const int l2, const int lambda, const int J)
{
	double result = pow(-1.0, j1 + j2 - J);

	result *= math_wigner_3j(l1, l2, lambda, 0, 0, 0);
	result *= math_wigner_3j(j1, j2, lambda, 0, 0, 0);
	result *= math_wigner_6j(j1, j2, lambda, l2, l1, J);
	result *= sqrt(as_double(2*j1 + 1)*as_double(2*j2 + 1)*as_double(2*l1 + 1)*as_double(2*l2 + 1));

	return result;
}

/******************************************************************************

 Function driver(): performs the calculation of a single task using one thread
 from one process.

******************************************************************************/

double simpson(const int grid_size,
               const double grid_step,
               const double potential[],
               const double wavef_b[],
               const double wavef_a[])
{
	const int n_max = ((grid_size - 2)%3 == 0? grid_size : grid_size - 1); // TODO

	double sum = potential[0]*wavef_a[0]*wavef_b[0]
	           + potential[n_max - 1]*wavef_a[n_max - 1]*wavef_b[n_max - 1];

	for (int n = 1; n < n_max; n += 3)
	{
		sum += 3.0*potential[n]*wavef_a[n]*wavef_b[n];

		sum += 3.0*potential[n + 1]*wavef_a[n + 1]*wavef_b[n + 1];

		sum += 2.0*potential[n + 2]*wavef_a[n + 2]*wavef_b[n + 2];
	}

	return grid_step*3.0*sum/8.0;
}

/******************************************************************************

 Function driver(): performs the calculation of a single task using one thread
 from one process.

******************************************************************************/

double ab_integral(const int J,
                   const multipole_set *m, const basis *a, const basis *b)
{
	ASSERT(m != NULL)
	ASSERT(a != NULL)
	ASSERT(b != NULL)

	double result = 0.0;
	for (int lambda = 0; lambda <= m->lambda_max; ++lambda)
	{
		if (m->set[lambda].value == NULL) continue;

		const double f = percival_seaton(a->j, b->j, a->l, b->l, lambda, J);

		if (f == 0.0) continue;

		const double v = simpson(m->grid_size, m->r_step,
		                         m->set[lambda].value, a->eigenvec, b->eigenvec);
		result += v*f;
	}

	return result;
}

/******************************************************************************

 Function driver(): performs the calculation of a single task using one thread
 from one process.

******************************************************************************/

void solve(const char arrang, const int J,
           const multipole_set *m, const struct tasks *t, matrix *c)
{
	const double start_time = wall_time();

	double result = 0.0;
	for (int lambda = 0; lambda <= m->lambda_max; ++lambda)
	{
		if (m->set[lambda].value == NULL) continue;

		const double f = percival_seaton(a->j, b->j, a->l, b->l, lambda, J);

		if (f == 0.0) continue;

		result += f*integral(J, arrang, lambda_max, lambda_step,
		                         t->basis_a, t->basis_b, m->R);
	}

	const double end_time = wall_time();

	if (t->a == t->b)
	{
		result += t->basis_a->eigenval + centr_term(t->basis_a->l, mass, m->R);

		matrix_diag_set(c, t->a, result);
	}
	else
	{
		matrix_symm_set(c, t->a, t->b, result);
	}
}

/******************************************************************************

 Function driver(): performs the calculation of a single task using one thread
 from one process.

******************************************************************************/

void driver(const char arrang,
            const int grid_index,
            const int J,
            const int max_channel,
            const tasks list,
            const bool use_omp)
{
	multipole_set m;
	multipole_read(arrang, grid_index, &m);

	const int max_task = max_channel*(max_channel + 1)/2;

	matrix *c = matrix_alloc(max_channel, max_channel, true);

	#pragma omp parallel for default(none) shared(list, m, c) schedule(static) if(use_omp)
	for (int n = 0; n < max_task; ++n)
	{
		solve(arrang, J, &m, &list[n], c);
	}

	multipole_free(&m);
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
 *	Arrangement (a = 1, b = 2, c = 3) and atomic masses:
 */

	const char arrang = 96 + (int) file_keyword(stdin, "arrang", 1.0, 3.0, 1.0);

	pes_init_mass(stdin, 'a');
	pes_init_mass(stdin, 'b');
	pes_init_mass(stdin, 'c');

/*
 *	Total angular momentum, J:
 */

	const int J = (int) file_keyword(stdin, "J", 0.0, INF, 0.0);

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
 *	OpenMP:
 */

	const bool use_omp = (bool) file_keyword(stdin, "use_omp", 0.0, 1.0, 1.0);

/*
 *	Read all n basis functions from the disk and sort them in an ordered list of
 *	n*(n + 1)/2 tasks, i.e. the upper triangular part of the symmetric coupling
 *	matrix. Thus, a better workload of tasks per thread is achieved if use_omp =
 *	true.
 */

	const int max_channel = basis_count(arrang, J);

	basis *b_list = allocate(max_channel, sizeof(basis), true);

	for (int n = 0; n < max_channel; ++n)
	{
		basis_read(arrang, n, J, &b_list[n]);
	}

	const int max_task = max_channel*(max_channel + 1)/2;

	tasks *list = allocate(max_task, sizeof(tasks), true);

	int counter = 0;
	for (int n = 0; n < max_channel; ++n)
	{
		for (int m = n; m < max_channel; ++m)
		{
			list[counter].a = n;
			list[counter].basis_a = &b_list[n];

			list[counter].b = m;
			list[counter].basis_b = &b_list[m];

			++counter;
		}
	}

	ASSERT(counter == max_task);

/*
 *	Whereas each OpenMP thread handle the integration of each matrix element, MPI
 *	processes are used to handle each scattering grid point from R_min to R_max. 
 */

	mpi_set_tasks(scatt_grid_size);

/*
 *	Resolve all tasks:
 */

	if (mpi_rank() == 0)
	{
		printf("# MPI CPUs = %d, OpenMP threads = %d, num. of tasks = %d\n", mpi_comm_size(), max_threads(), max_task);
		printf("# CPU    thread    lambda    r (a.u.)    R (a.u.)    multipole (a.u.)    wall time (s)\n");
		printf("# ------------------------------------------------------------------------------------\n");
	}

	for (int n = mpi_first_task(); n <= mpi_last_task(); ++n)
	{
		extra_step:
		driver(arrang, n, J, max_task, list);

		if (n == mpi_last_task() && mpi_extra_task() > 0)
		{
			n = mpi_extra_task();
			goto extra_step;
		}
	}

	for (int n = 0; n < max_channel; ++n)
	{
		if (b_list[n].eigenvec != NULL) free(b_list[n].eigenvec);
	}

	free(b_list);
	free(list);

	mpi_end();
	return EXIT_SUCCESS;
}
