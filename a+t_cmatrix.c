#include "modules/pes.h"
#include "modules/file.h"
#include "modules/math.h"
#include "modules/matrix.h"
#include "modules/mpi_lib.h"
#include "modules/globals.h"

#include "utils.h"

#define FORMAT "  %3d       %3d       %5d       %5d       %06f    %f\n"

struct tasks
{
	int a, b;
	basis *bra, *ket;
};

typedef struct tasks tasks;

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

double percival_seaton(const int j1,
                       const int j2,
                       const int l1,
                       const int l2,
                       const int lambda, const int J)
{
	double result = pow(-1.0, j1 + j2 - J);

	result *= math_wigner_3j(l1, l2, lambda, 0, 0, 0);
	result *= math_wigner_3j(j1, j2, lambda, 0, 0, 0);
	result *= math_wigner_6j(j1, j2, lambda, l2, l1, J);
	result *= sqrt(as_double(2*j1 + 1)*as_double(2*j2 + 1)*as_double(2*l1 + 1)*as_double(2*l2 + 1));

	return result;
}

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
	int n_max = 0;

	if ((grid_size - 1)%3 != 0 && (grid_size - 2)%3 == 0) n_max = grid_size - 1;
	if ((grid_size - 1)%3 != 0 && (grid_size - 2)%3 != 0) n_max = grid_size - 2;
	if ((grid_size - 1)%3 == 0) n_max = grid_size - 3;

	ASSERT(n_max > 0)

	double sum = potential[0]*wavef_a[0]*wavef_b[0]
	           + potential[n_max]*wavef_a[n_max]*wavef_b[n_max];

	for (int n = 1; n < n_max; n += 3)
	{
		sum += 3.0*potential[n]*wavef_a[n]*wavef_b[n];

		sum += 3.0*potential[n + 1]*wavef_a[n + 1]*wavef_b[n + 1];

		sum += 2.0*potential[n + 2]*wavef_a[n + 2]*wavef_b[n + 2];
	}

	return grid_step*3.0*sum/8.0;
}

/******************************************************************************

 Function bc_braket(): performs the integral between two states <a| and |b>,
 subjected to a potential expanded in a set of multipolar expansion
 coefficients.

******************************************************************************/

double bc_braket(const int J,
                 const multipole_set *m, const basis *a, const basis *b)
{
	ASSERT(a->r_step == b->r_step)
	ASSERT(b->r_step == m->r_step)

	ASSERT(a->grid_size == b->grid_size)
	ASSERT(b->grid_size == m->grid_size)

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

 Function bc_cmatrix(): returns the coupling matrix element <a|V|b> of the BC
 subsystem subjected to a potential expanded in a set of multipolar expansion
 coefficients computed at R = r_abcd.

******************************************************************************/

double bc_cmatrix(const int j, const int a, const int b, const multipole_set *m)
{
	ASSERT(m != NULL)

	basis bra;
	bra.eigenvec = NULL;
	basis_read('a', a, j, &bra, false);

	basis ket;
	ket.eigenvec = NULL;
	basis_read('a', b, j, &ket, false);

	const double result = bc_braket(j, m, &bra, &ket);

	free(bra.eigenvec);
	free(ket.eigenvec);

	if (a == b)
		return result + centr_term(ket.l, pes_mass_abc('a'), m->R);
	else
		return result;
}

/******************************************************************************

 Function bc_braket(): performs the integral between two states <a| and |b>,
 subjected to a potential expanded in a set of multipolar expansion
 coefficients.

******************************************************************************/

double bcd_braket(const int J,
                  const multipole_set *m, const basis *a, const basis *b)
{
	ASSERT(a->r_step == b->r_step)
	ASSERT(b->r_step == m->r_step)

	ASSERT(a->grid_size == b->grid_size)
	ASSERT(b->grid_size == m->grid_size)

	double result = 0.0;
	for (int lambda = 0; lambda <= m->lambda_max; ++lambda)
	{
		for (int eta = 0; eta <= eta_max; ++eta)
		{
			for (int m_eta = -eta; m_eta <= eta; ++eta)
			{
				if (m->set[lambda].value == NULL) continue;

				const double f = percival_seaton(a->j, b->j, a->l, b->l, lambda, J);

				if (f == 0.0) continue;

				const double v = simpson(m->grid_size, m->r_step,
				                         m->set[lambda].value, a->eigenvec, b->eigenvec);
				result += v*f;
			}
		}
	}

	return result;
}

/******************************************************************************

 Function driver(): compute one task corresponding to one element of the coupling

******************************************************************************/

void driver(const int J, const double mass,
            const multipole_set *m_list, const tasks *t, matrix *c)
{
	const double start_time = wall_time();

	double result = bcd_braket(J, m_list, t->bra, t->ket);

	const double end_time = wall_time();

	if (t->a == t->b)
	{
		const double mass = pes_mass_abcd();

		result += t->basis_a->eigenval + centr_term(t->basis_a->l, mass, m_list->R);
		matrix_set_diag(c, t->a, result);
	}
	else
	{
		matrix_set_symm(c, t->a, t->b, result);
	}

	#pragma omp critical
	printf(FORMAT, mpi_rank(), thread_id(), t->a, t->b, m_list->R, end_time - start_time);
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
	pes_init_mass(stdin, 'd');

/*
 *	Total angular momentum, J:
 */

	const int J = (int) file_keyword(stdin, "J", 0.0, INF, 0.0);

/*
 *	Scattering grid:
 */

	const int grid_size = multipole_count(arrang);

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
			list[counter].bra = &b_list[n];

			list[counter].b = m;
			list[counter].ket = &b_list[m];

			++counter;
		}
	}

	ASSERT(counter == max_task);

/*
 *	Whereas each OpenMP thread handle the integration of each matrix element, MPI
 *	processes are used to handle each scattering grid point from R_min to R_max.
 */

	mpi_set_tasks(grid_size);

/*
 *	Resolve all tasks:
 */

	if (mpi_rank() == 0)
	{
		printf("# MPI CPUs = %d", mpi_comm_size());
		printf(", OMP threads = %d", max_threads());
		printf(", channels = %d", max_channel);
		printf(", grid points = %d\n", grid_size);

		printf("# CPU    thread       Ch. a       Ch. b       R (a.u.)    wall time (s)\n");
		printf("# ---------------------------------------------------------------------\n");
	}

	multipole_set m_list;

	for (int n = mpi_first_task(); n <= mpi_last_task(); ++n)
	{
		extra_step:
		multipole_read(arrang, n, &m_list);

		matrix *c = matrix_alloc(max_channel, max_channel, true);

		#pragma omp parallel for default(none) shared(list, m_list, c) schedule(static) if(use_omp)
		for (int n = 0; n < max_task; ++n)
		{
			driver(J, &m_list, &list[n], c);
		}

		coupling_write(arrang, n, J, false, c);

		matrix_free(c);
		multipole_free(&m_list);

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
