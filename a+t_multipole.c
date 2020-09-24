#include "modules/pes.h"
#include "modules/file.h"
#include "modules/math.h"
#include "modules/mpi_lib.h"
#include "modules/globals.h"

#define FORMAT "  %3d       %3d       %3d       %3d       %3d    %06f    %06f    %06f       % -8e         %f\n"

struct tasks
{
	double r_bc, r_bcd, r_abcd, result;
	int lambda, eta, m_eta, grid_index;
	pes_multipole_set m;
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
				int n;
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

 Function legendre_integrand(): is an auxiliary routine which returns the atom-
 diatom Legendre multipole integrand for the inner integral in Eq. (22) of Ref.
 [1] at a given (r, R) Jacobi coordinate.

******************************************************************************/

double legendre_integrand(const double theta, void *params)
{
	struct tasks *p = (struct tasks *) params;

	const double v = pes_harmonics_multipole(p->eta,
	                                         p->m_eta,
	                                         p->r_bc,
	                                         p->r_bcd,
	                                         p->r_abcd,
	                                         theta*180.0/M_PI);

	return v*math_legendre_poly(p->lambda, cos(theta))*sin(theta);
}

int read_all()
{
	FILE *output = pes_multipole_file('?', n, "wb", false);

	const double r_abcd = abcd_rmin + as_double(n)*abcd_rstep;

	file_write(&r_abcd, sizeof(double), 1, output);

	file_write(&eta_max, sizeof(int), 1, output);

	for (int eta = eta_min; eta <= eta_max; eta += eta_step)
	{
		for (int m_eta = 0; m_eta <= eta; ++eta)
		{
			file_write(&eta, sizeof(int), 1, output);
			file_write(&m_eta, sizeof(int), 1, output);
			file_write(&bcd_grid_size, sizeof(int), 1, output);

			for (int p = 0; p < bcd_grid_size; ++p)
			{
				m.R = bcd_rmin + as_double(p)*bcd_rstep;

				pes_multipole_write(&m, output);
			}
		}
	}

	file_close(&output);
}

/******************************************************************************

 Function driver(): performs the calculation of a single task using one thread
 from one CPU.

******************************************************************************/

void driver(const int n)
{
	FILE *input = pes_multipole_file('?', n, "rb", false);

	double R;
	file_read(&R, sizeof(double), 1, input, 0);

	int eta_max;
	file_read(&eta_max, sizeof(int), 1, input, 0);

	int n_max;
	file_read(&n_max, sizeof(int), 1, input, 0);

	int grid_size;
	file_write(&grid_size, sizeof(int), 1, output);

	tasks *job = allocate(max_set, sizeof(tasks), true);

	for (int n = 0; n < n_max; ++n)
	{
		int eta;
		file_read(&eta, sizeof(int), 1, input, 0);

		int m_eta;
		file_read(&m_eta, sizeof(int), 1, input, 0);

		job[n].m = pes_multipole_read_all(grid_size, input);
	}

	file_close(&input);

	FILE *output = pes_multipole_file('?', n, "wb", false);

	file_write(&R, sizeof(double), 1, output);

	file_write(&eta_max, sizeof(int), 1, output);

	file_write(&n_max, sizeof(int), 1, output);

	file_write(&grid_size, sizeof(int), 1, output);

	for (int eta = eta_min; eta <= eta_max; eta += eta_step)
	{
		for (int m_eta = 0; m_eta <= eta; ++eta)
		{
			file_write(&eta, sizeof(int), 1, output);
			file_write(&m_eta, sizeof(int), 1, output);

			for (int p = 0; p < bcd_grid_size; ++p)
			{
				m.R = bcd_rmin + as_double(p)*bcd_rstep;

				pes_multipole_write(&m, output);
			}

			++counter;
		}
	}


	file_close(&input);
}

/*
void driver(struct tasks *job)
		
{
	const double start_time = wall_time();

	job->result = math_qags(0.0, M_PI, &job, legendre_integrand);
	job->result = as_double(2*job->lambda + 1)*job->result/2.0;

	const double end_time = wall_time();

	#pragma omp critical
	{
		printf(FORMAT, mpi_rank(), thread_id(), job->lambda, job->eta, job->m_eta,
		       job->r_bc, job->r_bcd, job->r_abcd, job->result, end_time - start_time);
	}
}
*/

/******************************************************************************

 Function sort_results(): CPU 0 writes the respective multipoles as function of
 r per lambda using binary format and filename labeled by the grid index, n (of
 R) and arrangement.

 TODO: working here, needs to call pes_multipole_write().

******************************************************************************/

void sort_results(const int max_task,
                  const int eta_max,
                  const int bcd_grid_size,
                  const double bcd_rmin,
                  const double bcd_rmax,
                  const double bcd_rstep,
                  const struct tasks list[],
                  const pes_multipole_set *m)
{
	if (mpi_rank() != 0) return;

	FILE *output = NULL;
	int last_index = -1, last_lambda = -1, last_eta = -1, last_m = -1;

	for (int n = 0; n < max_task; ++n)
	{
		if (n > 0 && list[n].grid_index != last_index)
		{
			file_close(&output);
			output = multipole_file(list[n].grid_index, "wb");

			file_write(&list[n].r_abcd, sizeof(double), 1, output);

			file_write(&bc_rmin, sizeof(double), 1, output);
			file_write(&bc_rmax, sizeof(double), 1, output);
			file_write(&bc_rstep, sizeof(double), 1, output);

			file_write(&bcd_rmin, sizeof(double), 1, output);
			file_write(&bcd_rmax, sizeof(double), 1, output);
			file_write(&bcd_rstep, sizeof(double), 1, output);

			file_write(&eta_max, sizeof(int), 1, output);
			file_write(&lambda_max, sizeof(int), 1, output);
			file_write(&bc_grid_size, sizeof(int), 1, output);
			file_write(&bcd_grid_size, sizeof(int), 1, output);

			last_lambda = -1;
			last_eta = -1;
		}

		if (list[n].lambda != last_lambda || list[n].eta != last_eta)
		{
			last_m = -1;
		}

		if (list[n].m_eta != last_m)
		{
			file_write(&list[n].lambda, sizeof(int), 1, output);
			file_write(&list[n].eta, sizeof(int), 1, output);
			file_write(&list[n].m_eta, sizeof(int), 1, output);
		}

		m->set[list[n].lambda].value = list[n].result;

		last_index = list[n].grid_index;
		last_lambda = list[n].lambda;
		last_eta = list[n].eta;
		last_m = list[n].m_eta;
	}

	file_close(&output);
}

/******************************************************************************
******************************************************************************/

int main(int argc, char *argv[])
{
	mpi_init(argc, argv);

	file_init_stdin(argv[1]);

/*
 *	Atomic masses:
 */

	pes_init_mass(stdin, 'a');
	pes_init_mass(stdin, 'b');
	pes_init_mass(stdin, 'c');
	pes_init_mass(stdin, 'd');

/*
 *	Vibrational grid for BC:
 */

	const int bc_grid_size = (int) file_keyword(stdin, "bc_grid_size", 1.0, INF, 500.0);

	const double bc_rmin = file_keyword(stdin, "bc_rmin", 0.0, INF, 0.5);

	const double bc_rmax = file_keyword(stdin, "bc_rmax", bc_rmin, INF, bc_rmin + 30.0);

	const double bc_rstep = (bc_rmax - bc_rmin)/as_double(bc_grid_size);

/*
 *	Vibrational grid for BC-D:
 */

	const int bcd_grid_size = (int) file_keyword(stdin, "bcd_grid_size", 1.0, INF, 500.0);

	const double bcd_rmin = file_keyword(stdin, "bcd_rmin", 0.0, INF, 0.5);

	const double bcd_rmax = file_keyword(stdin, "bcd_rmax", bcd_rmin, INF, bcd_rmin + 30.0);

	const double bcd_rstep = (bc_rmax - bc_rmin)/as_double(bc_grid_size);

/*
 *	Scattering grid:
 */

	const int abcd_grid_size = (int) file_keyword(stdin, "abcd_grid_size", 1.0, INF, 500.0);

	const double abcd_rmin = file_keyword(stdin, "abcd_rmin", 0.0, INF, 0.5);

	const double abcd_rmax = file_keyword(stdin, "abcd_rmax", abcd_rmin, INF, abcd_rmin + 30.0);

	const double abcd_rstep = (abcd_rmax - abcd_rmin)/as_double(abcd_grid_size);

/*
 *	Multipoles:
 */

	const int lambda_min = (int) file_keyword(stdin, "lambda_min", 0.0, INF, 0.0);

	const int lambda_max = (int) file_keyword(stdin, "lambda_max", 0.0, INF, 20.0);

	const int lambda_step = (int) file_keyword(stdin, "lambda_step", 1.0, INF, 2.0);

	const int lambda_grid_size = (lambda_max - lambda_min)/lambda_step + 1;

	ASSERT(lambda_max >= lambda_min)

	const int eta_min = (int) file_keyword(stdin, "eta_min", 0.0, INF, 0.0);

	const int eta_max = (int) file_keyword(stdin, "eta_max", 0.0, INF, 20.0);

	const int eta_step = (int) file_keyword(stdin, "eta_step", 1.0, INF, 2.0);

	const int eta_grid_size = (eta_max - eta_min)/eta_step + 1;

	ASSERT(eta_max >= eta_min)

	int counter = 0;
	for (int eta = eta_min; eta <= eta_max; eta += eta_step)
	{
		/* NOTE: only the positive m-projections are considered. */
		counter += eta + 1;
	}

	const int m_grid_size = counter;

/*
 *	OpenMP:
 */

	const bool use_omp = (bool) file_keyword(stdin, "use_omp", 0.0, 1.0, 1.0);

/*
 *	Sort each task (r_bc, r_bcd, r_abcd, lambda, eta, m_eta) in an ordered list:
 */

	pes_multipole_set m =
	{
		.r_min = bc_rmin,
		.r_max = bc_rmax,
		.r_step = bc_rstep,
		.lambda_max = lambda_max, 
		.grid_size = bc_grid_size
	};

	pes_multipole_init_set(&m, bc_grid_size, lambda_max, lambda_step);

	for (int n = 0; n < abcd_grid_size; ++n)
	{
		FILE *output = pes_multipole_file('?', n, "wb", false);

		const double r_abcd = abcd_rmin + as_double(n)*abcd_rstep;

		file_write(&r_abcd, sizeof(double), 1, output);

		file_write(&eta_max, sizeof(int), 1, output);

		counter = eta_grid_size*m_grid_size;
		file_write(&counter, sizeof(int), 1, output);

		file_write(&bcd_grid_size, sizeof(int), 1, output);

		counter = 0;
		for (int eta = eta_min; eta <= eta_max; eta += eta_step)
		{
			for (int m_eta = 0; m_eta <= eta; ++eta)
			{
				file_write(&eta, sizeof(int), 1, output);
				file_write(&m_eta, sizeof(int), 1, output);

				for (int p = 0; p < bcd_grid_size; ++p)
				{
					m.R = bcd_rmin + as_double(p)*bcd_rstep;

					pes_multipole_write(&m, output);
				}

				++counter;
			}
		}

		file_close(&output);

		ASSERT(counter == eta_grid_size*m_grid_size)
	}

/*
 *	MPI:
 */

	mpi_set_tasks(abcd_grid_size);

/*
 *	Resolve all tasks:
 */

	if (mpi_rank() == 0)
	{
		printf("# MPI CPUs = %d, OpenMP threads = %d, num. of tasks = %d\n", mpi_comm_size(), max_threads(), max_task);
		printf("# CPU    thread    lambda    eta    m_eta    r_bc (a.u.)    r_bcd (a.u.)    r_abcd (a.u.)    multipole (a.u.)    wall time (s)\n");
		printf("# ----------------------------------------------------------------------------------------------------------------------------\n");
	}

	/* FIXME: when using Intel icc, n must be declared with an OMP private clause. */
	int n = 0;

	#pragma omp parallel for default(none) private(n) shared(list) schedule(static) if(use_omp)
	for (n = mpi_first_task(); n <= mpi_last_task(); ++n)
	{
		extra_step:
		driver(n, eta_max, lambda_step, use_omp);

		if (n == mpi_last_task() && mpi_extra_task() > 0)
		{
			n = mpi_extra_task();
			goto extra_step;
		}
	}

/*	send_results(max_task, list);

	pes_multipole_set m;

	if (mpi_rank() == 0)
	{
		m.r_min = bc_rmin;
		m.r_max = bc_rmax;
		m.r_step = bc_rstep;
		m.lambda_max = lambda_max; 
		m.grid_size = bc_grid_size;
		m.set = allocate(lambda_max, sizeof(pes_multipole), true);

		for (int lambda = 0; lambda <= lambda_max; lambda += lambda_step)
		{
			m.set[lambda].value = allocate(bc_grid_size, sizeof(double), false);
		}
	}

	sort_results(max_task,
	             eta_max,
	             bcd_grid_size,
	             bcd_rmin,
	             bcd_rmax,
	             bcd_rstep,
	             list,
	             &m);


	free(list);
*/

	mpi_end();
	return EXIT_SUCCESS;
}
