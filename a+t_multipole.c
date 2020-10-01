#include "modules/pes.h"
#include "modules/file.h"
#include "modules/math.h"
#include "modules/mpi_lib.h"
#include "modules/globals.h"

#define FORMAT "  %3d       %3d       %3d       %3d       %3d    %06f    %06f    %06f       % -8e         %f\n"

struct task
{
	double r1, r2, r3;
	int eta, m_eta, lambda;
};

/******************************************************************************

 Function legendre_integrand(): is an auxiliary routine which returns the atom-
 diatom Legendre multipole integrand for the inner integral in Eq. (22) of Ref.
 [1] at a given (r, R) Jacobi coordinate.

******************************************************************************/

double legendre_integrand(const double theta, void *params)
{
	struct task *p = (struct task *) params;

	const double v = pes_harmonics_multipole(p->eta,
	                                         p->m_eta,
	                                         p->r1, p->r2, p->r3, theta*180.0/M_PI);

	return v*math_legendre_poly(p->lambda, cos(theta))*sin(theta);
}

/******************************************************************************

 Function driver(): performs the calculation of a single task using one thread
 from one CPU.

******************************************************************************/

void driver(const int eta,
            const int m_eta, const double r3, pes_multipole *m, FILE *output)
{
	ASSERT(m != NULL)
	ASSERT(output != NULL)

	struct task job =
	{
		.r2 = m->R,
		.r3 = r3,
		.eta = eta,
		.m_eta = m_eta
	};

	for (int lambda = m->lambda_min; lambda <= m->lambda_max; ++lambda)
	{
		for (int n = 0; n < m->grid_size; ++n)
		{
			job.lambda = lambda;
			job.r1 = m->r_min + as_double(n)*m->r_step;

			const double start_time = wall_time();

			m->value[lambda][n] = math_qags(0.0, M_PI, &job, legendre_integrand);

			m->value[lambda][n] *= as_double(2*lambda + 1)/2.0;

			const double end_time = wall_time();

			#pragma omp critical
			{
				printf(FORMAT, mpi_rank(), thread_id(), lambda, eta, m_eta,
				       job.r1, job.r2, job.r3, m->value[lambda][n], end_time - start_time);
			}

			pes_multipole_write(m, output);
		}
	}
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
 *	Vibrational grid for BC (r1):
 */

	const int bc_grid_size = (int) file_keyword(stdin, "bc_grid_size", 1.0, INF, 500.0);

	const double bc_min = file_keyword(stdin, "bc_min", 0.0, INF, 0.5);

	const double bc_max = file_keyword(stdin, "bc_max", bc_min, INF, bc_min + 30.0);

	const double bc_step = (bc_max - bc_min)/as_double(bc_grid_size);

/*
 *	Vibrational grid for BC-D (r2):
 */

	const int bcd_grid_size = (int) file_keyword(stdin, "bcd_grid_size", 1.0, INF, 500.0);

	const double bcd_min = file_keyword(stdin, "bcd_min", 0.0, INF, 0.5);

	const double bcd_max = file_keyword(stdin, "bcd_max", bcd_min, INF, bcd_min + 30.0);

	const double bcd_step = (bcd_max - bcd_min)/as_double(bcd_grid_size);

/*
 *	Scattering grid for A + BC-D (R or r3):
 */

	const int scatt_grid_size = (int) file_keyword(stdin, "scatt_grid_size", 1.0, INF, 500.0);

	const double R_min = file_keyword(stdin, "R_min", 0.0, INF, 0.5);

	const double R_max = file_keyword(stdin, "R_max", R_min, INF, R_min + 30.0);

	const double R_step = (R_max - R_min)/as_double(scatt_grid_size);

/*
 *	Multipoles:
 */

	const int lambda_min = (int) file_keyword(stdin, "lambda_min", 0.0, INF, 0.0);

	const int lambda_max = (int) file_keyword(stdin, "lambda_max", 0.0, INF, 20.0);

	const int lambda_step = (int) file_keyword(stdin, "lambda_step", 1.0, INF, 2.0);

	ASSERT(lambda_max >= lambda_min)

	const int eta_min = (int) file_keyword(stdin, "eta_min", 0.0, INF, 0.0);

	const int eta_max = (int) file_keyword(stdin, "eta_max", 0.0, INF, 20.0);

	const int eta_step = (int) file_keyword(stdin, "eta_step", 1.0, INF, 2.0);

	ASSERT(eta_max >= eta_min)

/*
 *	OpenMP:
 */

	const bool use_omp = (bool) file_keyword(stdin, "use_omp", 0.0, 1.0, 1.0);

/*
 *	MPI:
 */

	mpi_set_tasks(scatt_grid_size);

/*
 *	Resolve all tasks:
 */

	if (mpi_rank() == 0)
	{
		printf("# MPI CPUs = %d, OpenMP threads = %d\n", mpi_comm_size(), max_threads());
		printf("# CPU    thread    lambda    eta    m_eta    r1 (a.u.)    r2 (a.u.)      R (a.u.)    multipole (a.u.)    wall time (s)\n");
		printf("# --------------------------------------------------------------------------------------------------------------------\n");
	}

	pes_multipole job;

	job.value = NULL;

	job.r_min = bc_min;
	job.r_max = bc_max;
	job.r_step = bc_step;

	job.lambda_min = lambda_min;
	job.lambda_max = lambda_max;
	job.lambda_step = lambda_step;

	job.grid_size = bc_grid_size;

	/* FIXME: when using Intel icc, n must be declared with an OMP private clause. */
	int n = 0;

	#pragma omp parallel for default(none) private(job, n) schedule(static) if(use_omp)
	for (n = mpi_first_task(); n <= mpi_last_task(); ++n)
	{
		extra_step:
		pes_multipole_init(&job);

		const double r3 = R_min + as_double(n)*R_step;

		FILE *output = pes_multipole_file('a', n, "wb", false);

		file_write(&r3, sizeof(double), 1, output);

		file_write(&eta_min, sizeof(int), 1, output);

		file_write(&eta_max, sizeof(int), 1, output);

		file_write(&eta_step, sizeof(int), 1, output);

		file_write(&bcd_min, sizeof(double), 1, output);

		file_write(&bcd_max, sizeof(double), 1, output);

		file_write(&bcd_step, sizeof(double), 1, output);

		file_write(&bcd_grid_size, sizeof(int), 1, output);

		for (int eta = eta_min; eta <= eta_max; eta += eta_step)
		{
			for (int m_eta = 0; m_eta <= eta; ++m_eta)
			{
				for (int m = 0; m < bcd_grid_size; ++m)
				{
					job.R = bcd_min + as_double(m)*bcd_step;

					driver(eta, m_eta, r3, &job, output);
				}
			}
		}

		file_close(&output);

		pes_multipole_free(&job);

		if (n == mpi_last_task() && mpi_extra_task() > 0)
		{
			n = mpi_extra_task();
			goto extra_step;
		}
	}

	mpi_end();
	return EXIT_SUCCESS;
}
