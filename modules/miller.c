/******************************************************************************

 About
 -----

 This module is a collection of routines applying the algorithms published by
 William H. Miller and co-workers on regarding atom-molcule collisions.


 References
 ----------

 [1] W. H. Miller. J. Chem. Phys. 1, 50 (1969)
     doi:

******************************************************************************/

#include "pes.h"
#include "dvr.h"
#include "mass.h"
#include "matrix.h"
#include "miller.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_integration.h>

#define GSL_MAX_WORKSPACE 5000

struct vib_params
{
	const int lambda;
	const char arrang;
	const double r_min, r_max, R;
	const matrix *wavef_a, *wavef_b;
};

struct rot_params
{
	const int lambda;
	jacobi_coor x;
};

/******************************************************************************

 Function rot_integrand(): is an auxiliary routine which return the integrand
 for the inner integral in Eq. (22) of Ref. [1] at a given (r, R) Jacobi
 coordinate.

******************************************************************************/

static inline double rot_integrand(const double theta, void *params)
{
	struct rot_params *p = (struct rot_params *) params;

	p->x.theta = theta*180.0/M_PI;
	const double pot_energy = pes(&p->x) - pec(p->x.arrang, p->x.r);

	/* Eq. (22) of Ref. [1], inner integrand */
	return pot_energy*gsl_sf_legendre_Pl(p->lambda, cos(theta))*sin(theta);
}

/******************************************************************************

 Function miller_jcp69_rot_integral(): return the inner integral (in theta) of
 the PES at a given (r, R) Jacobi coordinate, as shown in Eq. (22) of Ref. [1].

 NOTE: Integration over theta in [0, pi] performed by the QAG algorithm.

******************************************************************************/

double miller_jcp69_rot_integral(const char arrang,
                                 const int lambda, const double r, const double R)
{
	struct rot_params p =
	{
		.lambda = lambda,
		.x.arrang = arrang,
		.x.r = r,
		.x.R = R
	};

	gsl_function f =
	{
		.params = &p,
		.function = &rot_integrand
	};

	double error = 0.0, result = 0.0, factor = 1.0, theta_max = M_PI;

	if (arrang == 'a' && mass(atom_b) == mass(atom_c))
	{
		 factor = 2.0;
		 theta_max = M_PI/2.0;
	}

	if (arrang == 'b' && mass(atom_c) == mass(atom_a))
	{
		 factor = 2.0;
		 theta_max = M_PI/2.0;
	}

	if (arrang == 'c' && mass(atom_a) == mass(atom_b))
	{
		 factor = 2.0;
		 theta_max = M_PI/2.0;
	}

	gsl_integration_workspace *work
		= gsl_integration_workspace_alloc(GSL_MAX_WORKSPACE);

	ASSERT(work != NULL)

	const int info = gsl_integration_qags(&f, 0.0, theta_max, 1.0e-6, 0.0,
	                                      GSL_MAX_WORKSPACE, work, &result, &error);

	if (info != 0)
	{
		PRINT_ERROR("QAG %s for r = %f, R = %f; error = %f\n", gsl_strerror(info), r, R, error)
	}

	return as_double(2*lambda + 1)*factor*result/2.0;
}

/******************************************************************************

 Function vib_integrand(): is an auxiliary routine which return the integrand
 for the outer (vibrational) integral in Eq. (22) of Ref. [1] at a given R
 Jacobi coordinate.

 NOTE: any factor r^2 of Eq. (22) is not included as it is absorbed in the
 definition of vibrational wavefunctions.

******************************************************************************/

static inline double vib_integrand(const double r, void *params)
{
	const struct vib_params *p
		= (struct vib_params *) params;

	const double wavef_a
		= dvr_fgh_wavef(p->wavef_a, 0, p->r_min, p->r_max, r);

	const double wavef_b
		= dvr_fgh_wavef(p->wavef_b, 0, p->r_min, p->r_max, r);

	const double rot_integral
		= miller_jcp69_rot_integral(p->arrang, p->lambda, r, p->R);

	/* Eq. (22) of Ref. [1], outer integrand */
	return wavef_a*rot_integral*wavef_b;
}

/******************************************************************************

 Function miller_jcp69_vib_integral(): return the lambda-dependent outer
 (vibrational) integral (in r) at a given R Jacobi coordinate as shown in Eq.
 (22) of Ref. [1]

 NOTE: integration over r in [r_min, r_max] performed by the QAG algorithm.

 NOTE: r_min and r_max are expected the same for both vib. wavefunctions.

******************************************************************************/

double miller_jcp69_vib_integral(const int lambda,
                                 const char arrang,
                                 const matrix *wavef_a,
                                 const matrix *wavef_b,
                                 const double r_min,
                                 const double r_max,
                                 const double R)
{
	struct vib_params p =
	{
		.R = R,
		.r_min = r_min,
		.r_max = r_max,
		.lambda = lambda,
		.arrang = arrang,
		.wavef_a = wavef_a,
		.wavef_b = wavef_b
	};

	gsl_function f =
	{
		.params = &p,
		.function = &vib_integrand
	};

	double error = 0.0, result = 0.0;

	gsl_integration_workspace *work
		= gsl_integration_workspace_alloc(GSL_MAX_WORKSPACE);

	ASSERT(work != NULL)

	gsl_set_error_handler_off();

	const int info = gsl_integration_qags(&f, r_min, r_max, 1.0e-6, 0.0,
	                                      GSL_MAX_WORKSPACE, work, &result, &error);

	gsl_integration_workspace_free(work);

	if (info != 0)
	{
		PRINT_ERROR("QAG %s for R = %f; error = %f\n", gsl_strerror(info), R, error)
	}

	return result;
}
