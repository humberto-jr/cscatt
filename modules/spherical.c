	#include "phys.h"
	#include "spherical.h"

	#include <gsl/gsl_sf_legendre.h>

	/******************************************************************************

	 Function spherical_harmonics(): returns the spherical harmonics for angular
	 momentum l, with projection m, at (theta, phi).

	******************************************************************************/

	double spherical_harmonics(const int l, const int m,
	                           const double theta, const double phi)
	{
		ASSERT(m >= 0)
		ASSERT(l >= m)

		const double x = cos(theta*M_PI/180.0);
		const double y = as_double(m)*(phi*M_PI/180.0);

		/* NOTE: see equation 1.43 (pag. 8) of Angular Momentum by Richard N. Zare. */
		return pow(-1.0, m)*gsl_sf_legendre_sphPlm(l, m, x)*pow(2.0*M_PI, -0.5)*exp(y);
	}

	/******************************************************************************

	 Function spherical_harmonics_ab(): returns the coupled spherical harmonics for
	 angular momentum l[0] + l[1], with projection m[0] + m[1], at (theta, phi);
	 built upon two spherical harmonics, |a> and |b>. Where, index 0 labels all
	 inputs for function a and index 1 those for function b.

	******************************************************************************/

	double spherical_harmonics_ab(const int l[], const int m[],
	                              const double theta[], const double phi[])
	{
		ASSERT(l != NULL)
		ASSERT(m != NULL)
		ASSERT(phi != NULL)
		ASSERT(theta != NULL)

		double result = 0.0;
		const int L = l[0] + l[1];

		for (int m_0 = 0; m_0 <= m[0]; ++m_0)
		{
			for (int m_1 = 0; m_1 <= m[1]; ++m_1)
			{
				const double y_0
					= spherical_harmonics(l[0], m_0, theta[0], phi[0]);

				const double y_1
					= spherical_harmonics(l[1], m_1, theta[1], phi[1]);

				const int M = m_0 + m_1;

				result += phys_clebsch_gordan(l[0], l[1], L, m_0, m_1, M)*y_0*y_1;
			}
		}

		return result;
	}
#endif
