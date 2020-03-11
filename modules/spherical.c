#include "phys.h"
#include "spherical.h"

#include <gsl/gsl_sf_legendre.h>

/******************************************************************************

 Function spherical_harmonics(): returns the spherical harmonics for angular
 momentum l, with projection m, at (theta, phi).

******************************************************************************/

double spherical_harmonics(const int l, const int m, const spherical *r)
{
	ASSERT(l >= abs(m))

	const double x
		= cos(r->theta*M_PI/180.0);

	const double m_phase
		= (m > 0? pow(-1.0, m) : 1.0);

	const double phi_wavef
		= exp(as_double(m)*r->phi*M_PI/180.0)/sqrt(2.0*M_PI);

	/* NOTE: see equation 1.43 (pag. 8) of Angular Momentum by Richard N. Zare. */
	return m_phase*gsl_sf_legendre_sphPlm(l, abs(m), x)*phi_wavef;
}

/******************************************************************************

 Function spherical_harmonics_coupl(): returns a coupled spherical harmonics
 built upon the combination of n uncoupled ones for a total angular momentum
 l[0] + l[1] + ... + l[n], with projection m[0] + m[1] + ... + m[n], at r[0],
 r[1], ..., r[n].

******************************************************************************/

double spherical_harmonics_coupl(const int n
                                 const int l[],
                                 const int m[],
                                 const spherical r[])
{
	ASSERT(n > 0)
	ASSERT(n < 3)
	ASSERT(l != NULL)
	ASSERT(m != NULL)
	ASSERT(r != NULL)

	double result = 0.0;

	switch (n)
	{
		case 1:
		{
			result = spherical_harmonics(l[0], m[0], &r[0]);
		}
		break;

		case 2:
		{
			const int L = l[0] + l[1];

			for (int a = -l[0]; a <= l[0]; ++a)
			{
				for (int b = -l[1]; b <= l[1]; ++b)
				{
					const int M = a + b;

					const double c
						= phys_clebsch_gordan(l[0], l[1], L, a, b, M);

					if (c == 0.0) continue;

					const double Ya
						= spherical_harmonics(l[0], a, &r[0]);

					const double Yb
						= spherical_harmonics(l[1], b, &r[1]);

					result += c*Ya*Yb;
				}
			}
		}
		break;
	}

	return result;
}
