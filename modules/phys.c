/******************************************************************************

 About
 -----

 This module is a collection of special math functions for problems of amo and
 quantum physics.


 References
 ----------

 [1] R. B. Bernstein et al. Proc. R. Soc. Lond. A 1963 274 , 427-442
     doi: https://doi.org/10.1098/rspa.1963.0142

 [2] A. E. DePristo et al. J. Phys. B: Atom. Molec. Phys., Vol. 9, No. 3 (1976)
     doi:

 [3] D. D. Lopez-Duran et al. Com. Phys. Comm. 179 (2008) 821-838
     doi: https://doi.org/10.1016/j.cpc.2008.07.017

******************************************************************************/

#include "phys.h"

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_coupling.h>
#include <gsl/gsl_sf_legendre.h>

/******************************************************************************

 Function wigner_3j():

******************************************************************************/

double phys_wigner_3j(const int a, const int b, const int c,
                      const int d, const int e, const int f)
{
	return gsl_sf_coupling_3j(2*a, 2*b, 2*c, 2*d, 2*e, 2*f);
}

/******************************************************************************

 Function wigner_6j():

******************************************************************************/

double phys_wigner_6j(const int a, const int b, const int c,
                      const int d, const int e, const int f)
{
	return gsl_sf_coupling_6j(2*a, 2*b, 2*c, 2*d, 2*e, 2*f);
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

double phys_percival_seaton(const int spin_mult, const int j1, const int j2,
                            const int l1, const int l2, const int lambda, const int J)
{
	double result = 0.0;

	switch (spin_mult)
	{
		case 1: /* Eq. (3) of Ref. [2] */
			result  = pow(-1.0, j1 + j2 - J);
			result *= phys_wigner_3j(l1, l2, lambda, 0, 0, 0);
			result *= phys_wigner_3j(j1, j2, lambda, 0, 0, 0);
			result *= phys_wigner_6j(j1, j2, lambda, l2, l1, J);
			result *= sqrt(as_double(2*j1 + 1)*as_double(2*j2 + 1)*as_double(2*l1 + 1)*as_double(2*l2 + 1));
		break;

		case 2: /* Eq. (C.21) of Ref. [3] */
		break;

		case 3: /* Eq. (C.22) of Ref. [3] */
		break;

		default:
			PRINT_ERROR("invalid spin multiplicity = %d\n", spin_mult)
			exit(EXIT_FAILURE);
	}

	return result;
}
