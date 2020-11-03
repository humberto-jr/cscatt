#include "modules/nist.h"
#include "modules/globals.h"

int main(int argc, char *argv[])
{
	if (argc < 2)
	{
		PRINT_ERROR("%d arguments given. Usage: ./mass.out [atom a] [atom b] [...]\n", argc - 1)
		return EXIT_FAILURE;
	}

	for (int n = 1; n < argc; ++n)
	{
		const enum isotope a = nist_isotope(argv[n]);

		/* NOTE: atomic mass unit. */
		const double amu = nist_atomic_mass(a);

		/* NOTE: atomic unit mass. */
		const double aum = nist_atomic_mass(a)*1822.8873843;

		printf("%6s\t %f\t %f\n", nist_atomic_symbol(a), amu, aum);
	}

	return EXIT_SUCCESS;
}
