/******************************************************************************

 About
 -----

 This module is a collection of routines to handle the many kinds of masses
 needed in atomic and molecular physics.

******************************************************************************/

#include "file.h"
#include "nist.h"
#include "mass.h"

static double mass_a, mass_b, mass_c,
              mass_ab, mass_bc, mass_ac, mass_abc, mass_bca, mass_cab;

/******************************************************************************

 Function read_atomic_mass(): an auxiliary routine that scan a given input file
 searching for the first occurrence of a keyword of format "[key] = [value]".
 Where, [value] is expected as an atomic symbol, e.g. "1H", "2H", "4He", ...,
 "12C", etc. The actual atomic mass is returned if the keyword is found.

 NOTE: Lines starting by '#' are ignored.

******************************************************************************/

static double read_atomic_mass(FILE *input, const char key[])
{
	ASSERT(input != NULL)

	char *line = file_find_string(input, key);
	char *token = strtok(line, "=");

	while (token != NULL)
	{
		if (strstr(line, key) != NULL)
		{
			token = strtok(NULL, "=");
			token = trim(token);

			const isotope s = nist_isotope(token);

			/* NOTE: masses are returned in atomic units. */
			return nist_atomic_mass(s)*1822.8873843;
		}

		token = strtok(NULL, "=");
	}

	PRINT_ERROR("no entry '%s' found\n", key)
	exit(EXIT_FAILURE);
}

/******************************************************************************

 Function mass_init(): reads from a given input the atomic symbols for atoms A,
 B and C using the following format:

 mass_a = [symbol]
 mass_b = [symbol]
 mass_c = [symbol]

 Where, [symbol] is expected as "1H", "2H", "4He", ..., "12C", etc.

 The routine will store the actual masses for a later use by mass().

 NOTE: Lines starting by '#' are ignored.

******************************************************************************/

void mass_init(FILE *input)
{
	ASSERT(input != NULL)

	mass_a = read_atomic_mass(input, "mass_a");
	mass_b = read_atomic_mass(input, "mass_b");
	mass_c = read_atomic_mass(input, "mass_c");

	mass_ab = mass_a*mass_b/(mass_a + mass_b);
	mass_bc = mass_b*mass_c/(mass_b + mass_c);
	mass_ac = mass_a*mass_c/(mass_a + mass_c);

	mass_abc = mass_a*(mass_b + mass_c)/(mass_a + mass_b + mass_c);
	mass_bca = mass_b*(mass_c + mass_a)/(mass_a + mass_b + mass_c);
	mass_cab = mass_c*(mass_a + mass_b)/(mass_a + mass_b + mass_c);
}

/******************************************************************************

 Function mass(): return the many kind of masses and reduced masses for a tri
 atomic system, previously initialized by mass_init().

 NOTE: see the enum mass_case definition for details.

******************************************************************************/

double mass(const mass_case m)
{
	switch (m)
	{
		case atom_a:   return mass_a;
		case atom_b:   return mass_b;
		case atom_c:   return mass_c;
		case pair_ab:  return mass_ab;
		case pair_bc:  return mass_bc;
		case pair_ac:  return mass_ac;
		case arrang_a: return mass_abc;
		case arrang_b: return mass_bca;
		case arrang_c: return mass_cab;
		case total:    return (mass_a + mass_b + mass_c);

		default:
			PRINT_ERROR("invalid choice of enum mass_case, m = %d\n", m)
			exit(EXIT_FAILURE);
	}
}
