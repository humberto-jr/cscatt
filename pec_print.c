#include "modules/pes.h"
#include "modules/file.h"
#include "modules/globals.h"

int main(int argc, char *argv[])
{
	ASSERT(argc > 1)
	file_init_stdin(argv[1]);

/*
 *	Rotational quantum numbers, j:
 */

	const int j_min = (int) file_keyword(stdin, "j_min", 0.0, INF, 0.0);

	const int j_max = (int) file_keyword(stdin, "j_max", 0.0, INF, 0.0);

	const int j_step = (int) file_keyword(stdin, "j_step", 1.0, INF, 1.0);

	ASSERT(j_max >= j_min)

/*
 *	Vibrational grid:
 */

	const int n_max = (int) file_keyword(stdin, "rovib_grid_size", 1.0, INF, 100.0);

	const double r_min = file_keyword(stdin, "r_min", 0.0, INF, 0.5);

	const double r_max = file_keyword(stdin, "r_max", r_min, INF, r_min + 30.0);

	const double r_step = (r_max - r_min)/as_double(n_max);

/*
 *	Energy scale and shifts:
 */

	const double shift = file_keyword(stdin, "energy_shift", -INF, INF, 0.0);

	const double scale = file_keyword(stdin, "energy_scale", -INF, INF, 1.0);

	const double inf = file_keyword(stdin, "inf", r_max, INF, 1000.0);

	pes_set_inf(inf);

/*
 *	Arrangement (a = 1, b = 2, c = 3), atomic masses and PES:
 */

	const char arrang = 96 + (int) file_keyword(stdin, "arrang", 1.0, 3.0, 1.0);

	pes_init_mass(stdin, 'a');
	pes_init_mass(stdin, 'b');
	pes_init_mass(stdin, 'c');

	pes_init();

	double (*pec)(const int, const double);

	switch (arrang)
	{
		case 'a': pec = pes_bc; break;
		case 'b': pec = pes_ac; break;
		case 'c': pec = pes_ab; break;

		default:
			PRINT_ERROR("invalid arrangement %c\n", arrang)
			exit(EXIT_FAILURE);
	}

	printf("# Arrangement = %c\n", arrang);
	printf("# j = (%d, %d, %d)\n", j_min, j_max, j_step);
	printf("# r = (%f, %f, %f)\n", r_min, r_max, r_step);

	printf("#       \t");

	for (int j = j_min; j <= j_max; j += j_step)
	{
		printf("  %12d", j);
	}

	printf("\n");

	for (int n = 0; n < n_max; ++n)
	{
		const double r = r_min + as_double(n)*r_step;

		printf("% 6f\t", r);

		for (int j = j_min; j <= j_max; j += j_step)
		{
			printf(" % -8e", (pec(j, r) + shift)*scale);
		}

		printf("\n");
	}

	return EXIT_SUCCESS;
}
