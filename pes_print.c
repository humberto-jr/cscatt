#include "modules/pes.h"
#include "modules/file.h"
#include "modules/globals.h"

#define FORMAT "% 6f\t % 6f\t % 6f\t % -8e % -8e % -8e\n"

int main(int argc, char *argv[])
{
	ASSERT(argc > 1)
	file_init_stdin(argv[1]);

/*
 *	Vibrational grid, r:
 */

	const int rovib_grid_size = (int) file_keyword(stdin, "rovib_grid_size", 1.0, INF, 100.0);

	const double r_min = file_keyword(stdin, "r_min", 0.0, INF, 0.5);

	const double r_max = file_keyword(stdin, "r_max", r_min, INF, r_min + 30.0);

	const double r_step = (r_max - r_min)/as_double(rovib_grid_size);

/*
 *	Scattering grid, R:
 */

	const int scatt_grid_size = (int) file_keyword(stdin, "scatt_grid_size", 1.0, INF, 100.0);

	const double R_min = file_keyword(stdin, "R_min", 0.0, INF, 0.5);

	const double R_max = file_keyword(stdin, "R_max", R_min, INF, R_min + 50.0);

	const double R_step = (R_max - R_min)/as_double(scatt_grid_size);

/*
 *	Angular grid, theta:
 */

	const int theta_grid_size = (int) file_keyword(stdin, "theta_grid_size", 1.0, INF, 36.0);

	const double theta_min = file_keyword(stdin, "theta_min", 0.0, 180.0, 0.0);

	const double theta_max = file_keyword(stdin, "theta_max", theta_min, 180.0, 180.0);

	const double theta_step = (theta_max - theta_min)/as_double(theta_grid_size);

/*
 *	Energy scale and shifts:
 */

	const double shift = file_keyword(stdin, "energy_shift", -INFINITY, INFINITY, 0.0);

	const double scale = file_keyword(stdin, "energy_scale", -INFINITY, INFINITY, 1.0);

/*
 *	Atomic masses and PES:
 */

	pes_init_mass(stdin, 'a');
	pes_init_mass(stdin, 'b');
	pes_init_mass(stdin, 'c');

	pes_init();

/*
 *	Print output:
 */

	printf("# r     = (%f, %f, %f)\n", r_min, r_max, r_step);
	printf("# R     = (%f, %f, %f)\n", R_min, R_max, R_step);
	printf("# theta = (%f, %f, %f)\n", theta_min, theta_max, theta_step);
	printf("# Shift = %f\n", shift);
	printf("# Scale = %f\n", scale);

	printf("#\n");

	for (int n = 0; n < rovib_grid_size; ++n)
	{
		const double r = r_min + as_double(n)*r_step;

		for (int m = 0; m < scatt_grid_size; ++m)
		{
			const double R = R_min + as_double(m)*R_step;

			for (int p = 0; p <= theta_grid_size; ++p)
			{
				const double theta = theta_min + as_double(p)*theta_step;

				const double a = (pes_abc('a', r, R, theta) + shift)*scale;
				const double b = (pes_abc('b', r, R, theta) + shift)*scale;
				const double c = (pes_abc('c', r, R, theta) + shift)*scale;

				printf(FORMAT, r, R, theta, a, b, c);

				if (p == 0 && theta_step == 0.0) break;
			}

			if (rovib_grid_size > 1 && scatt_grid_size > 1 && theta_step != 0.0) printf("\n");
			if (m == 0 && R_step == 0.0) break;
		}

		if (scatt_grid_size > 1 && R_step != 0.0) printf("\n");
		if (n == 0 && r_step == 0.0) break;
	}

	return EXIT_SUCCESS;
}
