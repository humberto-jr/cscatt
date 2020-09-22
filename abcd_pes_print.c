#include "modules/pes.h"
#include "modules/file.h"
#include "modules/globals.h"

#define FORMAT "% 6f\t % 6f\t % 6f\t % -8e\n"

int main(int argc, char *argv[])
{
	ASSERT(argc > 1)
	file_init_stdin(argv[1]);

/*
 *	BC-D geometry:
 */

	const double r1 = file_keyword(stdin, "r1", 0.0, INF, 2.0);

	const double r2 = file_keyword(stdin, "r2", 0.0, INF, 2.0);

	const double theta12 = file_keyword(stdin, "theta12", 0.0, INF, 2.0);

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
 *	Angular grid, theta:
 */

	const int phi_grid_size = (int) file_keyword(stdin, "phi_grid_size", 1.0, INF, 72.0);

	const double phi_min = file_keyword(stdin, "phi_min", 0.0, 360.0, 0.0);

	const double phi_max = file_keyword(stdin, "phi_max", phi_min, 360.0, 0.0);

	const double phi_step = (phi_max - phi_min)/as_double(phi_grid_size);

/*
 *	Energy scale and shifts:
 */

	const double shift = file_keyword(stdin, "energy_shift", -INF, INF, 0.0);

	const double scale = file_keyword(stdin, "energy_scale", -INF, INF, 1.0);

/*
 *	Atomic masses and PES:
 */

	pes_init_mass(stdin, 'a');
	pes_init_mass(stdin, 'b');
	pes_init_mass(stdin, 'c');
	pes_init_mass(stdin, 'd');

	pes_init();

/*
 *	Print output:
 */

	printf("# r1      = %f\n", r1);
	printf("# r2      = %f\n", r2);
	printf("# theta12 = %f\n", theta12);
	printf("# R       = (%f, %f, %f)\n", R_min, R_max, R_step);
	printf("# theta   = (%f, %f, %f)\n", theta_min, theta_max, theta_step);
	printf("# phi     = (%f, %f, %f)\n", phi_min, phi_max, phi_step);
	printf("# Shift   = %f\n", shift);
	printf("# Scale   = %f\n", scale);

	printf("#\n");

	for (int n = 0; n < scatt_grid_size; ++n)
	{
		const double R = R_min + as_double(n)*R_step;

		for (int m = 0; m < theta_grid_size; ++m)
		{
			const double theta = theta_min + as_double(m)*theta_step;

			for (int p = 0; p <= phi_grid_size; ++p)
			{
				const double phi = phi_min + as_double(p)*phi_step;

				const double a = (pes_abcd(r1, r2, R, theta12, theta, phi) + shift)*scale;

				printf(FORMAT, R, theta, phi, a);

				if (p == 0 && phi_step == 0.0) break;
			}

			if (scatt_grid_size > 1 && theta_grid_size > 1 && phi_step != 0.0) printf("\n");
			if (m == 0 && theta_step == 0.0) break;
		}

		if (theta_grid_size > 1 && theta_step != 0.0) printf("\n");
		if (n == 0 && R_step == 0.0) break;
	}

	return EXIT_SUCCESS;
}
