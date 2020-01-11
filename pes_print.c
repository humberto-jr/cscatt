#include "modules/pes.h"
#include "modules/coor.h"
#include "modules/file.h"
#include "modules/globals.h"

int main(int argc, char *argv[])
{
	ASSERT(argc > 1)
	file_init_stdin(argv[1]);

/*
 *	Vibrational grid, r:
 */

	const int rovib_grid_size
		= (int) file_get_key(stdin, "rovib_grid_size", 1.0, INF, 100.0);

	const double r_min
		= file_get_key(stdin, "r_min", 0.0, INF, 0.0);

	const double r_max
		= file_get_key(stdin, "r_max", r_min, INF, r_min + 100.0);

	const double r_step
		= (r_max - r_min)/as_double(rovib_grid_size);

/*
 *	Scattering grid, R:
 */

	const int scatt_grid_size
		= (int) file_get_key(stdin, "scatt_grid_size", 1.0, INF, 100.0);

	const double R_min
		= file_get_key(stdin, "R_min", 0.0, INF, 0.0);

	const double R_max
		= file_get_key(stdin, "R_max", R_min, INF, R_min + 100.0);

	const double R_step
		= (R_max - R_min)/as_double(scatt_grid_size);

/*
 *	Angular grid, theta:
 */

	const int theta_grid_size
		= (int) file_get_key(stdin, "theta_grid_size", 1.0, INF, 36.0);

	const double theta_min
		= file_get_key(stdin, "theta_min", 0.0, 180.0, 0.0);

	const double theta_max
		= file_get_key(stdin, "theta_max", theta_min, 180.0, 180.0);

	const double theta_step
		= (theta_max - theta_min)/as_double(theta_grid_size);

/*
 *	Energy scale and shifts:
 */

	const double shift
		= file_get_key(stdin, "energy_shift", -INFINITY, INFINITY, 0.0);

	const double scale
		= file_get_key(stdin, "energy_scale", -INFINITY, INFINITY, 1.0);

/*
 *	Arrangement:
 */

	const char arrang
		= 96 + (int) file_get_key(stdin, "arrang", 1.0, 3.0, 1.0);

	printf("# r            = (%f, %f, %f)\n", r_min, r_max, r_step);
	printf("# R            = (%f, %f, %f)\n", R_min, R_max, R_step);
	printf("# theta        = (%f, %f, %f)\n", theta_min, theta_max, theta_step);
	printf("# Arrangement  = %c\n", arrang);
	printf("# Energy shift = %f\n", shift);
	printf("# Energy scale = %f\n", scale);

	printf("# \n");

	jacobi_coor coor;
	coor.arrang = arrang;

	for (int n = 0; n < rovib_grid_size; ++n)
	{
		coor.r = r_min + as_double(n)*r_step;

		for (int m = 0; m < scatt_grid_size; ++m)
		{
			coor.R = R_min + as_double(m)*R_step;

			for (int p = 0; p <= theta_grid_size; ++p)
			{
				coor.theta = theta_min + as_double(p)*theta_step;

				printf("% 6f\t % 6f\t % 6f\t % -8e\n",
				coor.r, coor.R, coor.theta, (pes(&coor) + shift)*scale);
			}

			if (rovib_grid_size > 1 && scatt_grid_size > 1) printf("\n");
		}

		if (scatt_grid_size > 1) printf("\n");
	}

	return EXIT_SUCCESS;
}
