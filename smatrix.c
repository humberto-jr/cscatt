#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <stdbool.h>

#include "config.h"
#include "scatt_utils.h"

int main(int argc, char *argv[])
{
	init_stdin(argc, argv);

	const int J
		= (int) get_key(stdin, "J", 0.0, INFINITY, 0.0);

	const int coll_grid_size
		= (int) get_key(stdin, "coll_grid_size", 1.0, INFINITY, 100.0);

	const int scatt_grid_size
		= (int) get_key(stdin, "scatt_grid_size", 1.0, INFINITY, 500.0);

	const double R_min
		= get_key(stdin, "R_min", 0.0, INFINITY, 0.0);

	const double R_max
		= get_key(stdin, "R_max", R_min, INFINITY, R_min + 100.0);

	const double R_inc
		= (R_max - R_min)/as_double(scatt_grid_size);

	const double R /* NOTE: the last grid point */
		= R_min + as_double(scatt_grid_size - 1)*R_inc;

	const double E_min
		= get_key(stdin, "E_min", -INFINITY, INFINITY, 0.0);

	const double E_max
		= get_key(stdin, "E_max", E_min, INFINITY, E_min);

	const double E_inc
		= (E_max - E_min)/as_double(coll_grid_size);

	const char arrang /* ASCII: 96 + 1 = a, 96 + 2 = b, 96 + 3 = c */
		= 96 + (int) get_key(stdin, "arrangement", 1.0, 3.0, 1.0);

	const double rmass
		= init_arrang_mass(arrang);

	printf("#\n");
	printf("# BUILD:\n");
	printf("# %s: %s @%s\n", __FILE__, __DATE__, __TIME__);

	int *l = NULL;
	double *inf_energy = NULL;

	for (int m = 0; m < coll_grid_size; ++m)
	{
		const double tot_energy = E_min + as_double(m)*E_inc;

		matrix r =
		{
			.data = NULL,
			.max_row = 0,
			.max_col = 0
		};

/*
 *		Step 1: Load the Numerov ratio matrix, r, for the last R-grid point:
 */

		char filename[MAX_LINE_LENGTH];
		sprintf(filename, RMATRIX_BUFFER_FORMAT, arrang, scatt_grid_size - 1, m, J);

		load_matrix(filename, &r);
		printf("# %5d: reading %s, size = %dx%d", m, filename, r.max_row, r.max_col);

/*
 *		Step 2: If within the first energy-grid point, load the asymptotic
 *		energies for each channel as well as their angular momentum l:
 */

		if (m == 0)
		{
			inf_energy = malloc(sizeof(double)*r.max_row);
			ASSERT(inf_energy != NULL);

			l = malloc(sizeof(int)*r.max_row);
			ASSERT(l != NULL);

			for (int n = 0; n < r.max_row; ++n)
			{
				scatt_basis basis;
				basis.wavef = NULL;

				memset(filename, 0, sizeof(filename));
				sprintf(filename, BASIS_BUFFER_FORMAT, arrang, n, J);

				load_ch(filename, &basis);

				l[n] = basis.l;
				inf_energy[n] = basis.energy;

				free(basis.wavef);
			}
		}

/*
 *		Step 3: Resolve the reactance matrix, k, and its open-open block, k_oo:
 */

		double *k
			= react_matrix(r.max_row, l, rmass, R_inc, tot_energy, inf_energy, r.data, R);
/*
		int open_ch = 0;

		double *k_oo
			= open_open_block(r.max_row, tot_energy, inf_energy, k, &open_ch);

		printf(", open channels = %d\n", open_ch);
*/
/*
 *		Step 4: Resolve the scattering matrix, s, and write it in the disk:
 */

		s_matrix *s = scatt_matrix(r.max_row, inf_energy, k, tot_energy);

		memset(filename, 0, sizeof(filename));
		sprintf(filename, SMATRIX_BUFFER_FORMAT, arrang, m, J);

		save_scatt_matrix(filename, tot_energy, s);

		free(s->re_part);
		free(s->im_part);
//		free(k_oo);
		free(k);
		free(s);
	}

	free(l);
	free(inf_energy);
	return EXIT_SUCCESS;
}
