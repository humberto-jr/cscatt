#include "modules/pes.h"
#include "modules/file.h"
#include "modules/matrix.h"
#include "modules/mpi_lib.h"
#include "modules/globals.h"

/******************************************************************************

 Function fgh_multich_matrix(): the same as fgh_matrix(), except that the
 potential energy is the one of a problem with max_ch channels within grid_size
 points. For details see Eq. (19a), Eq. (19b) and Eq. (19c) of Ref. [2].

******************************************************************************/

void print(const int max_ch,
           const int grid_size,
           const double grid_step,
           const tensor pot_energy[],
           const double mass)
{
	ASSERT(max_ch > 0)
	ASSERT(pot_energy != NULL)

	const int size = grid_size*max_ch;

	const int non_zeros[2] = {1, 1};

	mpi_matrix *result = mpi_matrix_alloc(size, size, non_zeros);

	const double box_length = as_double(grid_size - 1)*grid_step;

	const double factor = (M_PI*M_PI)/(mass*box_length*box_length);

	const double nn_term = factor*as_double(grid_size*grid_size + 2)/6.0;

	int row = 0;
	for (int p = 0; p < max_ch; ++p)
	{
		for (int n = 0; n < grid_size; ++n)
		{
			int col = 0;
			for (int q = 0; q < max_ch; ++q)
			{
				for (int m = 0; m < grid_size; ++m)
				{
					double pnqm = 0.0;

					if (n == m && p == q)
					{
						pnqm = nn_term
						     + matrix_get(pot_energy[n].value, p, p);
					}

					if (n != m && p == q)
					{
						const double nm = as_double((n + 1) - (m + 1));

						const double nm_term = sin(nm*M_PI/as_double(grid_size));

						pnqm = pow(-1.0, nm)*factor/pow(nm_term, 2);
					}

					if (n == m && p != q)
					{
						pnqm = matrix_get(pot_energy[n].value, p, q);
					}

					if (pnqm != 0.0) printf("%d\t %d\t \n", row, col, pnqm);
					++col;
				}
			}

			++row;
		}
	}
}

int main(int argc, char *argv[])
{
	if (argc != 6)
	{
		PRINT_ERROR("%d arguments given. Usage: ./3j.out [arrang] [J] [R_min] [R_max] [mass]\n", argc - 1)
		return EXIT_FAILURE;
 	}

	ASSERT(argv != NULL)

	const char arrang = argv[1][0];

	const int J = atoi(argv[2]);

	const double R_min = atof(argv[3]);

	const double R_max = atof(argv[4]);

	const double mass = atof(argv[5]);

	const int n_max = fgh_basis_count(arrang, J);

	const double R_step = (R_max - R_min)/as_double(n_max);

	tensor *pot_energy = allocate(n_max, sizeof(tensor), true);

	for (int n = 0; n < n_max; ++n)
	{
		FILE *input = fgh_basis_file(arrang, n, J, "rb", false);

		fgh_basis_read(&pot_energy[n].value, input);

		fclose(input);
	}

	const int max_state = matrix_row(pot_energy[0].value);

	print(max_state, n_max, R_step, pot_energy, mass);

	free(pot_energy);
	return EXIT_SUCCESS;
}
