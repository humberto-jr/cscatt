#include "modules/math.h"
#include "modules/file.h"
#include "modules/matrix.h"
#include "modules/globals.h"

int main(int argc, char *argv[])
{
	if (argc != 5)
	{
		PRINT_ERROR("%d arguments given. Usage: ./simpson.out [a] [b] [col] [file]\n", argc - 1)
		return EXIT_FAILURE;
	}

	ASSERT(argv != NULL)

	const double a = atof(argv[1]);

	const double b = atof(argv[2]);

	const double col = atoi(argv[3]);

	ASSERT(col > -1)

	FILE *input = fopen(argv[4], "r");

	if (input == NULL)
	{
		PRINT_ERROR("unable to open %s\n", argv[4])
		exit(EXIT_FAILURE);
	}

	const int max_row = file_row_count(input);

	const int max_col = file_col_count(input);

	if (col >= max_col)
	{
		PRINT_ERROR("there are only %d columns in %s\n", max_col, argv[4])
		exit(EXIT_FAILURE);
	}

	matrix *data = matrix_read(input, max_row, max_col);

	fclose(input);

	const int n_max = matrix_row(data);

	double *integrand = matrix_raw_col(data, col, false);

	matrix_free(data);

	printf("% -8e\n", math_simpson_array(a, b, n_max, false, integrand));

	free(integrand);

	return EXIT_SUCCESS;
}
