#include "modules/file.h"
#include "modules/spline.h"
#include "modules/matrix.h"
#include "modules/globals.h"

int main(int argc, char *argv[])
{
	if (argc != 7)
	{
		PRINT_ERROR("%d arguments given. Usage: ./spline.out [x column] [f columns, start] [f columns, end] [grid size] [type] [file]\n", argc - 1)
		return EXIT_FAILURE;
	}

	ASSERT(argv != NULL)

	const double x_col = atoi(argv[1]);

	const double col_start = atoi(argv[2]);

	const double col_end = atoi(argv[3]);

	const double grid_size = atoi(argv[4]);

	const char type = argv[5][0];

	ASSERT(x_col > -1)
	ASSERT(grid_size > 0)
	ASSERT(col_start > -1)
	ASSERT(col_end >= col_start)

	FILE *input = file_open(argv[6], "r");

	const int max_row = file_row_count(input);

	const int max_col = file_col_count(input);

	if (x_col >= max_col || col_end >= max_col)
	{
		PRINT_ERROR("there are only %d columns in %s\n", max_col, argv[6])
		exit(EXIT_FAILURE);
	}

	matrix *data = matrix_read(input, max_row, max_col);

	fclose(input);

	double *x = matrix_raw_col(data, x_col, false);

	const double grid_step = (x[max_row - 1] - x[0])/as_double(grid_size);

	spline_handle *h = allocate(col_end - col_start + 1, sizeof(spline_handle), true);

	int counter = 0;
	for (int col = col_start; col <= col_end; ++col)
	{
		h[counter].x = x;
		h[counter].f = matrix_raw_col(data, col, false);
		h[counter].s = spline_alloc(max_row, x, h[counter].f, type);
		++counter;
	}

	matrix_free(data);

	for (int n = 0; n < grid_size; ++n)
	{
		const double x_value = x[0] + as_double(n)*grid_step;

		printf("% -8e", x_value);

		for (int col = 0; col < counter; ++col)
			printf("\t % -8e", spline_value(h[col].s, x_value));

		printf("\n");
	}

	for (int col = 0; col < counter; ++col)
	{
		spline_free(h[col].s);
		free(h[col].f);
	}

	free(x);
	free(h);

	return EXIT_SUCCESS;
}
