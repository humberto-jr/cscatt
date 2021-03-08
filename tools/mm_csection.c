#include "modules/math.h"
#include "modules/file.h"
#include "modules/globals.h"

int main(int argc, char *argv[])
{
	if (argc != 5)
	{
		PRINT_ERROR("%d arguments given. Usage: %s [j, in] [m, in] [beta] [filename]\n", argc - 1, argv[0])
		return EXIT_FAILURE;
	}

	const double j_in = atoi(argv[1]);
	const double m_in = atoi(argv[2]);
	const double beta = atof(argv[3]);
	const int j = as_int(j_in);

	FILE *input = file_open(argv[4], "r");

	char line[MAX_LINE_LENGTH];
	while (fgets(line, sizeof(line), input) != NULL)
	{
		if ((line[0] != '#') && (line[0] != '\n') && (line[0] != '\0'))
		{
			char *token = strtok(line, " \t");

			if (token == NULL)
			{
				PRINT_ERROR("invalid entry at line '%s' (1)\n", line)
				exit(EXIT_FAILURE);
			}

			const double coll_energy = atof(token);

			token = strtok(NULL, " \t");

			if (token == NULL)
			{
				PRINT_ERROR("invalid entry at line '%s' (2)\n", line)
				exit(EXIT_FAILURE);
			}

			const double sigma = atof(token);

			double sum = 0.0;
			for (int m = 0; m <= j; ++m)
			{
				double *d = math_wigner_d(m_in, as_double(m), j_in, beta);

				if (m == 0)
					sum += sigma*pow(d[j], 2);
				else
					sum += sigma*pow(d[j], 2)*2.0;

				free(d);
			}

			printf("% -8e\t % -8e\t % -8e\n", coll_energy*(1.160451812E4), sigma, sum);
		}
	}

	fclose(input);
	return EXIT_SUCCESS;
}
