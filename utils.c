#include "modules/file.h"
#include "modules/matrix.h"
#include "modules/globals.h"

#include "utils.h"

#if !defined(BASIS_FORMAT)
	#define BASIS_FORMAT "basis_arrang=%c_ch=%d_J=%d.%s"
#endif

#if !defined(MULTIPOLE_FORMAT)
	#define MULTIPOLE_FORMAT "multipole_arrang=%c_n=%d.%s"
#endif

#if !defined(COUPLING_MATRIX_FORMAT)
	#define COUPLING_MATRIX_FORMAT "cmatrix_arrang=%c_n=%d_J=%d.bin"
#endif

#if !defined(COUPLING_DATAFILE_FORMAT)
	#define COUPLING_DATAFILE_FORMAT "cmatrix_arrang=%c_n=%d_J=%d.dat"
#endif

/******************************************************************************

 Function basis_count(): counts how many basis functions are available in the
 disk for a given arrangement and total angular momentum, J.

******************************************************************************/

int basis_count(const char arrang, const int J)
{
	char filename[MAX_LINE_LENGTH];

	int counter = 0;
	sprintf(filename, BASIS_FORMAT, arrang, counter, J, "bin");

	while (file_exist(filename))
	{
		++counter;
		sprintf(filename, BASIS_FORMAT, arrang, counter, J, "bin");
	}

	return counter;
}

/******************************************************************************

 Function basis_file(): opens the file for the n-th channel, arrangement and
 total angular momentum, J. Where, mode is the file access mode of fopen()
 from the C library.

 NOTE: mode = "wb" for write + binary format and "rb" for read + binary format,
 extension used is .bin otherwise .dat is used assuming text mode ("w" or "r").

******************************************************************************/

FILE *basis_file(const char arrang, const int n, const int J, const char mode[])
{
	char filename[MAX_LINE_LENGTH], ext[4];

	if (strlen(mode) > 1 && mode[1] == 'b')
		sprintf(ext, "%s", "bin");
	else
		sprintf(ext, "%s", "dat");

	sprintf(filename, BASIS_FORMAT, arrang, n, J, ext);

	FILE *stream = fopen(filename, mode);

	ASSERT(stream != NULL)

//	if (stream != NULL && mode[0] == 'r') printf("# Reading %s\n", filename);
//	if (stream != NULL && mode[0] == 'w') printf("# Writing %s\n", filename);

	return stream;
}

/******************************************************************************

 Function basis_read(): loads from the disk the basis function for the n-th
 channel for a given arrangement and total angular momentum, J.

******************************************************************************/

void basis_read(const char arrang, const int n, const int J, basis *b)
{
	ASSERT(b != NULL)

	FILE *input = basis_file(arrang, n, J, "rb");

	file_read(&b->v, sizeof(int), 1, input, 0);

	file_read(&b->j, sizeof(int), 1, input, 0);

	file_read(&b->l, sizeof(int), 1, input, 0);

	file_read(&b->n, sizeof(int), 1, input, 0);

	file_read(&b->r_min, sizeof(double), 1, input, 0);

	file_read(&b->r_max, sizeof(double), 1, input, 0);

	file_read(&b->r_step, sizeof(double), 1, input, 0);

	file_read(&b->eigenval, sizeof(double), 1, input, 0);

	file_read(&b->grid_size, sizeof(int), 1, input, 0);

	b->eigenvec = allocate(b->grid_size, sizeof(double), false);

	file_read(b->eigenvec, sizeof(double), b->grid_size, input, 0);

	fclose(input);
}

/******************************************************************************

 Function multipole_count(): counts how many multipole files (one per grid
 point) are available in the disk for a given arrangement.

******************************************************************************/

int multipole_count(const char arrang)
{
	char filename[MAX_LINE_LENGTH];

	int counter = 0;
	sprintf(filename, MULTIPOLE_FORMAT, arrang, counter, "bin");

	while (file_exist(filename))
	{
		++counter;
		sprintf(filename, MULTIPOLE_FORMAT, arrang, counter, "bin");
	}

	return counter;
}

/******************************************************************************

 Function multipole_file(): opens the multipole file for a given grid point
 index, arrangement and total angular momentum, J. Where, mode is the file
 access mode of fopen() from the C library.

 NOTE: mode = "wb" for write + binary format and "rb" for read + binary format,
 extension used is .bin otherwise .dat is used assuming text mode ("w" or "r").

******************************************************************************/

FILE *multipole_file(const char arrang, const int grid_index, const char mode[])
{
	char filename[MAX_LINE_LENGTH], ext[4];

	if (strlen(mode) > 1 && mode[1] == 'b')
		sprintf(ext, "%s", "bin");
	else
		sprintf(ext, "%s", "dat");

	sprintf(filename, MULTIPOLE_FORMAT, arrang, grid_index, ext);

	FILE *stream = fopen(filename, mode);

	ASSERT(stream != NULL)

//	if (stream != NULL && mode[0] == 'r') printf("# Reading %s\n", filename);
//	if (stream != NULL && mode[0] == 'w') printf("# Writing %s\n", filename);

	return stream;
}

/******************************************************************************

 Function multipole_read(): loads from the disk a set of multipole terms for a
 given grid point index, arrangement and total angular momentum, J.

******************************************************************************/

void multipole_read(const char arrang, const int n, multipole_set *m)
{
	ASSERT(m != NULL)

	FILE *input = multipole_file(arrang, n, "rb");

	file_read(&m->R, sizeof(double), 1, input, 0);

	file_read(&m->r_min, sizeof(double), 1, input, 0);

	file_read(&m->r_max, sizeof(double), 1, input, 0);

	file_read(&m->r_step, sizeof(double), 1, input, 0);

	file_read(&m->lambda_max, sizeof(int), 1, input, 0);

	file_read(&m->grid_size, sizeof(int), 1, input, 0);

	m->set = allocate(m->lambda_max, sizeof(multipole), true);

	int lambda = -1;
	while (lambda != m->lambda_max)
	{
		file_read(&lambda, sizeof(int), 1, input, 0);

		m->set[lambda].value = allocate(m->grid_size, sizeof(double), false);

		file_read(m->set[lambda].value, sizeof(double), m->grid_size, input, 0);
	}

	fclose(input);
}

void multipole_free(multipole_set *m)
{
	if (m->set == NULL) return;

	for (int lambda = 0; lambda <= m->lambda_max; ++lambda)
	{
		if (m->set[lambda].value != NULL) free(m->set[lambda].value);
	}

	free(m->set);
}

/******************************************************************************

 Function coupling_count(): counts how many coupling matrices are available in
 the disk for a given arrangement and total angular momentum, J.

******************************************************************************/

int coupling_count(const char arrang, const int J)
{
	char filename[MAX_LINE_LENGTH];

	int counter = 0;
	sprintf(filename, COUPLING_MATRIX_FORMAT, arrang, counter, J);

	while (file_exist(filename))
	{
		++counter;
		sprintf(filename, COUPLING_MATRIX_FORMAT, arrang, counter, J);
	}

	return counter;
}

/******************************************************************************

 Function coupling_write(): save in the disk a set of multipole terms for a
 given grid point index, arrangement and total angular momentum, J.

******************************************************************************/

void coupling_write(const char arrang, const int grid_index,
                    const int J, const bool verbose, const matrix *c)
{
	char filename[MAX_LINE_LENGTH];
	sprintf(filename, COUPLING_MATRIX_FORMAT, arrang, grid_index, J);

	if (verbose) printf("# Writing %s\n", filename);

	matrix_save(c, filename);
}

/******************************************************************************

 Function coupling_read(): load from the disk a set of multipole terms for a
 given grid point index, arrangement and total angular momentum, J.

******************************************************************************/

matrix *coupling_read(const char arrang,
                      const int grid_index, const int J, const bool verbose)
{
	char filename[MAX_LINE_LENGTH];
	sprintf(filename, COUPLING_MATRIX_FORMAT, arrang, grid_index, J);

	if (verbose) printf("# Reading %s\n", filename);

	return matrix_load(filename);
}

/******************************************************************************

 Function coupling_datafile(): opens a datafile in the disk for the n-th scatt.
 channel's coupling matrix element, arrangement and total angular momentum, J.

 NOTE: mode = "w" for write + text format and "r" for read + text format, with
 extension .dat used.

******************************************************************************/

FILE *coupling_datafile(const char arrang, const int n,
                        const int J, const bool verbose, const char mode[])
{
	char filename[MAX_LINE_LENGTH];
	sprintf(filename, COUPLING_DATAFILE_FORMAT, arrang, n, J);

	if (verbose) printf("# Creating %s\n", filename);

	FILE *stream = fopen(filename, mode);

	ASSERT(stream != NULL)

	return stream;
}
