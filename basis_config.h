#if !defined(BASIS_CONFIG_HEADER)
	#define BASIS_CONFIG_HEADER

	#if !defined(BASIS_BUFFER_FORMAT)
		#define BASIS_BUFFER_FORMAT "basis_arrang=%c_ch=%d_J=%d.bin"
	#endif

	#if !defined(BASIS_DATAFILE_FORMAT)
		#define BASIS_DATAFILE_FORMAT "basis_arrang=%c_ch=%d_J=%d.dat"
	#endif

	/******************************************************************************

	 Type scatt_basis: represent one component of a scattering basis set based upon
	 asymptotic diatomic rovibrational states.

	******************************************************************************/

	struct scatt_basis
	{
		int v, j, l, grid_size, state;
		double r_min, r_max, grid_step, energy, *wavef;
	};

	typedef struct scatt_basis scatt_basis;

	/******************************************************************************

	 Function open_basis_file(): opens the basis function file for a given arrang.,
	 channel index and J. Where, mode is the file access mode of fopen() from the C
	 library.

	 NOTE: mode = "wb" for write + binary format and "rb" for read + binary format.

	******************************************************************************/

	inline static FILE *open_basis_file(const char mode[],
                                       const char arrang, const int ch, const int J)
	{
		char filename[MAX_LINE_LENGTH];
		sprintf(filename, BASIS_BUFFER_FORMAT, arrang, ch, J);

		FILE *basis = fopen(filename, mode);

		if (basis == NULL)
		{
			PRINT_ERROR("unable to open %s\n", filename)
			exit(EXIT_FAILURE);
		}

		return basis;
	}

	/******************************************************************************

	 Function check_basis_file(): checks whether a basis function file exist in the
	 disk for a given arrang., channel index and J.

	******************************************************************************/

	inline static bool check_basis_file(const char arrang,
	                                    const int ch, const int J)
	{
		char filename[MAX_LINE_LENGTH];
		sprintf(filename, BASIS_BUFFER_FORMAT, arrang, ch, J);

		return file_exist(filename);
	}

	/******************************************************************************

	 Function count_basis(): counts how many basis functions are available in the
	 disk for a given arrangement and J.

	******************************************************************************/

	inline static int count_basis_file(const char arrang, const int J)
	{
		int counter = 0;
		while (check_basis_file(arrang, counter, J) == true) ++counter;

		return counter;
	}

	/******************************************************************************

	 Function read_basis_file(): load from the disk the basis for a given arrang.,
	 channel index and J

	******************************************************************************/

	void read_basis_file(const char arrang,
	                     const int ch, const int J, scatt_basis *b)
	{
		ASSERT(b != NULL)

		FILE *input = open_basis_file("rb", arrang, ch, J);

		int status;

		status = fread(&b->v, sizeof(int), 1, input);
		status = fread(&b->j, sizeof(int), 1, input);
		status = fread(&b->l, sizeof(int), 1, input);
		status = fread(&b->state, sizeof(int), 1, input);

		status = fread(&b->r_min, sizeof(double), 1, input);
		status = fread(&b->r_max, sizeof(double), 1, input);
		status = fread(&b->grid_step, sizeof(double), 1, input);

		status = fread(&b->energy, sizeof(double), 1, input);
		status = fread(&b->grid_size, sizeof(int), 1, input);

		b->wavef = allocate(b->grid_size, sizeof(double), false);

		status = fread(b->wavef, sizeof(double), b->grid_size, input);

		ASSERT(status == b->grid_size)
		fclose(input);
	}
#endif
