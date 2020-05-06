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
		matrix *wavef;
		int v, j, l, spin_mult, max_state;
		double r_min, r_max, r_step, energy;
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
		while (check_basis_file(arrang, counter, J)) ++counter;

		return counter;
	}

	/******************************************************************************

	 Function load_basis(): load from the disk the basis for a given arrangement,
	 channel index and J

	******************************************************************************/

	void load_basis(const char arrang,
	                const int ch, const int J, scatt_basis *b)
	{
		ASSERT(b != NULL)

		char filename[MAX_LINE_LENGTH];
		sprintf(filename, BASIS_BUFFER_FORMAT, arrang, ch, J);

		b->wavef = matrix_load(filename);

		int file_offset = matrix_sizeof(b->wavef);
		file_read(filename, 1, sizeof(double), &b->r_min, file_offset);

		file_offset += sizeof(double);
		file_read(filename, 1, sizeof(double), &b->r_max, file_offset);

		file_offset += sizeof(double);
		file_read(filename, 1, sizeof(double), &b->r_step, file_offset);

		file_offset += sizeof(double);
		file_read(filename, 1, sizeof(double), &b->energy, file_offset);

		file_offset += sizeof(double);
		file_read(filename, 1, sizeof(int), &b->v, file_offset);

		file_offset += sizeof(int);
		file_read(filename, 1, sizeof(int), &b->j, file_offset);

		file_offset += sizeof(int);
		file_read(filename, 1, sizeof(int), &b->l, file_offset);

		file_offset += sizeof(int);
		file_read(filename, 1, sizeof(int), &b->spin_mult, file_offset);

		file_offset += sizeof(int);
		file_read(filename, 1, sizeof(int), &b->max_state, file_offset);
	}
#endif
