#if !defined(COUPL_CONFIG_HEADER)
	#define COUPL_CONFIG_HEADER

	#if !defined(CMATRIX_BUFFER_FORMAT)
		#define CMATRIX_BUFFER_FORMAT "cmatrix_arrang=%c_n=%d_J=%d.bin"
	#endif

	#if !defined(CMATRIX_DATAFILE_FORMAT)
		#define CMATRIX_DATAFILE_FORMAT "cmatrix_arrang=%c_n=%d_J=%d.dat"
	#endif

	/******************************************************************************

	 Function check_coupl(): checks whether the n-th coupling matrix file exist in
	 the disk for a given arrangement and J.

	******************************************************************************/

	inline static bool check_coupl(const char arrang, const int n, const int J)
	{
		char filename[MAX_LINE_LENGTH];
		sprintf(filename, CMATRIX_BUFFER_FORMAT, arrang, n, J);

		return file_exist(filename);
	}

	/******************************************************************************

	 Function count_coupl(): counts how many coupling matrices are available in the
	 disk for a given arrangement and J.

	******************************************************************************/

	inline static int count_coupl(const char arrang, const int J)
	{
		int counter = 0;
		while (check_coupl(arrang, counter, J)) ++counter;

		return counter;
	}

	/******************************************************************************

	 Function load_coupl(): load from the disk the n-th coupling matrix for a given
	 arrangement and J.

	******************************************************************************/

	inline static void load_coupl(const char arrang,
	                              const int n, const int J, matrix *c)
	{
		char filename[MAX_LINE_LENGTH];
		sprintf(filename, CMATRIX_BUFFER_FORMAT, arrang, n, J);

		c = matrix_load(filename);
	}
#endif
