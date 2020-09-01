#if !defined(UTILS_HEADER)
	#define UTILS_HEADER

	/******************************************************************************

	 Type basis: represent one component of a scattering basis set based upon
	 asymptotic diatomic rovibrational states.

	******************************************************************************/

	struct basis
	{
		int v, j, l, n, grid_size;
		double r_min, r_max, r_step, eigenval, *eigenvec;
	};

	typedef struct basis basis;

	/******************************************************************************

	 Type basis: represent one component of a scattering basis set based upon
	 asymptotic diatomic rovibrational states.

	******************************************************************************/

	struct multipole
	{
		double *value;
	};

	typedef struct multipole multipole;

	struct multipole_set
	{
		multipole *set;
		int lambda_max, grid_size;
		double R, r_min, r_max, r_step;
	};

	typedef struct multipole_set multipole_set;

	int basis_count(const char arrang, const int J);

	FILE *basis_file(const char arrang,
	                 const int n, const int J, const char mode[]);

	void basis_read(const char arrang, const int n, const int J, basis *b);

	int multipole_count(const char arrang);

	FILE *multipole_file(const char arrang,
	                     const int grid_index, const char mode[]);

	void multipole_read(const char arrang, const int n, multipole_set *m);

	void multipole_free(multipole_set *m);

	int coupling_count(const char arrang, const int J);

	void coupling_write(const char arrang, const int grid_index,
	                    const int J, const bool verbose, const matrix *c);

	matrix *coupling_read(const char arrang,
	                      const int grid_index, const int J, const bool verbose);

	FILE *coupling_datafile(const char arrang, const int n,
	                        const int J, const bool verbose, const char mode[]);

	void ratio_write(const char arrang, const int m,
	                 const int J, const bool verbose, const matrix *ratio);

	matrix *ratio_read(const char arrang,
	                   const int m, const int J, const bool verbose);
#endif
