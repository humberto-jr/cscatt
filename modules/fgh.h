#if !defined(FGH_HEADER)
	#define FGH_HEADER
	#include "globals.h"
	#include "mpi_lib.h"
	#include "matrix.h"

	/******************************************************************************

	 Type fgh_basis: represents one component of a scattering basis set based upon
	 asymptotic rovibrational FGH states (single channel, n = 0, or multichannel,
	 n > 0).

	******************************************************************************/

	struct fgh_basis
	{
		int v, j, l, n, grid_size;
		double r_min, r_max, r_step, eigenval, *eigenvec;
	};

	typedef struct fgh_basis fgh_basis;

	matrix *fgh_dense_single_channel(const int grid_size,
	                                 const double grid_step,
	                                 const double pot_energy[],
	                                 const double mass);

	matrix *fgh_dense_multi_channel(const int max_state,
	                                const int grid_size,
	                                const double grid_step,
	                                const tensor pot_energy[],
	                                const double mass);

	mpi_matrix *fgh_sparse_multi_channel(const int max_state,
	                                     const int grid_size,
	                                     const double grid_step,
	                                     const tensor pot_energy[],
	                                     const double mass);

	double *fgh_eigenvec(const matrix *fgh, const int v, const double grid_step);

	int fgh_basis_count(const char arrang, const int J);

	FILE *fgh_basis_file(const char arrang, const int n,
	                     const int J, const char mode[], const bool verbose);

	void fgh_basis_write(const fgh_basis *b, FILE *output);

	void fgh_basis_read(fgh_basis *b, FILE *input);
#endif
