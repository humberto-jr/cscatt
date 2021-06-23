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
		size_t v, j, l, n, grid_size;
		double r_min, r_max, r_step, eigenval, *eigenvec;
	};

	typedef struct fgh_basis fgh_basis;

	matrix *fgh_dense_single_channel(const size_t grid_size,
	                                 const double grid_step,
	                                 const double pot_energy[],
	                                 const double mass);

	matrix *fgh_dense_multi_channel(const size_t max_state,
	                                const size_t grid_size,
	                                const double grid_step,
	                                const tensor pot_energy[],
	                                const double mass);

	mpi_matrix *fgh_sparse_multi_channel(const size_t max_state,
	                                     const size_t grid_size,
	                                     const double grid_step,
	                                     const tensor pot_energy[],
	                                     const double mass);

	double fgh_interpolation(const size_t grid_size,
	                         const double eigenvec[],
	                         const double r_min,
	                         const double r_max,
	                         const double r_new);

	double *fgh_eigenvec(const matrix *fgh,
	                     const size_t v, const double grid_step);

	double *fgh_multi_channel_eigenvec(const matrix *fgh,
	                                   const double grid_step,
	                                   const size_t max_state,
	                                   const size_t v,
	                                   const size_t c);

	size_t fgh_basis_count(const char dir[], const char arrang, const size_t J);

	FILE *fgh_basis_file(const char dir[], const char arrang, const size_t n,
	                     const size_t J, const char mode[], const bool verbose);

	void fgh_basis_write(const fgh_basis *b, FILE *output);

	void fgh_basis_read(fgh_basis *b, FILE *input);
#endif
