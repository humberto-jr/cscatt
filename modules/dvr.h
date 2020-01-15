#if !defined(DVR_HEADER)
	#define DVR_HEADER
	#include "globals.h"
	#include "matrix.h"

	matrix *dvr_fgh(const int grid_size,
	                const double grid_step,
	                const double pot_energy[],
	                const double mass);

	double dvr_fgh_wavef(const matrix *fgh,
	                     const int a,
	                     const double r_min,
	                     const double r_max,
	                     const double r_new);

	double dvr_fgh_product(const matrix *fgh,
	                       const int a,
	                       const int b,
	                       const double r_min,
	                       const double r_max,
	                       const double r_new);

	void dvr_fgh_norm(matrix *fgh,
	                  const int a,
	                  const double grid_step,
	                  const bool use_omp);

	matrix *dvr_multich_fgh(const int max_ch,
	                        const int grid_size,
	                        const double grid_step,
	                        const tensor pot_energy[],
	                        const double mass);

	matrix *dvr_multich_fgh_comp(matrix *fgh,
	                             const int max_ch,
	                             const int n,
	                             const int v);

	void dvr_multich_fgh_norm(matrix *fgh,
	                          const int max_ch,
	                          const int v,
	                          const double grid_step,
	                          const bool use_omp);
#endif
