#if !defined(JOHNSON_HEADER)
	#define JOHNSON_HEADER
	#include "matrix.h"

	double johnson_riccati_bessel(const char type,
	                              const int l,
	                              const double wavenum,
	                              const double x);

	double johnson_modif_spher_bessel(const char type,
	                                  const int l,
	                                  const double wavenum,
	                                  const double x);

	double *johnson_jcp77_numerov(const int grid_size,
	                              const double grid_step,
	                              const double pot_energy[],
	                              const double trial_energy,
	                              const double mass,
	                              double *error,
	                              int *nodes);

	void johnson_jcp78_numerov(const double grid_step,
	                           matrix *pot_energy,
	                           const double tot_energy,
	                           const double mass,
	                           matrix *ratio);

	void johnson_jcp73_logd(const int n,
	                        const int grid_size,
	                        const double grid_step,
	                        const matrix *pot_energy,
	                        matrix *y);

	matrix *johnson_kmatrix(const int l[],
	                        const double grid_step,
	                        const double tot_energy,
	                        const double mass,
	                        const double level[],
	                        const matrix *ratio,
	                        const double R);

	matrix *johnson_smatrix(const matrix *k);
#endif
