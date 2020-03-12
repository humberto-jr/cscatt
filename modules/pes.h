#if !defined(PES_HEADER)
	#define PES_HEADER
	#include "globals.h"
	#include "jacobi.h"
	#include "coor.h"

	double pes(const jacobi_coor *x);

	double pes_sf(const jacobi_sf *set);

	void pes_init();

	double pec(const char arrang, const double r);

	double pes_asymptotic_min(const char arrang, const double scan_step);

	void pes_about(FILE *output);
#endif
