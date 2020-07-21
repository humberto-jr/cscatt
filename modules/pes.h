#if !defined(PES_HEADER)
	#define PES_HEADER
	#include "globals.h"
	#include "matrix.h"

	#define MAX_JACOBI_VECTOR 5
	#define MAX_ATOM (MAX_JACOBI_VECTOR + 1)
	#define MAX_INTERNUC_DISTANCE (3*MAX_ATOM - 6)

	void pes_init_mass(FILE *input, const char atom);

	double pes_mass(const char atom);

	double pes_mass_bc();

	double pes_mass_ac();

	double pes_mass_ab();

	double pes_mass_abc(const char arrang);

	double pes_abc(const char arrang,
	               const double r, const double R, const double theta);

	double pes_bc(const int j, const double r);

	double pes_ac(const int j, const double r);

	double pes_ab(const int j, const double r);

	double pes_legendre_multipole(const char arrang,
	                              const int lambda,
	                              const double r,
	                              const double R);

	void pes_init();

//	matrix *pes_olson_smith_model(const double x);

//	matrix *pes_tully_model(const int n, const double x);

	void pes_about(FILE *output);
#endif
