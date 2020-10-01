#if !defined(PES_HEADER)
	#define PES_HEADER
	#include "globals.h"
	#include "matrix.h"

	#define MAX_JACOBI_VECTOR 5
	#define MAX_ATOM (MAX_JACOBI_VECTOR + 1)
	#define MAX_INTERNUC_DISTANCE (3*MAX_ATOM - 6)

	/******************************************************************************

	 Type pes_multipole_set: value[lambda][n], n from r1

	******************************************************************************/

	struct pes_multipole
	{
		double R, r_min, r_max, r_step, **value;
		int lambda_min, lambda_max, lambda_step, grid_size;
	};

	typedef struct pes_multipole pes_multipole;

	typedef struct pes_multipole_set pes_multipole_set;

	void pes_set_inf(const double x_inf);

	void pes_init_mass(FILE *input, const char atom);

	void pes_init();

	double pes_mass(const char atom);

	double pes_mass_bc();

	double pes_mass_ac();

	double pes_mass_ab();

	double pes_mass_abc(const char arrang);

	double pes_mass_abcd();

	double pes_abc(const char arrang,
	               const double r, const double R, const double theta);

	double pes_abcd(const double r_bc,
	                const double r_bcd,
	                const double r_abcd,
	                const double theta_bc,
	                const double theta_a,
	                const double phi_a);

	double pes_bc(const int j, const double r);

	double pes_ac(const int j, const double r);

	double pes_ab(const int j, const double r);

	double pes_legendre_multipole(const char arrang,
	                              const int lambda,
	                              const double r,
	                              const double R);

	double pes_harmonics_multipole(const int eta,
	                               const int m_eta,
	                               const double r_bc,
	                               const double r_bcd,
	                               const double r_abcd,
	                               const double theta_bc);

	FILE *pes_multipole_file(const char arrang,
	                         const int n, const char mode[], const bool verbose);

	void pes_multipole_init(pes_multipole *m);

	void pes_multipole_init_all(const int n_max, pes_multipole m[]);

	void pes_multipole_write(const pes_multipole *m, FILE *output);

	void pes_multipole_write_all(const int n_max,
	                             const pes_multipole m[], FILE *output);

	void pes_multipole_read(pes_multipole *m, FILE *input);

	void pes_multipole_read_all(const int n_max, pes_multipole m[], FILE *input);

	void pes_multipole_free(pes_multipole *m);

	void pes_about(FILE *output);
#endif
