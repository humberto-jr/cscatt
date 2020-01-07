#if !defined(PHYS_HEADER)
	#define PHYS_HEADER

	double phys_wigner_3j(const int a, const int b, const int c,
	                      const int d, const int e, const int f);

	double phys_wigner_6j(const int a, const int b, const int c,
	                      const int d, const int e, const int f);

	double phys_percival_seaton(const int spin_mult, const int j1, const int j2,
	                            const int l1, const int l2, const int lambda, const int J);
#endif
