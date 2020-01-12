#if !defined(MASS_CONFIG_HEADER)
	#define MASS_CONFIG_HEADER

	/******************************************************************************

	 Function init_lambda_step(): returns the step for the lambda expansion index
	 for a given arrangement label using the following notation: 'a' (A + BC), 'b'
	 (B + CA) and 'c' (C + AB). Where, A, B and C represents three atoms.

	******************************************************************************/

	inline static int init_lambda_step(const char arrang)
	{
		switch (arrang)
		{
			case 'a': return (mass(atom_b) == mass(atom_c)? 2 : 1);
			case 'b': return (mass(atom_c) == mass(atom_a)? 2 : 1);
			case 'c': return (mass(atom_a) == mass(atom_b)? 2 : 1);
		}

		return 0;
	}

	/******************************************************************************

	 Function init_atomic_masses(): initialize atomic masses from an input file and
	 returns the enum mass_case value used by mass() from the mass module.

	 NOTE: mode = 'p' for pair_bc, pair_ac or pair_ab (reduced diatomic masses)
	 depending on the kind of arrangement. Otherwise, arrang_a, arrang_b or
	 arrang_c (arrangement reduced masses).

	******************************************************************************/

	mass_case init_atomic_masses(FILE *input,
	                             const char arrang, const char mode)
	{
		mass_init(input);

		printf("#\n");
		printf("# REDUCED MASSES:\n");
		printf("# Atom A = %f a.u.\n", mass(atom_a));
		printf("# Atom B = %f a.u.\n", mass(atom_b));
		printf("# Atom C = %f a.u.\n", mass(atom_c));

		switch (arrang)
		{
			case 'a':
				printf("# Diatom BC = %f a.u.\n", mass(pair_bc));
				printf("# Arrangement A + BC = %f a.u.\n", mass(arrang_a));

				return (mode == 'p'? pair_bc : arrang_a);

			case 'b':
				printf("# Diatom CA = %f a.u.\n", mass(pair_ac));
				printf("# Arrangement B + CA = %f a.u.\n", mass(arrang_b));

				return (mode == 'p'? pair_ac : arrang_b);

			case 'c':
				printf("# Diatom AB = %f a.u.\n", mass(pair_ab));
				printf("# Arrangement C + AB = %f a.u.\n", mass(arrang_c));

				return (mode == 'p'? pair_ab : arrang_c);

			default:
				PRINT_ERROR("invalid arrangement %c\n", arrang)
				exit(EXIT_FAILURE);
		}
	}
#endif
