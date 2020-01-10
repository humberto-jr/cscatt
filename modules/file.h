#if !defined(FILE_HEADER)
	#define FILE_HEADER
	#include "clib.h"

	bool file_exist(const char filename[]);

	FILE *file_open_input(const char filename[], const bool bin_format);

	FILE *file_open_output(const char filename[],
	                       const bool bin_format, const bool to_append);

	void file_init_stdin(const char filename[]);

	char *file_find_string(FILE *input, const char pattern[]);

	double file_get_key(FILE *input, const char key[], const double min,
	                    const double max, const double default_value);

	int file_row_count(FILE *input);

	int file_col_count(FILE *input);

	bool file_end(FILE *input);

	void file_write(const char filename[], const int n,
	                const int data_size, const void *data, const bool to_append);

	void file_read(const char filename[], const int n,
	               const int data_size, void *data, const int offset);
#endif
