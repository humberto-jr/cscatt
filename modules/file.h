#if !defined(FILE_HEADER)
	#define FILE_HEADER
	#include "globals.h"

	bool file_exist(const char filename[]);

	FILE *file_open_input(const char filename[], const bool bin_format);

	FILE *file_open_output(const char filename[],
	                       const bool bin_format, const bool to_append);

	void file_init_stdin(const char filename[]);

	void file_init_stdout(const char filename[], const bool to_append);

	char *file_find(FILE *input, const char pattern[]);

	double file_get_key(FILE *input, const char key[], const double min,
	                    const double max, const double default_value);

	double file_keyword(FILE *input, const char key[], const double min,
	                    const double max, const double default_value);

	int file_row_count(FILE *input);

	int file_col_count(FILE *input);

	bool file_end(FILE *stream);

	void file_write(void *data, const int data_size, const int n, FILE *stream);

	void file_read(void *data, const int data_size,
	               const int n, FILE *stream, const int offset);
#endif
