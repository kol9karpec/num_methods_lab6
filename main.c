#include <stdio.h>
#include <stdbool.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

/*
 * Reads matrix @m x @n from @file
 *
 * @return 0 if succeded, -1 otherwise
 */
int read_matrix_from_file(const char * file, gsl_matrix * matrix,
		const unsigned int m, const unsigned int n);

/*
 * Finds max/min eigenvalues and appropriate eigenvectors by power method.
 *
 * @res_val - result eigenvalue will be stored here
 * @res_vect - result vector will be stoter here
 * @max - true if max eigenvalue needed,
 *			false otherwise.
 * @return 0 if succeeds, -1 otherwise
 */
int power_method(const gsl_matrix * matrix, double * res_val,
		gsl_vector * res_vect, bool max);

/*
 * @argv[1] - M number
 * @argv[2] - N number, where
 *		m - rows count
 *		n - columns count
 * @argv[3] - filename with matrix MxN
 */
int main (int argc, const char * argv[])
{
	//gsl_matrix * A = gsl_matrix_alloc(1,1);

	return 0;
}

int read_matrix_from_file(const char * file, gsl_matrix * matrix,
		const unsigned int m, const unsigned int n) {


	return 0;
}

int power_method(const gsl_matrix * matrix, double * res_val,
		gsl_vector * res_vect, bool max) {


	return 0;
}
