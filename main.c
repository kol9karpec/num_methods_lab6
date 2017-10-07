#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#define DEF_VALUE 0.5
#define PRECISION 0.00001

#define str(a) #a
#define CHECK(a) \
	if(!(a)) { \
		printf("Wrong result value of an expression: (%s) in function %s()\n", \
				str(a),__func__); \
		return 1;\
	}

/*
 * Finds max/min eigenvalues and appropriate eigenvectors by power method.
 *
 * @res_val - result eigenvalue will be stored here
 * @res_vect - result vector will be stoter here
 * @max - true if max eigenvalue needed,
 *			false if min eigenvalue needed.
 * @return 0 if succeeds, -1 otherwise
 */
int power_method(const gsl_matrix * matrix, double * res_val,
		gsl_vector * res_vect, bool max);

/*
 * @argv[1] - M number
 * @argv[2] - N number, where
 *		M - rows count
 *		N - columns count
 * @argv[3] - filename with matrix MxN
 */
int main (int argc, const char * argv[])
{
	uint32_t m = 0;
	uint32_t n = 0;

	CHECK(argc > 3);
	CHECK(sscanf(argv[1],"%u",&m) == 1);
	CHECK(sscanf(argv[2],"%u",&n) == 1);

	gsl_matrix * A = gsl_matrix_alloc(m,n);
	FILE * in_file = fopen(argv[3],"r");

	CHECK(in_file);
	CHECK(!gsl_matrix_fscanf(in_file,A));

	gsl_vector * eigenvect = gsl_vector_alloc(n);
	CHECK(eigenvect);
	double eigenval = 0;
	CHECK(!power_method(A,&eigenval,eigenvect,1));

	printf("Eigenval: %lf\n",eigenval);
	return 0;
}

int power_method(const gsl_matrix * matrix, double * res_val,
		gsl_vector * res_vect, bool max) {
	//TODO: Handle max/min eigenvalue
	gsl_vector * ksi = gsl_vector_alloc(matrix->size2);
	CHECK(ksi);
	double mod_ksi = 0;
	double mod_x = 0;
	double lambda_old = 0;

	CHECK(res_vect);
	gsl_vector_set_all(res_vect,DEF_VALUE);
	*res_val = 0; //if not working - better to do here first iteration

	do {
		lambda_old = *res_val;

		/* ksi = A*x(i) */
		CHECK(!gsl_blas_dgemv( CblasNoTrans, 1.0, matrix, res_vect, 0.0, ksi ));

		/* |lambda|=||ksi||/||x(i)|| */
		mod_ksi = gsl_blas_dnrm2(ksi);
		mod_x = gsl_blas_dnrm2(res_vect);
		*res_val = mod_ksi / mod_x;

		/* x(i+1)=ksi/||ksi|| */
		CHECK(!gsl_vector_scale(ksi,(1 / mod_ksi)));
		CHECK(!gsl_vector_memcpy(res_vect, ksi));
	} while(abs(lambda_old - *res_val) > PRECISION);

	return 0;
}
