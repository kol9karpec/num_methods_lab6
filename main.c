#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#define DEF_VALUE 10
#define PRECISION 0.00001

#define str(a) #a
#define CHECK(a) \
	if(!(a)) { \
		printf("Wrong result value of an expression: (%s) in function %s()\n", \
				str(a),__func__); \
		return 1;\
	}

/*
 * Finds max eigenvalue and appropriate eigenvector by power method.
 *
 * @res_val - result eigenvalue will be stored here
 * @res_vect - result vector will be stoter here
 * @return 0 if succeeds, 1 otherwise
 */
int power_method(const gsl_matrix * matrix, double * res_val_max,
		gsl_vector * res_vect_max);

int power_method_min(const gsl_matrix * matrix, double * res_val_min,
		gsl_vector * res_vect_min);

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

	gsl_vector * max_eigenvect = gsl_vector_alloc(n);
	CHECK(max_eigenvect);
	double max_eigenval = 0;
	gsl_vector * min_eigenvect = gsl_vector_alloc(n);
	CHECK(min_eigenvect);
	double min_eigenval = 0;
	CHECK(!power_method(A,&max_eigenval,max_eigenvect,&min_eigenval,
				min_eigenvect));

	printf("Max eigenval: %lf\n",max_eigenval);
	printf("Eigenvect: \n");
	gsl_vector_fprintf(stdout,max_eigenvect,"%f");

	printf("Min eigenval: %lf\n",min_eigenval);
	printf("Eigenvect: \n");
	gsl_vector_fprintf(stdout,min_eigenvect,"%f");

	return 0;
}

int power_method(const gsl_matrix * matrix, double * res_val_max,
		gsl_vector * res_vect_max) {
	/* Counting max eigenvalue and appropriate eigenvector */
	gsl_vector * ksi = gsl_vector_alloc(matrix->size2);
	CHECK(ksi);
	double mod_ksi = 0;
	double mod_x = 0;
	double lambda_old = 0;

	CHECK(res_vect_max);
	CHECK(res_val_max);
	gsl_vector_set_all(res_vect_max,DEF_VALUE);
	*res_val_max = 0; //if not working - better to do here first iteration

	do {
		lambda_old = *res_val_max;

		/* ksi = A*x(i) */
		CHECK(!gsl_blas_dgemv( CblasNoTrans, 1.0, matrix, res_vect_max, 0.0, ksi ));

		/* |lambda|=||ksi||/||x(i)|| */
		mod_ksi = gsl_blas_dnrm2(ksi);
		mod_x = gsl_blas_dnrm2(res_vect_max);
		*res_val_max = mod_ksi / mod_x;

		/* x(i+1)=ksi/||ksi|| */
		CHECK(!gsl_vector_scale(ksi,(1 / mod_ksi)));
		CHECK(!gsl_vector_memcpy(res_vect_max, ksi));

		printf("Difference: %f\n",fabs(lambda_old - *res_val_max));
	} while(fabs(lambda_old - *res_val_max) >= PRECISION);

	gsl_vector_free(ksi);

	/* Determining the sign of founded lambda
	 * Av ?= |lambda|v
	 */

	// Av
	gsl_vector * A_mul_v = gsl_vector_alloc(matrix->size2);
	CHECK(A_mul_v);
	CHECK(!gsl_blas_dgemv( CblasNoTrans, 1.0, matrix, res_vect_max, 0.0, A_mul_v ));

	// |lambda|v
	gsl_vector * lambda_v = gsl_vector_alloc(matrix->size2);
	CHECK(lambda_v);

	CHECK(!gsl_vector_memcpy(lambda_v, res_vect_max));
	CHECK(!gsl_vector_scale(lambda_v,*res_val_max));

	/*
	printf("Comparing vectors: \n");
	gsl_vector_fprintf(stdout,A_mul_v,"%f");
	printf("\n");
	gsl_vector_fprintf(stdout,lambda_v,"%f");
	*/

	//TODO: right compare signs of vectors
	if(gsl_vector_get(A_mul_v,0)*gsl_vector_get(lambda_v,0) < 0) {
		*res_val_max *= -1;
	}

	gsl_vector_free(A_mul_v);
	gsl_vector_free(lambda_v);

	return 0;
}

int power_method_min(const gsl_matrix * matrix, double * res_val_min,
		gsl_vector * res_vect_min) {
	CHECK(res_val_min);
	CHECK(res_vect_min);




	return 0;
}
