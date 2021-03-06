#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#define DEF_VALUE 0.5
#define PRECISION 0.00001
#define DEF_OUTPUT stdout
#define CUR_OUTPUT DEF_OUTPUT
#define DEF_EPS 0.000001

#define str(a) #a
#define CHECK(a) \
	if(!(a)) { \
		printf("Wrong result value of an expression: (%s) in function %s()\n", \
				str(a),__func__); \
		return 1;\
	}
#define PRINT(a) printf("%s: %0.6lf\n",str(a),a);

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
		gsl_vector * res_vect_min, const double max_eigenval);

int LR_method(FILE * stream, const gsl_matrix * matrix, gsl_vector * eigenvals,
		gsl_matrix * eigenvectors);

void gsl_matrix_print(FILE * stream, gsl_matrix * a);

void gsl_vector_print(FILE * stream, gsl_vector * v);

void iteration_printf(FILE * stream, gsl_matrix * A, gsl_matrix * L,
		gsl_matrix * R, const uint32_t iteration_number);

double frobenius_norm(const gsl_matrix * matr);

int find_eigenvects_triang_matr(const gsl_matrix * A, gsl_matrix * eigenvects_dest);

int diff_vect(const gsl_matrix * A, const gsl_vector * eigenvals,
		const gsl_matrix * eigenvects, gsl_matrix * res);

/*
 * Returns false if |aij-bij| > eps
 */
bool check_diff_mm(const gsl_matrix * A, const gsl_matrix * B,
		const double eps);

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
	CHECK(!power_method(A,&max_eigenval,max_eigenvect));

	fprintf(CUR_OUTPUT,"----------Power method---------\n");

	fprintf(CUR_OUTPUT,"Max eigenval: %lf\n",max_eigenval);
	fprintf(CUR_OUTPUT,"Eigenvect: \n");
	gsl_vector_fprintf(CUR_OUTPUT,max_eigenvect,"%f");

	CHECK(!power_method_min(A,&min_eigenval, min_eigenvect, max_eigenval));
	fprintf(CUR_OUTPUT,"Min eigenval: %lf\n",min_eigenval);
	fprintf(CUR_OUTPUT,"Eigenvect: \n");
	gsl_vector_fprintf(CUR_OUTPUT,min_eigenvect,"%f");
	fprintf(CUR_OUTPUT,"\n");

	gsl_vector * eigenvals = gsl_vector_alloc(n);
	gsl_matrix * eigenvects = gsl_matrix_alloc(n,n);

	fprintf(CUR_OUTPUT,"----------LR method---------\n");
	LR_method(CUR_OUTPUT,A,eigenvals,eigenvects);

	fprintf(CUR_OUTPUT,"\nSolution:\n\tEigenvavules\n");
	gsl_vector_print(CUR_OUTPUT,eigenvals);
	fprintf(CUR_OUTPUT,"\n\tAppropriate eigenvectors:\n");
	gsl_matrix_print(CUR_OUTPUT,eigenvects);

	gsl_matrix * diff_vects = gsl_matrix_alloc(n,n);
	diff_vect(A,eigenvals,eigenvects,diff_vects);
	fprintf(CUR_OUTPUT,"\tAppropriate |Av-lv|:\n");
	gsl_matrix_print(CUR_OUTPUT,diff_vects);

	fclose(in_file);
	gsl_matrix_free(A);
	gsl_vector_free(max_eigenvect);
	gsl_vector_free(min_eigenvect);
	gsl_vector_free(eigenvals);
	gsl_matrix_free(eigenvects);
	gsl_matrix_free(diff_vects);

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

		//printf("Difference: %f\n",fabs(lambda_old - *res_val_max));
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
		gsl_vector * res_vect_min, const double max_eigenval) {
	//int i = 0;
	CHECK(res_val_min);
	CHECK(res_vect_min);

	gsl_matrix * B = gsl_matrix_alloc(matrix->size1, matrix->size2);
	CHECK(B);
	CHECK(!gsl_matrix_memcpy(B,matrix));

	gsl_matrix * E = gsl_matrix_alloc(matrix->size1, matrix->size2);
	CHECK(E);
	gsl_matrix_set_identity(E);
	/*
	gsl_matrix_set_all(E,0);
	for(i=0;i<matrix->size1;i++) {
		gsl_matrix_set(E,i,i,1.0l);
	}
	*/

	CHECK(!gsl_matrix_scale(E,max_eigenval));
	CHECK(!gsl_matrix_sub(B,E));

	gsl_vector * B_max_eigenvect = gsl_vector_alloc(matrix->size2);
	CHECK(B_max_eigenvect);
	double B_max_eigenval = 0.0l;

	CHECK(!power_method(B,&B_max_eigenval,B_max_eigenvect));
	*res_val_min = max_eigenval + B_max_eigenval;
	gsl_vector_memcpy(res_vect_min,B_max_eigenvect);

	gsl_matrix_free(B);
	gsl_matrix_free(E);
	gsl_vector_free(B_max_eigenvect);

	return 0;
}

int LR_method(FILE * stream, const gsl_matrix * matrix, gsl_vector * eigenvals,
		gsl_matrix * eigenvectors) {
	gsl_matrix * A = gsl_matrix_alloc(matrix->size1, matrix->size2);
	CHECK(A);
	gsl_permutation * p = gsl_permutation_alloc(matrix->size1);
	CHECK(p);
	int signum = 0;
	uint32_t index = 0;
	int i = 0, j = 0;

	gsl_matrix_memcpy(A,matrix);

	gsl_matrix * L = gsl_matrix_alloc(matrix->size1, matrix->size2);
	gsl_matrix * R = gsl_matrix_alloc(matrix->size1, matrix->size2);
	CHECK(L);
	CHECK(R);

	gsl_matrix * A_old = gsl_matrix_alloc(matrix->size1, matrix->size2);
	gsl_matrix * L_mult = gsl_matrix_alloc(matrix->size1, matrix->size2);
	gsl_matrix * L_buf = gsl_matrix_alloc(matrix->size1, matrix->size2);
	gsl_matrix_set_identity(L_mult);
	gsl_matrix_set_zero(A_old);

	while (!check_diff_mm(A,A_old,DEF_EPS)) {
		gsl_matrix_memcpy(A_old,A);
		CHECK(!gsl_linalg_LU_decomp(A,p,&signum));
		gsl_matrix_set_identity(L);
		gsl_matrix_set_zero(R);

		/* Copying matrix L*/
		for(i=0;i<A->size1;i++) {
			for(j=0;j<i;j++) {
				gsl_matrix_set(L,i,j,gsl_matrix_get(A,i,j));
			}
		}

		/* Copying matrix R*/
		for(j=0;j<A->size1;j++) {
			for(i=0;i<=j;i++) {
				gsl_matrix_set(R,i,j,gsl_matrix_get(A,i,j));
			}
		}
		iteration_printf(stream,A_old,L,R,++index);

		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, R, L,0.0, A);
		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, L_mult, L,0.0, L_buf);
		gsl_matrix_memcpy(L_mult,L_buf);
	}

	/* Returning eigenvalues*/
	for(i=0;i<A->size1;i++)
		gsl_vector_set(eigenvals,i,gsl_matrix_get(A,i,i));

	/* Computing eigenvectors of diagonalized A*/
	gsl_matrix * eigenvects_A_diag = gsl_matrix_alloc(A->size1, A->size2);
	CHECK(!find_eigenvects_triang_matr(A, eigenvects_A_diag));

	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, L_mult, eigenvects_A_diag,
			0.0, eigenvectors);

	gsl_matrix_free(A);
	gsl_permutation_free(p);
	gsl_matrix_free(L);
	gsl_matrix_free(R);
	gsl_matrix_free(A_old);
	gsl_matrix_free(L_mult);
	gsl_matrix_free(L_buf);
	gsl_matrix_free(eigenvects_A_diag);

	return 0;
}

void gsl_matrix_print(FILE * stream, gsl_matrix * A) {
	int i,j;

	for(i=0;i<A->size1;i++) {
		for(j=0;j<A->size2;j++) {
			fprintf(stream,"%0.6lf\t",gsl_matrix_get(A,i,j));
		}
		printf("\n");
	}
}

void iteration_printf(FILE * stream, gsl_matrix * A, gsl_matrix * L,
		gsl_matrix * R, const uint32_t iteration_number) {
	fprintf(stream,"Iteration #%u\n",iteration_number);
	fprintf(stream,"A: \n");
	fprintf(stream,"Frobenius norm: %0.6lf\n",frobenius_norm(A));
	gsl_matrix_print(stdout,A);
	fprintf(stream,"L: \n");
	gsl_matrix_print(stdout,L);
	fprintf(stream,"R: \n");
	gsl_matrix_print(stdout,R);
	fprintf(stream,"-------------------------------------------------\n");
}

double frobenius_norm(const gsl_matrix * matr) {
	uint32_t i=0,j=0;
	double sum = 0;

	for(i=0;i<matr->size1;i++) {
		for(j=0;j<matr->size2;j++) {
			sum += gsl_matrix_get(matr,i,j);
		}
	}

	return sqrt(sum);
}

int find_eigenvects_triang_matr(const gsl_matrix * A, gsl_matrix * eigenvects_dest) {
	CHECK(A);
	CHECK(eigenvects_dest);
	gsl_matrix_set_zero(eigenvects_dest);

	int i = 0, j = 0, k = 0;
	double a_1 = 0.l,
		   a_2 = 0.l,
		   a_3 = 0.l,
		   u_1 = 0.l,
		   sum = 0.l;
	for(i=0;i<A->size1;i++) {
		for(j=i+1;j<A->size2;j++)
			gsl_matrix_set(eigenvects_dest, j, i, 0.l);

		gsl_matrix_set(eigenvects_dest, i, i, 1.l);

		for(j=i-1;j>=0;j--) {
			sum = 0.l;
			for(k=j+1;k<=i;k++) {
				a_1 = gsl_matrix_get(A,j,k);
				u_1 = gsl_matrix_get(eigenvects_dest,k,i);
				sum += a_1*u_1;
			}

			a_2 = gsl_matrix_get(A,i,i);
			a_3 = gsl_matrix_get(A,j,j);

			gsl_matrix_set(eigenvects_dest,j,i,(sum/(a_2-a_3)));
		}
	}

	return 0;
}

bool check_diff_mm(const gsl_matrix * A, const gsl_matrix * B,
		const double eps) {
	int i = 0, j = 0;

	for(i=0;i<A->size1;i++)
		for(j=0;j<A->size2;j++)
			if(fabs(gsl_matrix_get(A,i,j) - gsl_matrix_get(B,i,j)) >= eps)
				return false;

	return true;
}

void gsl_vector_print(FILE * stream, gsl_vector * v) {
	 int i = 0;
	 
	 for(;i<v->size;i++)
		 fprintf(stream,"%0.6lf\t",gsl_vector_get(v,i));
}

int diff_vect(const gsl_matrix * A, const gsl_vector * eigenvals,
		const gsl_matrix * eigenvects, gsl_matrix * res) {
	int i = 0, j = 0;
	gsl_vector * a_scaled = gsl_vector_alloc(A->size1);
	gsl_vector * v = gsl_vector_alloc(A->size1);

	for(;i<eigenvals->size;i++) {
		gsl_matrix_get_col(v,eigenvects,i);
		gsl_blas_dgemv(CblasNoTrans, 1.0, A, v, 0.0, a_scaled);

		gsl_vector_scale(v,gsl_vector_get(eigenvals,i));

		gsl_vector_sub(a_scaled,v);
		for(j=0;j<a_scaled->size;j++)
			if(gsl_vector_get(a_scaled,j) < 0)
				gsl_vector_set(a_scaled,j,(-1)*gsl_vector_get(a_scaled,j));
		gsl_matrix_set_col(res,i,a_scaled);
	}

	return 0;
}
