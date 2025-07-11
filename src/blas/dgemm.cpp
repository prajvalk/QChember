#include "newscf/lapackinterface.hpp"

extern "C" {
	void dgemm_ (char* TRANSA, char* TRANSB, int* M, int* N, int* K, double* alpha,
	double* A, int* LDA,double* B, int* LDB, double* beta, double* C, int* LDC);
}

void newscf_dgemm (char* TRANSA, char* TRANSB, int* M, int* N, int* K, double* alpha,
	double* A, int* LDA,double* B, int* LDB, double* beta, double* C, int* LDC) {
	dgemm_ (TRANSA, TRANSB, M, N, K, alpha, A, LDA, B, LDB, beta, C, LDC);
}