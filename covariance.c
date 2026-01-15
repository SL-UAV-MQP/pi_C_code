/**
 * @file covariance.c
 * @brief Spatial covariance matrix computation
 */

#include "music_uca_6.h"
#include "music_config.h"
#include <string.h>
#include <complex.h>

/* BLAS/LAPACK function declarations (external linkage) */
/* ZHERK: Hermitian rank-k update: C := alpha*A*A^H + beta*C */
extern void zherk_(
    const char* uplo,    /* 'U' or 'L' */
    const char* trans,   /* 'N' or 'C' (conjugate transpose) */
    const int* n,        /* Order of matrix C */
    const int* k,        /* Number of columns of A */
    const double* alpha, /* Scalar alpha */
    const cdouble_t* A,  /* Matrix A */
    const int* lda,      /* Leading dimension of A */
    const double* beta,  /* Scalar beta */
    cdouble_t* C,        /* Output matrix C */
    const int* ldc       /* Leading dimension of C */
);

/* ============================================================================
 * Covariance Matrix Computation
 * ============================================================================ */

music_status_t compute_covariance_matrix(
    const signal_matrix_t* signals,
    bool use_forward_backward,
    covariance_matrix_t* covariance)
{
    // compute_covariance_matrix
    // Computes spatial covariance using BLAS ZHERK with optional FB averaging

    if (!signals || !covariance) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    int M = signals->M;
    int N = signals->N;

    if (covariance->M != M) {
        return MUSIC_ERROR_INVALID_SIZE;
    }

    // Compute forward covariance using BLAS ZHERK
    // R_forward = (1/N) * X * X^H
    const char uplo = 'U';    // Upper triangular
    const char trans = 'N';   // No transpose
    double alpha = 1.0 / N;   // Scaling factor
    double beta = 0.0;        // Initialize output

    zherk_(&uplo, &trans, &M, &N, &alpha, signals->data, &M, &beta, covariance->data, &M);

    // Fill lower triangular part (ZHERK only computes upper)
    for (int j = 0; j < M; j++) {
        for (int i = j + 1; i < M; i++) {
            // R[i,j] = conj(R[j,i])
            covariance->data[i + j * M] = conj(covariance->data[j + i * M]);
        }
    }

    if (use_forward_backward) {
        // Forward-backward averaging: R = 0.5 * (R_f + J*conj(R_f)*J)
        cdouble_t* R_backward = (cdouble_t*)malloc(M * M * sizeof(cdouble_t));
        if (!R_backward) {
            return MUSIC_ERROR_MEMORY;
        }

        // Apply exchange matrix: R_b[i,j] = conj(R_f[M-1-i, M-1-j])
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < M; j++) {
                int i_flip = M - 1 - i;
                int j_flip = M - 1 - j;
                R_backward[i + j * M] = conj(covariance->data[i_flip + j_flip * M]);
            }
        }

        // Average: R = 0.5 * (R_f + R_b)
        for (int i = 0; i < M * M; i++) {
            covariance->data[i] = 0.5 * (covariance->data[i] + R_backward[i]);
        }

        free(R_backward);
    }

    return MUSIC_SUCCESS;
}
