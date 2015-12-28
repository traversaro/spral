#pragma once

namespace spral {
namespace ics {

/***************************************************************************
 * BLAS
 ***************************************************************************/

/* GEMM */
template<typename T> void gemm(char transa, char transb, int m, int n, int k,
      T alpha, T const* a, int lda, T const* b, int ldb, T beta, T* c, int ldc);

/* SYRK */
template<typename T> void syrk(char uplo, char trans, int n, int k, T alpha,
      T const* a, int lda, T beta, T* c, int ldc);

/* TRSM */
template<typename T> void trsm(char side, char uplo, char trans, char diag,
      int m, int n, T alpha, T const* a, int lda, T* b, int ldb);

/***************************************************************************
 * LAPACK
 ***************************************************************************/

/* POTRF */
template<typename T> int potrf(char uplo, int n, T* a, int lda);

}
} /* spral */
