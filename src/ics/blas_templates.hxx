#pragma once

namespace spral {
namespace ics {

/***************************************************************************
 * BLAS
 ***************************************************************************/

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
