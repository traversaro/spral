#pragma once

extern "C" {
   /* BLAS */
   void dsyrk_(char const* uplo, char const* trans, int const* n, int const* k,
         double const* alpha, double const* a, int const* lda,
         double const* beta, double* c, int const* ldc);
   void dtrsm_(char const* side, char const* uplo, char const* trans,
         char const *diag, int const* m, int const* n, double const* alpha,
         double const* a, int const* lda, double* b, int const* ldb);

   /* LAPACK */
   void dpotrf_(char const* uplo, int const* n, double* a, int const* lda,
         int *info);
}

namespace spral {
namespace ics {

/***************************************************************************
 * BLAS
 ***************************************************************************/

/* SYRK */
template <typename T> void syrk(char uplo, char trans, int n, int k, T alpha,
      T const* a, int lda, T beta, T* c, int ldc);
template <> void syrk<double>(char uplo, char trans, int n, int k, T alpha,
      T const* a, int lda, T beta, T* c, int ldc) {
   dsyrk_(&uplo, &trans, &n, &k, &alpha, a, &lda, &beta, c, &ldc);
}

/* TRSM */
template <typename T> void trsm(char side, char uplo, char trans, char diag,
      int m, int n, T alpha, T const* a, int lda, T* b, int ldb);
template <> void trsm<double>(char side, char uplo, char trans, char diag,
      int m, int n, T alpha, T const* a, int lda, double* b, int ldb) {
   dtrsm_(&side, &uplo, &trans, &diag, &m, &n, &alpha, a, &lda, b, &ldb);
}

/***************************************************************************
 * LAPACK
 ***************************************************************************/

/* POTRF */
template <typename T> int potrf(char uplo, int n, T* a, int lda);
template <> int potrf<double> (char uplo, int n, double* a, int lda) {
   int info;
   dpotrf_(&uplo, &n, a, &lda, &info);
   return info;
}

}
} /* spral */
