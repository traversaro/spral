#include "blas_templates.hxx"

namespace spral {
namespace ics {

/***************************************************************************
 * BLAS
 ***************************************************************************/

/* GEMM */
extern "C" {
   void dgemm_(char const* transa, char const* transb, int const*m,
         int const* n, int const* k, double const* alpha, double const* a,
         int const* lda, double const* b, int const* ldb, double const* beta,
         double* c, int const* ldc);
}
template<> void gemm<double>(char transa, char transb, int m, int n, int k,
      double alpha, double const* a, int lda, double const* b, int ldb,
      double beta, double* c, int ldc) {
   dgemm_(&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c,
         &ldc);
}

/* SYRK */
extern "C" {
   void dsyrk_(char const* uplo, char const* trans, int const* n, int const* k,
         double const* alpha, double const* a, int const* lda,
         double const* beta, double* c, int const* ldc);
}
template<> void syrk<double>(char uplo, char trans, int n, int k, double alpha,
      double const* a, int lda, double beta, double* c, int ldc) {
   dsyrk_(&uplo, &trans, &n, &k, &alpha, a, &lda, &beta, c, &ldc);
}

/* TRSM */
extern "C" {
   void dtrsm_(char const* side, char const* uplo, char const* trans,
         char const *diag, int const* m, int const* n, double const* alpha,
         double const* a, int const* lda, double* b, int const* ldb);
}
template<> void trsm<double>(char side, char uplo, char trans, char diag,
      int m, int n, double alpha, double const* a, int lda, double* b,
      int ldb) {
   dtrsm_(&side, &uplo, &trans, &diag, &m, &n, &alpha, a, &lda, b, &ldb);
}

/***************************************************************************
 * LAPACK
 ***************************************************************************/

/* POTRF */
extern "C" {
   void dpotrf_(char const* uplo, int const* n, double* a, int const* lda,
         int *info);
}
template<> int potrf<double> (char uplo, int n, double* a, int lda) {
   int info;
   dpotrf_(&uplo, &n, a, &lda, &info);
   return info;
}

}
} /* spral */
