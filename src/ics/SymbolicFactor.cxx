#include "SymbolicFactor.hxx"

#include <stdexcept>

#include "spral_core_analyse.h"
#include "spral_metis_wrapper.h"

namespace { /*anonymous namespace for internal fns */

/* Convert a lower triangle only matrix to full format.
 * On entry ptr_out[] must have dimension n+3, and row_out[] must be of
 * sufficient size to hold the matrix (easiest: make it twice size of row[]) */
void lwr_to_full(int n, const int ptr[], const int row[],
      int ptr_out[], int row_out[]) {
   /* Count entries in each column at offset of +2 */
   for(int i=0; i<n+3; i++) ptr_out[i] = 0;
   for(int i=0; i<n; i++) {
      for(int j=ptr[i]; j<ptr[i+1]; j++) {
         int k = row[j];
         ptr_out[i+2]++;
         if (i!=k) ptr_out[k+2]++;
      }
   }
   /* Determine column starts at offset of +1 */
   ptr_out[0] = 0; ptr_out[1] = 0;
   for(int i=1; i<n+1; i++) ptr_out[i+1] += ptr_out[i];
   /* Drop entries into place */
   for(int i=0; i<n; i++) {
      for(int j=ptr[i]; j<ptr[i+1]; j++) {
         int k = row[j];
         row_out[ptr_out[i+1]] = k; ptr_out[i+1]++;
         if (i!=k) {
            row_out[ptr_out[k+1]] = i; ptr_out[k+1]++;
         }
      }
   }
}

} /* anonymous namespace */

namespace spral {
namespace ics {

/* Constructs symbolic factorization from matrix data */
SymbolicFactor::SymbolicFactor (int n, int ptr[], int row[], int nemin) :
      nemin(nemin), nfact(0), nflop(0), n_(n), nnodes_(0), perm_(nullptr),
      sptr_(nullptr), sparent_(nullptr), rptr_(nullptr), rlist_(nullptr) {

   /* Perform METIS ordering */
   perm_ = new int[n];
   int *invp = new int[n];
   int flag = spral_metis_order(n, ptr, row, perm_, invp, 0);
   delete[] invp;
   if(flag)
      throw std::runtime_error("spral_metis_order() failed");

   /* Find assembly tree - requires full matrix, not just lwr triangle */
   int *ptr_full = new int[n+3];
   int *row_full = new int[2*ptr[n]];
   lwr_to_full(n, ptr, row, ptr_full, row_full);
   /* NB following call allocates sptr, sparent, rptr, rlist using malloc() */
   flag = spral_core_analyse_basic_analyse(n_, ptr_full, row_full, perm_,
      &nnodes_, &sptr_, &sparent_, &rptr_, &rlist_, nemin, &nfact, &nflop, 0);
   if(flag)
      throw std::runtime_error("spral_core_analyse_basic_analyse() failed");
   delete[] ptr_full;
   delete[] row_full;
}

} /* namespace ics */
} /* namespace spral */
