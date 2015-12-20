#include "symbolic.hxx"

#include <stdexcept>

#include "spral_core_analyse.h"
#include "spral_metis_wrapper.h"

namespace spral {
namespace ics {

/* Constructs symbolic factorization from matrix data */
symbolic::symbolic (int n, int ptr[], int row[], int nemin) :
      n_(n), nnodes_(0), perm_(nullptr), sptr_(nullptr), sparent_(nullptr),
      rptr_(nullptr), rlist_(nullptr), nemin(nemin), nfact(0), nflop(0) {

   /* Perform METIS ordering */
   perm_ = new int[n];
   int *invp = new int[n];
   int flag = spral_metis_order(n, ptr, row, perm_, invp, 0);
   delete[] invp;
   if(flag)
      throw std::runtime_error("spral_metis_order() failed");

   /* Find assembly tree */
   flag = spral_core_analyse_basic_analyse(n_, ptr, row, perm_, &nnodes_,
      &sptr_, &sparent_, &rptr_, &rlist_, nemin, &nfact, &nflop, 0);
   if(flag)
      throw std::runtime_error("spral_core_analyse_basic_analyse() failed");
}

} /* namespace ics */
} /* namespace spral */
