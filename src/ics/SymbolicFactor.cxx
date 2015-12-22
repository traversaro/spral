#include "SymbolicFactor.hxx"

#include <stdexcept>

#include "AssemblyTree.hxx"
#include "spral_metis_wrapper.h"

namespace spral {
namespace ics {

/* Constructs symbolic factorization from matrix data */
SymbolicFactor::SymbolicFactor (int n, int ptr[], int row[], int nemin) :
      nemin(nemin), nfact(0), nflop(0), n_(n), nnodes_(0), perm_(nullptr),
      factor_mem_size_(0) {

   /* Perform METIS ordering */
   perm_ = new int[n];
   int *invp = new int[n];
   int flag = spral_metis_order(n, ptr, row, perm_, invp, 0);
   delete[] invp;
   if(flag)
      throw std::runtime_error("spral_metis_order() failed");

   /* Construct AssemblyTree */
   AssemblyTree tree(n, ptr, row, perm_, nemin);

   /* Construct chunk buckets */
   for(auto nitr=tree.leaf_first_begin(); nitr!=tree.leaf_first_end(); ++nitr) {
      AssemblyTree::Node node = *nitr;
      printf("Node %d is %d x %d\n", node.idx, node.get_nrow(), node.get_ncol());
   }

}

} /* namespace ics */
} /* namespace spral */
