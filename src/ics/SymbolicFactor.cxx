#include "SymbolicFactor.hxx"

#include <stdexcept>

#include "AssemblyTree.hxx"
#include "Chunker.hxx"
#include "spral_metis_wrapper.h"

namespace spral {
namespace ics {

/* Constructs symbolic factorization from matrix data */
SymbolicFactor::SymbolicFactor (int n, int ptr[], int row[], int nemin)
: nemin(nemin), n_(n), perm_(nullptr), factor_mem_size_(0),
  max_workspace_size_(0), tree_(n)
{

   /* Perform METIS ordering */
   perm_ = new int[n];
   int *invp = new int[n];
   int flag = spral_metis_order(n, ptr, row, perm_, invp, 0);
   delete[] invp;
   if(flag)
      throw std::runtime_error("spral_metis_order() failed");

   /* Construct AssemblyTree */
   tree_.construct_tree(ptr, row, perm_, nemin);

   /* Find assignment of nodes to chunks */
   Chunker chunker(tree_);

   /* Construct chunks */
   chunks_.reserve(chunker.get_nchunks());
   for(auto ci=chunker.begin(); ci!=chunker.end(); ++ci) {
      chunks_.emplace_back(*this);
   }
   int idx=0;
   for(auto ci=chunker.begin(); ci!=chunker.end(); ++ci, ++idx) {
      long sz;
      chunks_[idx].setup(ci->begin(), ci->end(), factor_mem_size_, sz);
      max_workspace_size_ = std::max(max_workspace_size_, sz);
   }
   max_workspace_size_ = max_workspace_size_*sizeof(double) +
      n_*sizeof(int);

   /* Build parent/child relations between chunks */
   std::vector<int> seen(chunks_.size(), -1);
   idx=0;
   for(auto chunk=chunker.begin(); chunk!=chunker.end(); ++chunk, ++idx) {
      for(auto node = chunk->begin(); node!=chunk->end(); ++node) {
         if(!node->has_parent()) continue; // is a root
         int parent = chunker[node->get_parent_node().idx];
         if(seen[parent] >= idx) continue; // Already handled
         seen[parent] = idx;
         add_relation(chunks_[idx], chunks_[parent]);
      }
   }

   /* Now we know where chunks are in memory and their relation, build map */
   for(auto chunk=chunks_.begin(); chunk!=chunks_.end(); ++chunk)
      chunk->build_contribution_map(tree_.get_nnodes());
}

} /* namespace ics */
} /* namespace spral */
