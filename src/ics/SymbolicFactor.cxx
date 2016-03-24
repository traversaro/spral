#include "SymbolicFactor.hxx"

#include <stdexcept>

#include "AssemblyTree.hxx"
#include "Chunk.hxx"
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

   /* Construct list of chunks */
   Chunker chunker(tree_);
   /*for(auto node=tree_.begin(); node!=tree_.end(); ++node) {
      printf("Node %d (%d x %d) in chunk %d\n", node->idx, node->get_nrow(),
            node->get_ncol(), chunker[*node]);
   }*/

   /* Construct list of nodes */
   int max_contrib_size = 0;
   for(auto node=tree_.begin(); node!=tree_.end(); ++node) {
      int m = node->get_nrow();
      int n = node->get_ncol();
      nodes_.push_back(SingleNode<double>(*node));
      max_contrib_size = std::max(max_contrib_size, m-n);
   }
   max_workspace_size_ = max_contrib_size*max_contrib_size*sizeof(double) +
      n_*sizeof(int);

   /* Assign memory locations so chunks are contigous */
   factor_mem_size_ = 0;
   for(auto ci=chunker.begin(); ci!=chunker.end(); ++ci) {
      for(auto node = ci->begin(); node!=ci->end(); ++node) {
         int m = node->get_nrow();
         int n = node->get_ncol();
         nodes_[node->idx].set_memloc(factor_mem_size_, m);
         factor_mem_size_ += m*((long) n);
         if(node->has_parent())
            nodes_[node->get_parent_node().idx].add_child();
      }
   }

   /* Now we know where nodes are in memory, build map */
   int *map = new int[n_];
   for(auto node=nodes_.begin(); node!=nodes_.end(); ++node) {
      node->build_contribution_map(
            get_ancestor_iterator(*node), get_ancestor_iterator_root(), map
            );
   }
   for(auto ci=chunker.begin(); ci!=chunker.end(); ++ci) {
      if(ci->size() == 1) {
         chunks_.push_back(&(nodes_[ci->front().idx]));
      } else {
         for(auto n = ci->begin(); n!=ci->end(); ++n) {
            chunks_.push_back(&(nodes_[n->idx]));
         }
      }
   }
   delete[] map;
}

} /* namespace ics */
} /* namespace spral */
