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
   std::vector<int> node_to_chunk(tree_.get_nnodes(), -1);
   for(auto ci=chunker.begin(); ci!=chunker.end(); ++ci) {
      if(ci->size() == 1) {
         chunks_.emplace_back(*this);
         chunks_.back().add_node(&(nodes_[ci->front().idx]));
         node_to_chunk[ci->front().idx] = chunks_.size()-1;
      } else {
         for(auto n = ci->begin(); n!=ci->end(); ++n) {
            chunks_.emplace_back(*this);
            chunks_.back().add_node(&(nodes_[n->idx]));
            node_to_chunk[n->idx] = chunks_.size()-1;
         }
      }
   }

   /* Build parent/child relations between chunks */
   std::vector<int> seen(chunks_.size(), -1);
   int idx=0;
   for(auto chunk=chunker.begin(); chunk!=chunker.end(); ++chunk) {
      for(auto node = chunk->begin(); node!=chunk->end(); ++node, ++idx) {
         if(!node->has_parent()) continue; // is a root
         int parent = node_to_chunk[node->get_parent_node().idx];
         if(seen[parent] >= idx) continue; // Already handled
         add_relation(chunks_[idx], chunks_[parent]);
      }
   }

   for(auto node=tree_.begin(); node!=tree_.end(); ++node) {
      if(node_to_chunk[node->idx] != chunker[node->idx])
         printf("hmm %d %d\n", node_to_chunk[node->idx], chunker[node->idx]);
   }

   for(auto chunk=chunks_.begin(); chunk!=chunks_.end(); ++chunk)
      chunk->build_contribution_map();
}

} /* namespace ics */
} /* namespace spral */
