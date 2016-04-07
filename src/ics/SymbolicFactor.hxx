#pragma once

#include "AssemblyTree.hxx"
#include "Chunk.hxx"
#include "SingleNode.hxx"
#include "MultiNode.hxx"

#include <cstdlib>

namespace spral {
namespace ics {

class SymbolicFactor {
   typedef double T;
public:
   /** Performs a symbolic factorization as part of the construction */
   SymbolicFactor (int n, int ptr[], int row[], int nemin);
   ~SymbolicFactor () {
      if(perm_) delete[] perm_;
   }

   /** Returns length of factor array to be used */
   long get_factor_mem_size(void) const {
      return factor_mem_size_;
   }

   /** Returns maximum size of workspace to be allocated */
   long get_max_workspace_size(void) const {
      return max_workspace_size_;
   }

   /** Returns iterator to beginning of chunk list */
   std::vector< Chunk<T> >::const_iterator chunk_begin(void) const {
      return chunks_.cbegin();
   }
   /** Returns iterator to end of chunk list */
   std::vector< Chunk<T> >::const_iterator chunk_end(void) const {
      return chunks_.cend();
   }

   /** Returns iterator to reverse beginning of chunk list */
   std::vector< Chunk<T> >::const_reverse_iterator chunk_rbegin(void) const {
      return chunks_.crbegin();
   }
   /** Returns iterator to reverse end of chunk list */
   std::vector< Chunk<T> >::const_reverse_iterator chunk_rend(void) const {
      return chunks_.crend();
   }

   /* Information */
   const int nemin;
   long get_nfact() const { return tree_.get_nfact(); }
   long get_nflop() const { return tree_.get_nflop(); }
   int get_n() const { return n_; }
   int get_nnodes() const { return tree_.get_nnodes(); }
   int const* get_perm() const { return perm_; }
   int get_nchunks() const { return chunks_.size(); }

private:
   /* Core data */
   int n_;
   int *perm_;
   long factor_mem_size_;
   long max_workspace_size_;
   AssemblyTree tree_;
   std::vector< Chunk<T> > chunks_;
};

} /* namespace ics */
} /* namespace spral */
