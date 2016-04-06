#pragma once

#include <cstdlib>

#include "AssemblyTree.hxx"
#include "SingleNode.hxx"

namespace spral {
namespace ics {

class SymbolicFactor {
   typedef double T;
public:

   class Chunk {
   public:
      Chunk(SymbolicFactor const& sf, SingleNode<T> *sn)
      : sf_(sf), sn_(sn)
      {}

      friend
      void add_relation(Chunk &child, Chunk &parent) {
         child.parents_.push_back(&parent);
         parent.children_.push_back(&child);
      }

      /** Return true if this chunk has a parent */
      bool has_parent() const {
         return (parents_.size() > 0);
      }

      /** Return reference to parent Chunk */
      Chunk const& get_parent() const {
         return *parents_.front();
      }

      /** Get chunk index of this node */
      int get_idx() const {
         return sn_->get_idx();
      }

      /** Factorize all nodes in this chunk */
      void factor(T const* aval, T* lval, WorkspaceManager &memhandler) const {
         sn_->factor(aval, lval, memhandler);
      }

      void forward_solve(int nrhs, T* x, int ldx, T const* lval,
            WorkspaceManager &memhandler) const {
         sn_->forward_solve(nrhs, x, ldx, lval, memhandler);
      }

      void backward_solve(int nrhs, T* x, int ldx, T const* lval,
            WorkspaceManager &memhandler) const {
         sn_->backward_solve(nrhs, x, ldx, lval, memhandler);
      }

      void print(T const* lval) const {
         printf("CHUNK %d\n", get_idx());
         sn_->print(lval);
      }

      /** Build contribution map */
      void build_contribution_map() {
         if(!has_parent()) return; // no parent, no map

         Chunk const* parent = &get_parent();
         auto anc_end = sf_.get_ancestor_iterator_root();
         for(auto anc_itr = sf_.get_ancestor_iterator(*sn_); anc_itr != anc_end; ++anc_itr) {
            sn_->build_contribution_map(*anc_itr);
            if(parent->has_parent()) {
               parent = &parent->get_parent();
            }
         }
      }

      SingleNode<T> const* node_begin() {
         return sn_;
      }
      SingleNode<T> const* node_end() {
         return sn_+1;
      }

   private:
      SymbolicFactor const& sf_;
      std::vector<Chunk *> parents_;
      std::vector<Chunk *> children_;
      SingleNode<T> *sn_;
   };

   class ancestor_iterator
   : public boost::iterator_facade<
     ancestor_iterator,
     SingleNode<T> const,
     boost::forward_traversal_tag
     >
   {
   public:
      explicit ancestor_iterator(
         SymbolicFactor const& sfact,
         SingleNode<T> const* node
         )
      : sfact_(sfact), node_(node)
      {}
   private:
      friend class boost::iterator_core_access;

      void increment() {
         if(!node_) return; // Can't increment a root, so don't try
         int parent = node_->get_parent_idx();
         if(parent >= sfact_.get_nnodes()) {
            node_ = nullptr;
            return;
         }
         node_ = &sfact_.nodes_[parent];
      }
      bool equal(ancestor_iterator const& other) const {
         return (node_ == other.node_);
      }
      SingleNode<T> const& dereference() const {
         return *node_;
      }

      SymbolicFactor const& sfact_;
      SingleNode<T> const* node_;
   };

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
   std::vector< Chunk >::const_iterator chunk_begin(void) const {
      return chunks_.cbegin();
   }
   /** Returns iterator to end of chunk list */
   std::vector< Chunk >::const_iterator chunk_end(void) const {
      return chunks_.cend();
   }

   /** Returns iterator to reverse beginning of chunk list */
   std::vector< Chunk >::const_reverse_iterator chunk_rbegin(void) const {
      return chunks_.crbegin();
   }
   /** Returns iterator to reverse end of chunk list */
   std::vector< Chunk >::const_reverse_iterator chunk_rend(void) const {
      return chunks_.crend();
   }

   /** Returns iterator to node's ancestors (starts at parent) */
   ancestor_iterator get_ancestor_iterator(SingleNode<T> const& node) const {
      return std::next(ancestor_iterator(*this, &node), 1);
   }
   ancestor_iterator get_ancestor_iterator_root() const {
      return ancestor_iterator(*this, nullptr);
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
   std::vector< SingleNode<T> > nodes_;
   std::vector< Chunk > chunks_;

};

} /* namespace ics */
} /* namespace spral */
