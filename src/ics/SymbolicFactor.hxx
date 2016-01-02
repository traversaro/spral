#pragma once

#include <cstdlib>

#include "AssemblyTree.hxx"
#include "Node.hxx"

namespace spral {
namespace ics {

class SymbolicFactor {
   typedef double T;

public:
   class ancestor_iterator
   : public boost::iterator_facade<
     ancestor_iterator,
     Node<T> const,
     boost::forward_traversal_tag
     >
   {
   public:
      explicit ancestor_iterator(
         SymbolicFactor const& sfact,
         Node<T> const* node
         )
      : sfact_(sfact), node_(node)
      {}
   private:
      friend class boost::iterator_core_access;

      void increment() {
         if(!node_) return; // Can't increment a root, so don't try
         int parent = node_->get_parent_node_idx();
         if(parent >= sfact_.get_nnodes()) {
            node_ = nullptr;
            return;
         }
         node_ = &sfact_.nodes_[parent];
      }
      bool equal(ancestor_iterator const& other) const {
         return (node_ == other.node_);
      }
      Node<T> const& dereference() const {
         return *node_;
      }

      SymbolicFactor const& sfact_;
      Node<T> const* node_;
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

   /** Returns iterator to beginning of node list (in assembly tree order) */
   std::vector< Node<T> >::const_iterator node_begin(void) const {
      return nodes_.cbegin();
   }
   /** Returns iterator to end of node list (in assembly tree order) */
   std::vector< Node<T> >::const_iterator node_end(void) const {
      return nodes_.cend();
   }

   /** Returns iterator to beginning of chunk list */
   std::vector< Factorizable<T>* >::const_iterator chunk_begin(void) const {
      return chunks_.cbegin();
   }
   /** Returns iterator to end of chunk list */
   std::vector< Factorizable<T>* >::const_iterator chunk_end(void) const {
      return chunks_.cend();
   }

   /** Returns iterator to reverse beginning of node list (in assembly tree order) */
   std::vector< Node<T> >::const_reverse_iterator node_rbegin(void) const {
      return nodes_.crbegin();
   }
   /** Returns iterator to reverse end of node list (in assembly tree order) */
   std::vector< Node<T> >::const_reverse_iterator node_rend(void) const {
      return nodes_.crend();
   }

   /** Returns iterator to node's ancestors (starts at parent) */
   ancestor_iterator get_ancestor_iterator(Node<T> const& node) const {
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

private:
   /* Core data */
   int n_;
   int *perm_;
   long factor_mem_size_;
   long max_workspace_size_;
   AssemblyTree tree_;
   std::vector< Node<T> > nodes_;
   std::vector< Factorizable<T>* > chunks_;

};

} /* namespace ics */
} /* namespace spral */
