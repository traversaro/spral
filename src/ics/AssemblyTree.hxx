#pragma once

#include <vector>

#include "boost/iterator/iterator_facade.hpp"

namespace spral {
namespace ics {

class AssemblyTree {
public:
   /** Represents a single node of the assembly tree */
   class Node {
   public:
      /** Constructor */
      Node(AssemblyTree const& tree, int idx)
      : idx(idx), tree_(tree)
      {}

      /***********************
       * Node specific getters
       ***********************/
      /** Returns number of rows in node */
      int get_nrow(void) const {
         return tree_.rptr_[idx+1] - tree_.rptr_[idx];
      }
      /** Returns number of columns in node */
      int get_ncol(void) const {
         return tree_.sptr_[idx+1] - tree_.sptr_[idx];
      }
      /** Returns prerequisite node in leaf-first order */
      int get_leaf_prereq(void) const {
         return tree_.leaf_prereq_[idx];
      }
      /** Return parent node */
      Node getParentNode(void) const {
         return Node(tree_, tree_.sparent_[idx]);
      }

      /***********************
       * Tree-wide getters
       ***********************/
      /** Returns largest row index that can occur in tree */
      int get_max_tree_row_index(void) const {
         return tree_.n_;
      }

      /***********************
       * General routines
       ***********************/
      /** Returns true if node contains given column, or false otherwise. */
      bool contains_column(int col) const {
         return (col >= tree_.sptr_[idx] && col < tree_.sptr_[idx+1]);
      }
      /** Constructs a lookup array from relevant entries of rlist[].
       *  Does not zero any non-present entries! */
      void construct_row_map(int *map) const;

      /***********************
       * Public members
       ***********************/
      int const idx; //< Index within tree (elimination post-order)

   private:
      /* Private members */
      AssemblyTree const& tree_;
   };

   /** Iterator over tree */
   class node_iterator
      : public boost::iterator_facade<
        node_iterator,
        Node,
        boost::forward_traversal_tag,
        Node
        >
   {
   public:
      explicit node_iterator(
            AssemblyTree const& tree,
            std::vector<int>::iterator const it
            )
      : tree_(tree), it_(it)
      {}
   private:
      friend class boost::iterator_core_access;

      void increment() { ++it_; }
      bool equal(node_iterator const& other) const {
         return (this->it_ == other.it_);
      }
      Node dereference() const {
         return Node(tree_, *it_);
      }

      AssemblyTree const& tree_;
      std::vector<int>::iterator it_;
   };

   /** Construct assembly tree. IMPORTANT: Changes perm */
   AssemblyTree(int n, int const ptr[], int const row[], int perm[],
         int nemin);
   ~AssemblyTree () {
      if(sptr_) free(sptr_); // Allocated by malloc in C interface fn
      if(sparent_) free(sparent_); // Allocated by malloc in C interface fn
      if(rptr_) free(rptr_); // Allocated by malloc in C interface fn
      if(rlist_) free(rlist_); // Allocated by malloc in C interface fn
   }

   long getNfact() { return nfact_; }
   long getNflop() { return nflop_; }

   /* Iterators for leaf first ordering */
   node_iterator leaf_first_begin() {
      return node_iterator(*this, leaf_first_order_.begin());
   }
   node_iterator leaf_first_end() {
      return node_iterator(*this, leaf_first_order_.end());
   }

private:
   /** Constructs the leaf first ordering */
   void build_leaf_first_order();

   /* Core tree data */
   int const n_;
   int nnodes_;
   int *sptr_;
   int *sparent_;
   long *rptr_;
   int *rlist_;

   /* Iteration data */
   std::vector<int> leaf_first_order_; // Ordering that goes from leaves up
   std::vector<int> leaf_prereq_; // Latest child idx in leaf_first_order

   /* Stats */
   long nfact_;
   long nflop_;
};

} /* namespace ics */
} /* namespace spral */
