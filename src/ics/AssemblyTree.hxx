#pragma once

#include <limits>
#include <vector>

#include "boost/iterator/iterator_facade.hpp"

namespace spral {
namespace ics {

class AssemblyTree {
public:
   /** Maps a single entry from A to within a node of L */
   struct ALMap {
      int src;
      int dest_row;
      int dest_col;

      ALMap(int src, int dest_row, int dest_col)
      : src(src), dest_row(dest_row), dest_col(dest_col)
      {}
   };

   /** Represents a single node of the assembly tree */
   class Node {
   public:
      /** Constructor */
      Node(AssemblyTree const& tree, int idx)
      : idx(idx), tree_(tree)
      {
         if(idx > tree_.n_) idx = std::numeric_limits<int>::max();
      }

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
      Node get_parent_node(void) const {
         return Node(tree_, tree_.sparent_[idx]);
      }
      /** Returns true if node has a parent */
      bool has_parent(void) const {
         return (tree_.sparent_[idx] < tree_.nnodes_);
      }
      /** Return first child node: descendants are exactly those nodes in range
       *  [first_child:idx-1] */
      Node get_first_child(void) const {
         return Node(tree_, tree_.first_child_[idx]);
      }
      /** Return iterator to start of a_to_l list for node */
      std::vector<struct ALMap>::const_iterator a_to_l_begin(void) const {
         return std::next(
               tree_.a_to_l_map_.cbegin(),
               tree_.a_to_l_ptr_[idx]
               );
      }
      /** Return iterator to end of a_to_l list for node */
      std::vector<struct ALMap>::const_iterator a_to_l_end(void) const {
         return std::next(
               tree_.a_to_l_map_.cbegin(),
               tree_.a_to_l_ptr_[idx+1]
               );
      }
      /** Return iterator to start of row list for node */
      int const* row_begin(void) const {
         return &tree_.rlist_[tree_.rptr_[idx]];
      }
      /** Return iterator to end of row list for node */
      int const* row_end(void) const {
         return &tree_.rlist_[tree_.rptr_[idx+1]];
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
      /** Return true if other is a descendant of this node, false o'wise. If
       *  other == *this, returns true. */
      bool has_descendant(Node const& other) const {
         return (other.idx >= tree_.first_child_[idx] &&
                 other.idx <= idx);
      }

      /***********************
       * Public members
       ***********************/
      int const idx; //< Index within tree (elimination post-order)

   private:
      /* Private members */
      AssemblyTree const& tree_;
   };

   AssemblyTree(int n)
   : n_(n), nnodes_(0), sptr_(nullptr), sparent_(nullptr), rptr_(nullptr),
     rlist_(nullptr), a_to_l_map_(), nfact_(0), nflop_(0)
   {}
   /** Construct assembly tree. IMPORTANT: Changes perm */
   AssemblyTree(int n, int const ptr[], int const row[], int perm[],
         int nemin)
   : n_(n), nnodes_(0), sptr_(nullptr), sparent_(nullptr), rptr_(nullptr),
     rlist_(nullptr), a_to_l_map_(), nfact_(0), nflop_(0)
   {
      construct_tree(ptr, row, perm, nemin);
   }
   ~AssemblyTree () {
      if(sptr_) free(sptr_); // Allocated by malloc in C interface fn
      if(sparent_) free(sparent_); // Allocated by malloc in C interface fn
      if(rptr_) free(rptr_); // Allocated by malloc in C interface fn
      if(rlist_) free(rlist_); // Allocated by malloc in C interface fn
   }

   /** Initializer. IMPORTANT: Changes perm */
   void construct_tree(int const ptr[], int const row[], int perm[], int nemin);

   int get_nnodes() const { return nnodes_; }
   long get_nfact() const { return nfact_; }
   long get_nflop() const { return nflop_; }

   /* Iterators for standard ordering (bottom to top) */
   std::vector<Node>::const_iterator begin() const {
      return nodes_.cbegin();
   }
   std::vector<Node>::const_iterator end() const {
      return nodes_.cend();
   }

   /* Iterators for standard reverse ordering (top to bottom) */
   std::vector<Node>::const_reverse_iterator rbegin() const {
      return nodes_.crbegin();
   }
   std::vector<Node>::const_reverse_iterator rend() const {
      return nodes_.crend();
   }

   /* Iterators for leaf first ordering */
   std::vector<Node>::const_iterator leaf_first_begin() const {
      return leaf_first_order_.cbegin();
   }
   std::vector<Node>::const_iterator leaf_first_end() const {
      return leaf_first_order_.cend();
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

   /* Data wrt matrix A */
   std::vector<int> a_to_l_ptr_;
   std::vector<struct ALMap> a_to_l_map_;

   /* Iteration order data */
   std::vector<Node> nodes_; // In standard ordering
   std::vector<Node> leaf_first_order_; // Ordering that goes from leaves up
   std::vector<int> leaf_prereq_; // Latest child idx in leaf_first_order
   std::vector<int> first_child_; // Index of lowest numbered child

   /* Stats */
   long nfact_;
   long nflop_;
};

} /* namespace ics */
} /* namespace spral */
