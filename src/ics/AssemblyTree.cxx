#include "AssemblyTree.hxx"

#include <stdexcept>

#include "spral_core_analyse.h"

namespace { /*anonymous namespace for internal fns */

/* Convert a lower triangle only matrix to full format. Produce map to trace
 * back where indices came from.
 *
 * On entry ptr_out[] must have dimension n+3, and row_out[], map[] must be of
 * sufficient size to hold the matrix (easiest: make it twice size of row[]) */
void lwr_to_full(int n, int const ptr[], int const row[],
      int ptr_out[], int row_out[], int map[]) {
   /* Count entries in each column at offset of +2 */
   for(int i=0; i<n+3; i++) ptr_out[i] = 0;
   for(int i=0; i<n; i++) {
      for(int j=ptr[i]; j<ptr[i+1]; j++) {
         int k = row[j];
         ptr_out[i+2]++;
         if (i!=k) ptr_out[k+2]++;
      }
   }
   /* Determine column starts at offset of +1 */
   ptr_out[0] = 0; ptr_out[1] = 0;
   for(int i=1; i<n+1; i++) ptr_out[i+1] += ptr_out[i];
   /* Drop entries into place */
   for(int i=0; i<n; i++) {
      for(int j=ptr[i]; j<ptr[i+1]; j++) {
         int k = row[j];
         row_out[ptr_out[i+1]] = k;
         map[ptr_out[i+1]] = j;
         ptr_out[i+1]++;
         if (i!=k) {
            row_out[ptr_out[k+1]] = i;
            map[ptr_out[k+1]] = j;
            ptr_out[k+1]++;
         }
      }
   }
}

} /* anonymous namespace */

namespace spral {
namespace ics {

void AssemblyTree::Node::construct_row_map (int *map) const {
   for(long ridx = tree_.rptr_[idx]; ridx < tree_.rptr_[idx+1]; ++ridx) {
      int row = tree_.rlist_[ridx];
      map[row] = ridx - tree_.rptr_[idx];
   }
}

void AssemblyTree::construct_tree (int const ptr[], int const row[],
      int perm[], int nemin) {

   /* Construct full matrix from lwr triangle */
   int *ptr_full = new int[n_+3];
   int *row_full = new int[2*ptr[n_]];
   int *a_half_map = new int[2*ptr[n_]];
   lwr_to_full(n_, ptr, row, ptr_full, row_full, a_half_map);

   /* Call analysis routine to construct tree
    * NB following call allocates sptr, sparent, rptr, rlist using malloc() */
   int flag = spral_core_analyse_basic_analyse(n_, ptr_full, row_full, perm,
      &nnodes_, &sptr_, &sparent_, &rptr_, &rlist_, nemin, &nfact_, &nflop_, 0);
   switch(flag) {
   case 0: break; // Normal return
   case 1: // Analyse detected matrix is structurally singular [only some cases]
      printf("Warning: matrix is structurally singular\n");
      break;
   default:
      throw std::runtime_error("spral_core_analyse_basic_analyse() failed");
   }

   /* Construct inverse permutation */
   int *inverse_perm = new int[n_];
   for(int i=0; i<n_; ++i)
      inverse_perm[ perm[i] ] = i;

   /* Construct lists (idx in aval -> idx in destination node) */
   int *map = new int[n_];
   a_to_l_ptr_.clear(); a_to_l_ptr_.reserve(nnodes_+1);
   a_to_l_map_.clear(); a_to_l_map_.reserve(ptr[n_]);
   for(int node=0; node<nnodes_; ++node) {
      a_to_l_ptr_.push_back( a_to_l_map_.size() );
      Node(*this, node).construct_row_map(map);
      for(int col=sptr_[node]; col<sptr_[node+1]; ++col) {
         int src_col = inverse_perm[col];
         for(int i=ptr[src_col]; i<ptr[src_col+1]; ++i) {
            int r = perm[ row[i] ];
            if(r < col) continue; // In upper triangle of L: ignore
            a_to_l_map_.push_back(ALMap(
                  a_half_map[i], // Source in A
                  map[r],        // Destination row of node in L
                  col            // Destination col of node in L
                  ));
         }
      }
   }
   a_to_l_ptr_.push_back( a_to_l_map_.size() );
   delete[] map;

   /* Done with full matrix and inverse perm: release memory */
   delete[] inverse_perm;
   delete[] a_half_map;
   delete[] ptr_full;
   delete[] row_full;
   
   /* Build leaf first ordering */
   build_leaf_first_order();
}

void AssemblyTree::build_leaf_first_order() {
   /* Find maximum child level from bottom - leaves count as level 0.
    * Note: As tree is stored in post order, guarunteed all children occur
    *       prior to their parents. */
   std::vector<int> level(nnodes_, 0);
   for(int i=0; i<nnodes_; i++) {
      int parent = sparent_[i];
      if(parent >= nnodes_) continue; // No parent: is a root
      if(level[parent] < level[i]+1) level[parent] = level[i]+1;
   }

   /* Construct ordering by looping over each potential level until all nodes
    * are assigned (probably inefficient, but is simple) */
   leaf_first_order_.reserve(nnodes_); // Ensure we don't need to resize
   int current_level = 0;
   while(leaf_first_order_.size() != (unsigned int) nnodes_) {
      for(int i=0; i<nnodes_; i++) {
         if(level[i]==current_level) leaf_first_order_.push_back(i);
      }
      current_level++;
   }

   /* Construct inverse temporarily for lookup purposes */
   std::vector<int> inverse_lookup(nnodes_);
   for(int i=0; i<nnodes_; i++) inverse_lookup[ leaf_first_order_[i] ] = i;

   /* Build prerequisite list: posn of latest child node in leaf_first_order_ */
   leaf_prereq_.clear();
   leaf_prereq_.resize(nnodes_, -1); // Init to no dependency for all nodes
   for(int i=0; i<nnodes_; i++) {
      if(sparent_[i] >= nnodes_) continue; // This is a root, no parent
      int parent = inverse_lookup[ sparent_[i] ];
      leaf_prereq_[ parent ] = i;
   }
}

} /* namespace ics */
} /* namespace spral */
