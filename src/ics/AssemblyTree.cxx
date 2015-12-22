#include "AssemblyTree.hxx"

#include <stdexcept>

#include "spral_core_analyse.h"

namespace { /*anonymous namespace for internal fns */

/* Convert a lower triangle only matrix to full format.
 * On entry ptr_out[] must have dimension n+3, and row_out[] must be of
 * sufficient size to hold the matrix (easiest: make it twice size of row[]) */
void lwr_to_full(int n, int const ptr[], int const row[],
      int ptr_out[], int row_out[]) {
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
         row_out[ptr_out[i+1]] = k; ptr_out[i+1]++;
         if (i!=k) {
            row_out[ptr_out[k+1]] = i; ptr_out[k+1]++;
         }
      }
   }
}

} /* anonymous namespace */

namespace spral {
namespace ics {

AssemblyTree::AssemblyTree (int n, int const ptr[], int const row[],
      int perm[], int nemin)
      : nnodes_(0), sptr_(nullptr), sparent_(nullptr), rptr_(nullptr),
        rlist_(nullptr), nfact_(0), nflop_(0) {
   /* Construct full matrix from lwr triangle */
   int *ptr_full = new int[n+3];
   int *row_full = new int[2*ptr[n]];
   lwr_to_full(n, ptr, row, ptr_full, row_full);

   /* Call analysis routine to construct tree
    * NB following call allocates sptr, sparent, rptr, rlist using malloc() */
   int flag = spral_core_analyse_basic_analyse(n, ptr_full, row_full, perm,
      &nnodes_, &sptr_, &sparent_, &rptr_, &rlist_, nemin, &nfact_, &nflop_, 0);
   switch(flag) {
   case 0: break; // Normal return
   case 1: // Analyse detected matrix is structurally singular [only some cases]
      printf("Warning: matrix is structurally singular\n");
      break;
   default:
      throw std::runtime_error("spral_core_analyse_basic_analyse() failed");
   }

   /* We're now done with our full matrix copy - release memory */
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
