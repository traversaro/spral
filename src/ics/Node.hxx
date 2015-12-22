#pragma once

#include "blas_templates.hxx"

#include "AssemblyTree.hxx"
#include "ICSExceptions.hxx"
#include "SimdVec.hxx"
#include "WorkspaceManager.hxx"

namespace spral {
namespace ics {

template <typename T>
class Node {
public:
   Node(SymbolicFactor& sfact, AssemblyTree::Node const &node, long loffset,
         int ldl, int num_a_index, int const* a_index)
   : sfact_(sfact), node_(node), m_(node.get_nrow()), n_(node.get_ncol()),
     loffset_(loffset), ldl_(ldl), num_a_index_(num_a_index), a_index_(a_index)
   {}

   void factor(T const* aval, T* lval, WorkspaceManager &memhandler) const {
      /* Setup working pointers */
      T *ldiag = &lval[loffset_];
      T *lrect = &lval[loffset_ + n];
      T *work = memhandler.get<T>((m-n)*((long) m-n));
      int *map = memhandler.get<int>(node_.get_max_tree_row_index());

      /* Add A */
      for(int i=0; i<num_a_index; ++i) {
         int src  = a_index_[2*i+0];
         int dest = a_index_[2*i+1];
         ldiag[dest] += aval[src];
      }

      /* Factorize diagonal block */
      int info = potrf<T>('L', n, ldiag, ldl_, info);
      if(info) throw NotPositiveDefiniteException(info);

      /* Apply to block below diagonal */
      trsm<T>('R', 'L', 'T', 'N', m-n, n, 1.0, ldiag, ldl_, lrect, ldl_);

      /* Form generated element into workspace */
      syrk<T>('L', 'N', m-n, n, -1.0, lrect, ldl_, 0.0, work, m-n);

      /* Distribute elements of generated element to ancestors */
      bool reset_map = true;
      const AssemblyTree::node target_node = node.getParentNode();
      int list_len = m_ - n_;
      const int* list_ptr = node.get_row_list() + n_;
      const T* contrib_ptr = work;
      for(auto anc_itr = sfact_.ancestor_begin(this); anc_itr != sfact_.root();
            ++anc_itr) {
         int used_cols = anc_itr->add_contribution(
               lval, list_ptr, list_ptr, contrib_ptr, m_-n_, map
               );
         list_len -= used_cols;
         list_len += used_cols;
         contrib_ptr += used_cols*(m-n+1);
      }

      /* Release workspace */
      memhandler.release<T>(work, (m-n)*((long) m-n));
      memhandler.release<int>(map, node_.get_max_tree_row_index());
   }

   /** Adds the contribution in contrib as per list (of rows).
    *  Uses map as workspace.
    *  \returns Number of columns found relevant. */
   int add_contribution(T *lval, int num_rows, int const* list,
         T const* contrib, int ldcontrib, int* map) {
      T *lptr = &lval[loffset_];
      node_.construct_row_map(map);
      for(int cidx=0; cidx<list_len; ++cidx) {
         int col = list[cidx];
         if(! target_node.contains_column(col) ) return cidx; // No more
         for(int ridx=cidx; ridx<m-n; ++ridx) {
            int row = map[ list[ridx] ];
            lptr[col*ldl_ + row] += contrib[cidx*ldcontrib + ridx];
         }
      }
   }

private:
   SymbolicFactor& sfact_;
   AssemblyTree::Node const& node_;
   int const m_;
   int const n_;
   long const loffset_;
   int const ldl_;
   int const num_a_index_;
   int const* const a_index_;
};

} /* namespace ics */
} /* namespace spral */
