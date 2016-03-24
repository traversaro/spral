#pragma once

#include "SimdVec.hxx"
#include "SingleNode.hxx"

namespace spral {
namespace ics {

/** Class that represents a chunk. That is a set of independent nodes that are
 *  handled as separate entries of a vector. */
template <
   typename T, //< Type of value
   int NVEC, //< Number of vectors worth of source columns
   int MVEC  //< Number of vectors worth of source rows to handle at a time
   >
class Chunk {
   static int const vector_length = SimdVec<T>::vector_length;
public:
   Chunk(int m, int n, int loffset, int *aidx)
   : nnodes_(0), m_(m), n_(n), loffset_(loffset), ldl_(m*vector_length),
     aidx_(aidx)
   {}

   void add_node(SingleNode<double> &node) {
      nodes_.push_back(&node);
   }

   void add_col(int idx, int const* aidx) {
      idx_[nnodes_] = idx;
      int *aidxptr = &aidx_[m_*n_*nnodes_];
      for(int i=0; i<m_*n_; i++) aidxptr[i] = aidx[i];
      nnodes_++;
   }

   void factor(T const* aval, T* lval, WorkspaceManager &memhandler) const {
      for(auto node = nodes_.begin(); node != nodes_.end(); ++node)
         (*node)->factor(aval, lval, memhandler);
   }

#if 0
   void factor(
         T const aval[],   //< Entries of A
         T lval[]          //< Entries of L
         ) {

      for(int col=0; col<n_; col++) {
         /* Each iteration of this loop handles ONE column PER NODE */
         T *lptr = &lval[loffset_ + col*NVEC*MVEC*vector_length];
         int const* aidxptr = &aidx_[col*NVEC*MVEC*vector_length];

         /* Load diagonal entries, add A, sqrt, store. */
         SimdVec<T> diag[NVEC]; // Each entry from a different node
         for(int j=0; j<NVEC; ++j) {
            diag[j] = SimdVec<T>::load_aligned( lptr[j*vector_length] );
            diag[j] += SimdVec<T>::gather(aval, aidxptr[j*vector_length], 1);
            diag[j] = sqrt(diag[j]);
            diag[j].store_aligned( lptr[j*vector_length] );
         }

         /* Load off diagonal entries, add A, divide by diag, store. */
         lptr += NVEC*vector_length;
         SimdVec<T> work[NVEC*MVEC]; // Each entry from a different node
         for(int i=0; i<MVEC; ++i) {
            for(int j=0; j<NVEC; ++j) {
               work[j*MVEC+i] =
                  SimdVec<T>::load_aligned( lptr[j*vector_length] );
               work[j*MVEC+i] += SimdVec<T>::gather(aval, aidxptr, 1);
               work[j*MVEC+i] /= diag[j];
               work[j*MVEC+i].store_aligned( lptr[j*vector_length] );
            }
         }

         /* Update rest of node */
      }
   }
#endif

   int get_nnodes() const { return nnodes_; }

private:
   int nnodes_;
   int idx_[NVEC*vector_length];
   int const m_; //< Number of rows (per node)
   int const n_; //< Number of columns (per node)
   int const loffset_; // Offset into lval of this node
   int const ldl_; // Leading dimension of lval
   int *aidx_; // Indices of aval to gather
   std::vector<SingleNode<double>*> nodes_;
};

} /* namespace ics */
} /* namespace spral */

