#pragma once

#include "SimdVec.cxx"

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
   const int vector_length = SimdVec<T>::vector_length;
public:
   Chunk() {
   }

   void factor(
         const T aval[],   //< Entries of A
         T lval[]          //< Entries of L
         ) {

      for(int col=0; col<n_; col++) {
         /* Each iteration of this loop handles ONE column PER NODE */
         T *lptr = &lval[loffset_ + col*NVEC*MVEC*vector_length];
         const int *aidxptr = *aidx_[col*NVEC*MVEC*vector_length];

         /* Load diagonal entries, add A, sqrt, store. */
         SimdVec<T> diag[NVEC]; // Each entry from a different node
         for(int j=0; j<NVEC; ++j) {
            diag[j] = SimdVec<T>::load_aligned( lptr[j*vector_length] );
            diag[j] += SimdVec<T>::gather(aval, aidxptr[j*vector_length], 1);
            diag[j] = sqrt(didx[j]);
            diag[j].store_aligned( lptr[j*vector_length] );
         }

         /* Load off diagonal entries, add A, divide by diag, store. */
         lptr += NVEC*vector_length;
         SimdVec<T> work[NVEC*MVEC]; // Each entry from a different node
         for(int i=0; i<MVEC; ++i) {
            for(int j=0; j<NVEC; ++j) {
               work[j*MVEC+i] =
                  SimdVec<T>::load_aligned( lptr[j*vector_length] );
               work[j*MVEC+i] += SimdVec<T>::gather(aval, aidx, 1);
               work[j*MVEC+i] /= diag[j];
               work[j*MVEC+i].store_aligned( lptr[j*vector_length] );
            }
         }

         /* Update rest of node */
      }
   }

private:
   int loffset; // Offset into lval of this node
   int ldl; // Leading dimension of lval
   int *aidx; // Indices of aval to gather
};

} /* namespace ics */
} /* namespace spral */

