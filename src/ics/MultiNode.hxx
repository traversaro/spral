#pragma once

#include "SingleNode.hxx"

#include <vector>

namespace spral {
namespace ics {

template<typename T>
class MultiNode {
public:
   MultiNode() {} // FIXME: remove

   template <typename node_itr>
   MultiNode(int matrix_n, node_itr const& nbegin, node_itr const& nend, long &loffset, long &max_work_size)
   : matrix_n_(matrix_n)
   {
      contrib_size_ = 0;
      for(auto node=nbegin; node!=nend; ++node) {
         int m = node->get_nrow();
         long n = node->get_ncol();
         nodes_.push_back(new SingleNode<T>(*node));
         nodes_.back()->set_memloc(loffset, m);
         loffset += m*n;
         coffset_.push_back(contrib_size_);
         contrib_size_ += (m-n)*(m-n);
         ldcontrib_.push_back(int(m-n));
      }
      max_work_size = contrib_size_;
   }

   ~MultiNode() {
      for(auto map : maps_)
         delete map;
      for(auto node : nodes_)
         delete node;
   }

   void factor(T const* aval, T* lval, WorkspaceManager &memhandler) const {
      T *contrib = memhandler.get<T>(contrib_size_);
      int *work = memhandler.get<int>(matrix_n_);

      int idx;
      for(auto node : nodes_) {
         long coffset = coffset_[idx];
         int ldcontrib = ldcontrib_[idx];
         node->factor_local(aval, lval, &contrib[coffset], ldcontrib);
         idx++;
      }

      for(auto map : maps_)
         map->apply(lval, contrib, work);

      memhandler.release<int>(work, matrix_n_);
      memhandler.release<T>(contrib, contrib_size_);
   }

#if 0
   void factor_local(
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

   void forward_solve(int nrhs, T* x, int ldx, T const* lval,
         WorkspaceManager &memhandler) const {
      for(auto node : nodes_)
         node->forward_solve(nrhs, x, ldx, lval, memhandler);
   }

   void backward_solve(int nrhs, T* x, int ldx, T const* lval,
         WorkspaceManager &memhandler) const {
      for(auto node : nodes_)
         node->backward_solve(nrhs, x, ldx, lval, memhandler);
   }

   void print(T const* lval) const {
      for(auto node : nodes_)
         node->print(lval);
   }

   /** Return a representative index in range [0,nnodes-1] */
   int get_idx() const {
      return nodes_.front()->get_idx();
   }

   void build_contribution_map(SingleNode<T> const& ancestor) {
      maps_.push_back( new MultiMap<T>(ancestor, *this) );
   }

   void build_contribution_map(MultiNode<T> const& ancestor) {
      maps_.push_back( new MultiMap<T>(ancestor, *this) );
   }

private:
   /* Friends */
   friend class NodeToMultiMap<T>;
   friend class MultiMap<T>;

   /* Members */
   int matrix_n_;
   long contrib_size_;
   std::vector<SingleNode<T>*> nodes_;
   std::vector<long> coffset_; //< Offset into contrib for each node
   std::vector<int> ldcontrib_; //< Leading dimensions of contrib blocks
   std::vector<MultiMap<T>*> maps_;
};

} /* namespace ics */
} /* namespace spral */
