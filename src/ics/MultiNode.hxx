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
         nodes_.emplace_back(new SingleNode<T>(*node), contrib_size_, int(m-n));
         nodes_.back().node->set_memloc(loffset, m);
         loffset += m*n;
         contrib_size_ += (m-n)*(m-n);
      }
      max_work_size = contrib_size_;
   }

   ~MultiNode() {
      for(auto map : maps_)
         delete map;
      for(auto node : nodes_)
         delete node.node;
   }

   void factor(T const* aval, T* lval, WorkspaceManager &memhandler) const {
      T *contrib = memhandler.get<T>(contrib_size_);
      int *work = memhandler.get<int>(matrix_n_);

      for(auto nodeinfo : nodes_) {
         nodeinfo.node->factor_local(
               aval, lval, &contrib[nodeinfo.coffset], nodeinfo.ldcontrib
               );
      }

      for(auto map : maps_)
         map->apply(lval, contrib, work);

      memhandler.release<int>(work, matrix_n_);
      memhandler.release<T>(contrib, contrib_size_);
   }

   void factor_local(T const *aval, T *lval, T *contrib) {
   }

   void forward_solve(int nrhs, T* x, int ldx, T const* lval,
         WorkspaceManager &memhandler) const {
      for(auto node : nodes_)
         node.node->forward_solve(nrhs, x, ldx, lval, memhandler);
   }

   void backward_solve(int nrhs, T* x, int ldx, T const* lval,
         WorkspaceManager &memhandler) const {
      for(auto node : nodes_)
         node.node->backward_solve(nrhs, x, ldx, lval, memhandler);
   }

   void print(T const* lval) const {
      for(auto node : nodes_)
         node.node->print(lval);
   }

   /** Return a representative index in range [0,nnodes-1] */
   int get_idx() const {
      return nodes_.front().node->get_idx();
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

   struct node_info {
      SingleNode<T> *node;
      long coffset;
      int ldcontrib;
      node_info(SingleNode<T>* node, long coffset, int ldcontrib)
      : node(node), coffset(coffset), ldcontrib(ldcontrib)
      {}
   };

   /* Members */
   int matrix_n_;
   long contrib_size_;
   std::vector<node_info> nodes_;
   std::vector<MultiMap<T>*> maps_;
};

} /* namespace ics */
} /* namespace spral */
