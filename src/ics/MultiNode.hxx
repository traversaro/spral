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
   MultiNode(node_itr const& nbegin, node_itr const& nend, long &loffset, long &max_work_size) {
      max_work_size = 0;
      for(auto node=nbegin; node!=nend; ++node) {
         int m = node->get_nrow();
         long n = node->get_ncol();
         nodes_.push_back(new SingleNode<T>(*node));
         nodes_.back()->set_memloc(loffset, m);
         loffset += m*n;
         max_work_size = std::max(max_work_size, (m-n)*(m-n));
      }
   }

   ~MultiNode() {
      for(auto node : nodes_)
         delete node;
   }

   void factor(T const* aval, T* lval, WorkspaceManager &memhandler) const {
      for(auto node : nodes_)
         node->factor(aval, lval, memhandler);
   }

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
      for(auto node : nodes_)
         node->build_contribution_map(ancestor);
   }

   void build_contribution_map(MultiNode<T> const& ancestor) {
      for(auto node : nodes_)
         node->build_contribution_map(ancestor.nodes_);
   }

//private:
   std::vector<SingleNode<T>*> nodes_;
};

} /* namespace ics */
} /* namespace spral */
