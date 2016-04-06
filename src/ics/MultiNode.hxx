#pragma once

#include "SingleNode.hxx"

#include <vector>

namespace spral {
namespace ics {

template<typename T>
class MultiNode {
public:
   template <typename node_itr>
   MultiNode(node_itr const& nbegin, node_itr const& nend, long &loffset, long &max_work_size) {
      for(auto node=nbegin; node!=nend; ++node) {
         int m = node.get_nrow();
         long n = node.get_ncol();
         nodes_.push_back(SingleNode<T>(node));
         nodes_.back().set_memloc(loffset, m);
         loffset += m*n;
         max_work_size = std::max(max_work_size, (m-n)*(m-n));
      }
   }

private:
   std::vector<SingleNode<T>> nodes_;
};

} /* namespace ics */
} /* namespace spral */
