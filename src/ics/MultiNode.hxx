#pragma once

#include "SingleNode.hxx"

#include <vector>

namespace spral {
namespace ics {

template<typename T>
class MultiNode {
   class MultiMap {
   public:
      MultiMap(SingleNode<T> const &ancestor, MultiNode const& from)
      {
         int idx=0;
         for(auto node : from.nodes_) {
            typename SingleNode<T>::NodeToNodeMap *map =
               SingleNode<T>::n2n_factory(*node, ancestor);
            long coffset = from.coffset_[idx];
            int ldcontrib = from.ldcontrib_[idx];
            if(map) maps_.emplace_back(coffset, ldcontrib, map);
            idx++;
         }
      }
      MultiMap(MultiNode const &ancestor, MultiNode const& from)
      {
         int idx=0;
         for(auto fnode : from.nodes_) {
            for(auto anode: ancestor.nodes_) {
               typename SingleNode<T>::NodeToNodeMap *map =
                  SingleNode<T>::n2n_factory(*fnode, *anode);
               long coffset = from.coffset_[idx];
               int ldcontrib = from.ldcontrib_[idx];
               if(map) maps_.emplace_back(coffset, ldcontrib, map);
            }
            idx++;
         }
      }
      ~MultiMap() {
         for(auto map : maps_)
            delete map.map;
      }
      void apply(T *lval, T const* contrib, int* work) {
         for(auto map : maps_)
            map.apply(lval, contrib, work);
      }
   private:
      struct map_tripple {
         map_tripple(long coffset, int ldcontrib, typename SingleNode<T>::NodeToNodeMap *map) : coffset(coffset), ldcontrib(ldcontrib), map(map) {}
         void apply(T *lval, T const* contrib, int* work) {
            map->apply(lval, &contrib[coffset], ldcontrib, work);
         }
         long coffset;
         int ldcontrib;
         typename SingleNode<T>::NodeToNodeMap *map;
      };
      std::vector<map_tripple> maps_;
   };
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
      maps_.push_back( new MultiMap(ancestor, *this) );
   }

   void build_contribution_map(MultiNode<T> const& ancestor) {
      maps_.push_back( new MultiMap(ancestor, *this) );
   }

//private:
   std::vector<SingleNode<T>*> nodes_;
private:
   int matrix_n_;
   long contrib_size_;
   std::vector<long> coffset_; //< Offset into contrib for each node
   std::vector<int> ldcontrib_; //< Leading dimensions of contrib blocks
   std::vector<MultiMap*> maps_;
};

} /* namespace ics */
} /* namespace spral */
