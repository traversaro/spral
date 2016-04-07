#pragma once

namespace spral {
namespace ics {

template<typename T>
class SingleNode;
template<typename T>
class MultiNode;

/** Virtual base class for NodeToNodeMap and NodeToChunkMap */
template <typename T>
class MapBase {
public:
   virtual
   void apply(T *lval, T const* contrib, int ldcontrib, int* map) const=0;
};

template <typename T>
class NodeToNodeMap: public MapBase<T> {
public:
   NodeToNodeMap(SingleNode<T> const& ancestor, int const* row_start, int const* row_end, int offset)
   : ancestor_(ancestor), row_start_(row_start), row_end_(row_end), offset_(offset)
   {}

   /** Adds the contribution in contrib as per list (of rows).
    *  Uses map as workspace.
    *  \returns Number of columns found relevant. */
   void apply(T *lval, T const* contrib, int ldcontrib, int* map) const {
      contrib += offset_*(ldcontrib+1);
      T *lptr = &lval[ancestor_.loffset_];
      ancestor_.node_.construct_row_map(map);
      for(auto src_col=row_start_; src_col!=row_end_; ++src_col) {
         int cidx = std::distance(row_start_, src_col);
         if(! ancestor_.node_.contains_column(*src_col) )
            return; // done
         int col = map[ *src_col ]; // NB: Equal to *src_col - sptr[node]
         T const* src = &contrib[cidx*(ldcontrib+1)]; // Start on diagonal
         T *dest = &lptr[col * ancestor_.ldl_];
         for(auto src_row=src_col; src_row!=row_end_; ++src_row) {
            int row = map[ *src_row ];
            dest[row] += *(src++);
         }
      }
   }

   /** Return offset into node index list (i.e. how many rows to skip) */
   int get_offset() const { return offset_; }
private:
   int const* row_start_;
   int const* row_end_;
   long offset_;
   SingleNode<T> const& ancestor_;
};

template <typename T>
class NodeToMultiMap: public MapBase<T> {
public:
   template <typename U>
   NodeToMultiMap(U const& ancestor_nodes, SingleNode<T> const &from) {
      for(auto anode : ancestor_nodes) {
         NodeToNodeMap<T> *map = n2n_factory(from, *anode);
         if(map) maps_.push_back(map);
      }
   }
   ~NodeToMultiMap() {
      for(auto map : maps_) delete map;
   }

   void apply(T *lval, T const* contrib, int ldcontrib, int* work) const {
      for(auto map : maps_) {
         map->apply(lval, contrib, ldcontrib, work);
      }
   }
private:
   std::vector<NodeToNodeMap<T>*> maps_;
};

template <typename T>
class MultiMap {
public:
   MultiMap(SingleNode<T> const &ancestor, MultiNode<T> const& from)
   {
      int idx=0;
      for(auto node : from.nodes_) {
         NodeToNodeMap<T> *map = n2n_factory(*node, ancestor);
         long coffset = from.coffset_[idx];
         int ldcontrib = from.ldcontrib_[idx];
         if(map) maps_.emplace_back(coffset, ldcontrib, map);
         idx++;
      }
   }
   MultiMap(MultiNode<T> const &ancestor, MultiNode<T> const& from)
   {
      int idx=0;
      for(auto fnode : from.nodes_) {
         for(auto anode: ancestor.nodes_) {
            NodeToNodeMap<T> *map = n2n_factory(*fnode, *anode);
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
      map_tripple(long coffset, int ldcontrib, NodeToNodeMap<T> *map) : coffset(coffset), ldcontrib(ldcontrib), map(map) {}
      void apply(T *lval, T const* contrib, int* work) {
         map->apply(lval, &contrib[coffset], ldcontrib, work);
      }
      long coffset;
      int ldcontrib;
      NodeToNodeMap<T> *map;
   };
   std::vector<map_tripple> maps_;
};

/** Returns a map from a descendant from to an ancestor to. */
template <typename T>
static
NodeToNodeMap<T>* n2n_factory(SingleNode<T> const &from, SingleNode<T> const &to) {
   if(!to.is_ancestor_of(from)) return nullptr; // no upd
   int afirst = *(to.node_.row_begin());
   auto row_start = std::next(from.node_.row_begin(), from.n_);
   for(; *row_start < afirst && row_start!=from.node_.row_end(); ++row_start);
   if(!to.node_.contains_column(*row_start)) return nullptr; // no upd
   int offset = std::distance(from.node_.row_begin(), row_start) - from.n_;
   return new NodeToNodeMap<T>(to, row_start, from.node_.row_end(), offset);
}

} /* namespace ics */
} /* namespace spral */
