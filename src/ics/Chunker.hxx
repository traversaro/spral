#pragma once

#include <unordered_map>
#include <vector>

#include "AssemblyTree.hxx"

namespace spral {
namespace ics {

class Chunker {
private:
   static int const ROW_CHUNK_SIZE = 4; //< waste up to ROW_CHUNK_SIZE-1 rows
   static int const COL_CHUNK_SIZE = 1; //< waste up to COL_CHUNK_SIZE-1 cols
   static int const MAX_NROW = 20*ROW_CHUNK_SIZE; //< maximum number of rows
   static int const MAX_NCOL = 8*COL_CHUNK_SIZE; //< maximum number of cols
   static int const NUM_PER_CHUNK = 4; //< Number of nodes to group per chunk
   typedef std::pair<int,int> Coord;
public:
   Chunker(AssemblyTree const& tree)
   : node_to_chunk_(tree.get_nnodes(), -1)
   {
      /* Find depth for each node (distance from root) */
      int maxdepth = 0;
      std::vector<int> depth(tree.get_nnodes(), 0);
      for(auto node=tree.rbegin(); node!=tree.rend(); ++node) {
         int parent = node->get_parent_node().idx;
         if(parent >= tree.get_nnodes()) continue; // No parent for roots
         depth[node->idx] = depth[parent] + 1;
         maxdepth = std::max(maxdepth, depth[node->idx]);
      }

      /* Sort nodes into buckets based on their sizes */
      std::unordered_multimap< Coord, AssemblyTree::Node, CoordHash > buckets;
      for(auto node=tree.rbegin(); node!=tree.rend(); ++node)
         buckets.emplace( get_coord(*node), *node );

      /* Starting with largest buckets, try and fill out chunks */
      int next_chunk=0;
      for(int c=MAX_NCOL-1; c>=0; --c) {
         for(int r=MAX_NROW-1; r>=0; --r) {
            while(buckets.count( Coord(r,c) ) > 0) {
               auto node = buckets.find( Coord(r, c) ); // Any will do!
               std::vector<AssemblyTree::Node> chunk = gather_chunk(buckets, node);
               for(auto i = chunk.begin(); i!=chunk.end(); ++i)
                  node_to_chunk_[i->idx] = next_chunk;
               ++next_chunk;
            }
         }
      }
   }

   int operator[](AssemblyTree::Node const& node) const {
      return node_to_chunk_[node.idx];
   }

private:
   /** Given a first element, gather a chunk to fill it out. Update buckets
    *  and node_to_chunk to refelect this. */
   template<typename map_type>
   std::vector<AssemblyTree::Node> gather_chunk(map_type &buckets, typename map_type::iterator first_node_itr) {
      /* Init chunk membership list to supplied entry */
      std::vector<AssemblyTree::Node> chunk;
      chunk.push_back(first_node_itr->second);
      int const rorig = first_node_itr->first.first;
      int const c = first_node_itr->first.second;
      buckets.erase(first_node_itr);
      if(chunk.size() >= NUM_PER_CHUNK) return chunk;
      
      /* NB we can reduce #rows, but not #cols */
      for(int r=rorig; r>=0; --r) {
         auto range = buckets.equal_range( Coord(r, c) );
         for(auto i=range.first; i!=range.second;) {
            AssemblyTree::Node const& candidate = i->second;
            // check if candidate is independent of rest of chunk
            bool valid = true;
            for(auto node=chunk.begin(); node!=chunk.end(); ++node) {
               if(candidate.has_descendant(*node)) valid = false;
               if(node->has_descendant(candidate)) valid = false;
            }
            if(valid) {
               chunk.push_back(candidate);
               i = buckets.erase(i);
               if(chunk.size() >= NUM_PER_CHUNK) return chunk;
            } else {
               ++i; // NB: only increment if we've not erased [which increments]
            }
         }
      }

      /* Otherwise, return non-full chunk! */
      return chunk;
   }

   Coord get_coord(AssemblyTree::Node const& node) {
      return Coord(
         std::min( (node.get_nrow()-1) / ROW_CHUNK_SIZE, MAX_NROW ),
         std::min( (node.get_ncol()-1) / COL_CHUNK_SIZE, MAX_NCOL )
         );
   }

   class CoordHash {
   public:
      std::size_t operator() (Coord const& coord) const {
         return coord.first + coord.second;
      }
   };

   std::vector<int> node_to_chunk_;
};

} /* namespace ics */
} /* namespace spral */
