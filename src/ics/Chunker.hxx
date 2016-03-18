#pragma once

#include <unordered_map>
#include <queue>
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
      /* Find number of children for each node; build list of nodes with no
       * children. */
      std::vector<int> nchild(tree.get_nnodes(), 0);
      std::queue<AssemblyTree::Node> ready;
      for(auto node=tree.begin(); node!=tree.end(); ++node) {
         if(node->has_parent())
            nchild[node->get_parent_node().idx]++;
         if(nchild[node->idx] == 0) ready.push(*node);
      }

      /* Loop over ready nodes adding them to chunks until we run out. */
      std::unordered_map< Coord, std::vector<AssemblyTree::Node>, CoordHash > buckets;
      auto check_itr = tree.leaf_first_begin();
      int next_chunk = 0;
      for(int nremain=tree.get_nnodes(); nremain>0; ) {
         /* Whilst we have ready nodes, keep assigning them */
         while(ready.size()) {
            auto node = ready.front(); ready.pop();
            auto& chunk = buckets[get_coord(node)];
            chunk.push_back(node);
            if(chunk.size() >= NUM_PER_CHUNK) {
               // Chunk is ready
               for(auto i=chunk.begin(); i!=chunk.end(); ++i) {
                  node_to_chunk_[i->idx] = next_chunk;
                  if(i->has_parent()) {
                     auto parent = i->get_parent_node();
                     if(--nchild[parent.idx] == 0)
                        ready.push(parent);
                  }
                  nremain--;
               }
               chunks_.push_back(chunk);
               next_chunk++;
               buckets.erase(get_coord(node));
            }
         }
         // FIXME: Try merging smaller nodes into bigger ones with large cols
         //        before we start closing off tiny ones
         /* If we run out, chase up from bottom of tree finding unassinged
          * nodes and closing them */
         // Find unallocated node
         for(; node_to_chunk_[check_itr->idx]!=-1; ++check_itr);
         if(check_itr==tree.leaf_first_end()) break;
         // Close it out
         auto& chunk = buckets[get_coord(*check_itr)];
         for(auto i=chunk.begin(); i!=chunk.end(); ++i) {
            node_to_chunk_[i->idx] = next_chunk;
            if(i->has_parent()) {
               auto parent = i->get_parent_node();
               if(--nchild[parent.idx] == 0)
                  ready.push(parent);
            }
            nremain--;
         }
         next_chunk++;
         chunks_.push_back(chunk);
         buckets.erase(get_coord(*check_itr));
      }
   }

   int operator[](AssemblyTree::Node const& node) const {
      return node_to_chunk_[node.idx];
   }
   int operator[](int idx) const {
      return node_to_chunk_[idx];
   }

   std::vector<AssemblyTree::Node> const& get_chunk(int idx) const {
      return chunks_[idx];
   }

   std::vector< std::vector<AssemblyTree::Node> >::const_iterator begin() const {
      return chunks_.begin();
   }

   std::vector< std::vector<AssemblyTree::Node> >::const_iterator end() const {
      return chunks_.end();
   }

private:
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
   std::vector< std::vector<AssemblyTree::Node> > chunks_;
};

} /* namespace ics */
} /* namespace spral */
