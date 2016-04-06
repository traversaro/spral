#pragma once

#include <cstdlib>

#include "AssemblyTree.hxx"
#include "SingleNode.hxx"

namespace spral {
namespace ics {

class SymbolicFactor {
   typedef double T;
public:

   class Chunk {
   public:
      /** Constructor */
      Chunk(SymbolicFactor const& sf)
      : sf_(sf), sn_(nullptr)
      {}

      /** Destructor */
      ~Chunk() {
         if(sn_) {
            delete sn_;
         } else {
            for(auto node=nodes_.begin(); node!=nodes_.end(); ++node)
               delete (*node);
         }
      }

      SingleNode<T>* emplace_node(AssemblyTree::Node const& node) {
         SingleNode<T> *sn = new SingleNode<T>(node);
         add_node(sn);
         return sn;
      }

      /** Add a node to the chunk */
      void add_node(SingleNode<T> *sn) {
         if(nodes_.size() > 0) {
            nodes_.push_back(sn);
         }
         else if(sn_) {
            nodes_.push_back(sn_);
            sn_ = nullptr;
            nodes_.push_back(sn);
         }
         else {
            sn_ = sn;
         }
      }

      /** Adds parent/child relation between chunks */
      friend
      void add_relation(Chunk &child, Chunk &parent) {
         child.parents_.push_back(&parent);
         parent.children_.push_back(&child);
      }

      /** Return true if this chunk has a parent */
      bool has_parent() const {
         return (parents_.size() > 0);
      }

      /** Factorize all nodes in this chunk */
      void factor(T const* aval, T* lval, WorkspaceManager &memhandler) const {
         if(sn_)
            sn_->factor(aval, lval, memhandler);
         else
            for(auto node=nodes_.begin(); node<nodes_.end(); ++node)
               (*node)->factor(aval, lval, memhandler);
      }

      /** Perform forward solve for all nodes in this chunk */
      void forward_solve(int nrhs, T* x, int ldx, T const* lval,
            WorkspaceManager &memhandler) const {
         if(sn_)
            sn_->forward_solve(nrhs, x, ldx, lval, memhandler);
         else
            for(auto node=nodes_.begin(); node<nodes_.end(); ++node)
               (*node)->forward_solve(nrhs, x, ldx, lval, memhandler);
      }

      /** Perform backwards solve for all nodes in this chunk */
      void backward_solve(int nrhs, T* x, int ldx, T const* lval,
            WorkspaceManager &memhandler) const {
         if(sn_)
            sn_->backward_solve(nrhs, x, ldx, lval, memhandler);
         else
            for(auto node=nodes_.begin(); node<nodes_.end(); ++node)
               (*node)->backward_solve(nrhs, x, ldx, lval, memhandler);
      }

      /** Print the data in this chunk */
      void print(T const* lval) const {
         printf("CHUNK\n");
         if(sn_)
            sn_->print(lval);
         else
            for(auto node=nodes_.begin(); node<nodes_.end(); ++node)
               (*node)->print(lval);
      }

      /** Build contribution map */
      void build_contribution_map(int nnodes) {
         if(!has_parent()) return; // no parent, no map

         std::vector<int> seen(nnodes, false);
         std::vector<Chunk const*> stack;
         for(auto pchunk=parents_.begin(); pchunk!=parents_.end(); ++pchunk) {
            int pfn = (*pchunk)->sn_ ? (*pchunk)->sn_->get_idx() : (*pchunk)->nodes_.front()->get_idx();
            seen[pfn] = true;
            stack.push_back(*pchunk);
         }
         
         while(stack.size()) {
            auto chunk = stack.back(); stack.pop_back();
            int first_node = chunk->sn_ ? chunk->sn_->get_idx() : chunk->nodes_.front()->get_idx();
            for(auto pchunk=chunk->parents_.begin(); pchunk!=chunk->parents_.end(); ++pchunk) {
               int pfn = (*pchunk)->sn_ ? (*pchunk)->sn_->get_idx() : (*pchunk)->nodes_.front()->get_idx();
               if(!seen[pfn]) {
                  seen[pfn] = true;
                  stack.push_back(*pchunk);
               }
            }
            if(sn_) {
               if(chunk->sn_) {
                  // Both this chunk and ancestor are single nodes
                  sn_->build_contribution_map(*chunk->sn_);
               } else {
                  // This chunk is a single node, ancestor is multi
                  for(auto pnode=chunk->nodes_.begin(); pnode!=chunk->nodes_.end(); ++pnode) {
                     sn_->build_contribution_map(**pnode);
                  }
               }
            } else {
               if(chunk->sn_) {
                  // This chunk is multi but ancestor is single node
                  for(auto node=nodes_.begin(); node!=nodes_.end(); ++node) {
                     (*node)->build_contribution_map(*(chunk->sn_));
                  }
               } else {
                  // Both this chunk and ancestor are multi
                  for(auto node=nodes_.begin(); node!=nodes_.end(); ++node) {
                     for(auto pnode=chunk->nodes_.begin(); pnode!=chunk->nodes_.end(); ++pnode) {
                        (*node)->build_contribution_map(**pnode);
                     }
                  }
               }
            }
         }
      }

   private:
      SymbolicFactor const& sf_;
      std::vector<Chunk *> parents_;
      std::vector<Chunk *> children_;
      SingleNode<T> *sn_;
      std::vector<SingleNode<T>*> nodes_;
   };

   /** Performs a symbolic factorization as part of the construction */
   SymbolicFactor (int n, int ptr[], int row[], int nemin);
   ~SymbolicFactor () {
      if(perm_) delete[] perm_;
   }

   /** Returns length of factor array to be used */
   long get_factor_mem_size(void) const {
      return factor_mem_size_;
   }

   /** Returns maximum size of workspace to be allocated */
   long get_max_workspace_size(void) const {
      return max_workspace_size_;
   }

   /** Returns iterator to beginning of chunk list */
   std::vector< Chunk >::const_iterator chunk_begin(void) const {
      return chunks_.cbegin();
   }
   /** Returns iterator to end of chunk list */
   std::vector< Chunk >::const_iterator chunk_end(void) const {
      return chunks_.cend();
   }

   /** Returns iterator to reverse beginning of chunk list */
   std::vector< Chunk >::const_reverse_iterator chunk_rbegin(void) const {
      return chunks_.crbegin();
   }
   /** Returns iterator to reverse end of chunk list */
   std::vector< Chunk >::const_reverse_iterator chunk_rend(void) const {
      return chunks_.crend();
   }

   /* Information */
   const int nemin;
   long get_nfact() const { return tree_.get_nfact(); }
   long get_nflop() const { return tree_.get_nflop(); }
   int get_n() const { return n_; }
   int get_nnodes() const { return tree_.get_nnodes(); }
   int const* get_perm() const { return perm_; }
   int get_nchunks() const { return chunks_.size(); }

private:
   /* Core data */
   int n_;
   int *perm_;
   long factor_mem_size_;
   long max_workspace_size_;
   AssemblyTree tree_;
   std::vector< Chunk > chunks_;
};

} /* namespace ics */
} /* namespace spral */
