#pragma once

#include "blas_templates.hxx"

#include "AssemblyTree.hxx"
#include "ICSExceptions.hxx"
#include "Maps.hxx"
#include "SimdVec.hxx"
#include "WorkspaceManager.hxx"

namespace spral {
namespace ics {

template <typename T>
class SingleNode {
public:
   explicit SingleNode(AssemblyTree::Node const& node)
   : node_(node), m_(node.get_nrow()), n_(node.get_ncol()), loffset_(0),
     ldl_(0)
   {}
   ~SingleNode() {
      for(auto cmap : contribution_map_)
         delete cmap;
   }

   void set_memloc(long loffset, int ldl) {
      loffset_ = loffset;
      ldl_ = ldl;
   }

   /** Builds a contribution map for a single node. Perform no-op if not an actual ancestor. */
   void build_contribution_map(SingleNode<T> const& ancestor) {
      NodeToNodeMap<T> *map = n2n_factory(*this, ancestor);
      if(map) contribution_map_.push_back(map);
   }

   /** Builds a contribution map for a chunk. */
   void build_contribution_map(MultiNode<T> const& ancestor) {
      contribution_map_.push_back(
            new NodeToMultiMap<T>(ancestor, *this)
            );
   }

   bool has_parent() const {
      return node_.has_parent();
   }
   int get_parent_idx() const {
      return node_.get_parent_node().idx;
   }
   /** Return a representative index in range [0,nnodes-1] */
   int get_idx() const {
      return node_.idx;
   }

   bool is_ancestor_of(SingleNode<T> const& node) const {
      return node_.has_descendant(node.node_);
   }

   void print(T const*lval) const {
      printf("NODE %d is %d x %d (parent %d)\n", node_.idx, m_, n_,
            node_.get_parent_node().idx);
      int idx=0;
      for(auto row=node_.row_begin(); row!=node_.row_end(); ++row, ++idx) {
         printf("%d:", *row);
         for(int col=0; col<n_; ++col)
            printf(" %e", lval[loffset_ + col*ldl_ + idx]);
         printf("\n");
      }
      printf("\n");
   }

   void factor(T const* aval, T* lval, WorkspaceManager &memhandler) const {
      //printf("Factor %d: %dx%d\n", node_.idx, m_, n_);

      int ldcontrib = m_ - n_;
      T *contrib =
         (ldcontrib > 0) ? memhandler.get<T>((m_-n_)*((long) ldcontrib))
                         : nullptr;
      int *map =
         (ldcontrib > 0) ? memhandler.get<int>(node_.get_max_tree_row_index())
                         : nullptr;

      factor_local(aval, lval, contrib, ldcontrib);

      /* Distribute elements of generated element to ancestors */
      for(auto n2n_map = contribution_map_.begin(); n2n_map != contribution_map_.end(); ++n2n_map) {
         (*n2n_map)->apply(lval, contrib, ldcontrib, map);
      }

      /* Release workspace */
      if(map) memhandler.release<int>(map, node_.get_max_tree_row_index());
      if(contrib) memhandler.release<T>(contrib, (m_-n_)*((long) ldcontrib));
   }

   /* Perform factor operations local to this node */
   void factor_local(T const *aval, T *lval, T *contrib, int ldcontrib) const {
      /* Setup working pointers */
      T *ldiag = &lval[loffset_];
      T *lrect = &lval[loffset_ + n_];

      /* Add A */
      for(auto ait=node_.a_to_l_begin(); ait!=node_.a_to_l_end(); ++ait) {
         int dest = ait->dest_col*ldl_ + ait->dest_row;
         ldiag[dest] += aval[ait->src];
      }

      /* Factorize diagonal block */
      int info = potrf<T>('L', n_, ldiag, ldl_);
      if(info) throw NotPositiveDefiniteException(info);

      /* Only work with generated element if it exists! */
      if(m_ - n_ > 0) { 
         /* Apply to block below diagonal */
         trsm<T>('R', 'L', 'T', 'N', m_-n_, n_, 1.0, ldiag, ldl_, lrect, ldl_);

         /* Form generated element into workspace */
         syrk<T>('L', 'N', m_-n_, n_, -1.0, lrect, ldl_, 0.0, contrib,
               ldcontrib);
      }
   }

   void forward_solve(int nrhs, T* x, int ldx, T const* lval,
         WorkspaceManager &memhandler) const {
      /* Perform solve with diagonal block */
      T const* ldiag = &lval[loffset_];
      T* xdiag = &x[*node_.row_begin()];
      trsm<T>('L', 'L', 'N', 'N', n_, nrhs, 1.0, ldiag, ldl_, xdiag, ldx);

      if(m_-n_>0) { /* Entries below diagonal block exist */
         /* Get some workspace to calculate result in */
         int ldxlocal = m_ - n_;
         T *xlocal = memhandler.get<T>(nrhs*ldxlocal);

         /* Calculate contributions to ancestors */
         T const* lrect = &ldiag[n_];
         gemm<T>('N', 'N', m_-n_, nrhs, n_, -1.0, lrect, ldl_, xdiag, ldx, 0.0,
               xlocal, ldxlocal);

         /* Distribute contributions */
         int idx=0;
         for(auto row=node_.row_begin()+n_; row!=node_.row_end(); ++row, ++idx)
            for(int r=0; r<nrhs; r++)
               x[r*ldx + *row] += xlocal[r*ldxlocal + idx];

         /* Release workspace */
         memhandler.release<T>(xlocal, nrhs*ldxlocal);
      }
   }

   void backward_solve(int nrhs, T* x, int ldx, T const* lval,
         WorkspaceManager &memhandler) const {
      /* Establish useful pointers */
      T const* ldiag = &lval[loffset_];
      T* xdiag = &x[*node_.row_begin()];

      if(m_-n_>0) { /* Entries below diagonal block exist */
         /* Get some workspace to calculate result in */
         int ldxlocal = m_ - n_;
         T *xlocal = memhandler.get<T>(nrhs*ldxlocal);

         /* Gather values from ancestors */
         int idx=0;
         for(auto row=node_.row_begin()+n_; row!=node_.row_end(); ++row, ++idx)
            for(int r=0; r<nrhs; r++)
               xlocal[r*ldxlocal + idx] = x[r*ldx + *row];

         /* Apply update to xdiag[] */
         T const* lrect = &ldiag[n_];
         gemm<T>('T', 'N', n_, nrhs, m_-n_, -1.0, lrect, ldl_, xlocal,
               ldxlocal, 1.0, xdiag, ldx);

         /* Release workspace */
         memhandler.release<T>(xlocal, nrhs*ldxlocal);
      }

      /* Perform solve with diagonal block */
      trsm<T>('L', 'L', 'T', 'N', n_, nrhs, 1.0, ldiag, ldl_, xdiag, ldx);
   }

private:
   /* Friends */
   friend class NodeToNodeMap<T>;
   friend NodeToNodeMap<T>* n2n_factory<T>(SingleNode<T> const&, SingleNode<T> const&);

   /* Members */
   AssemblyTree::Node const node_;
   int const m_;
   int const n_;
   long loffset_;
   int ldl_;
   std::vector< MapBase<T>* > contribution_map_;
};

} /* namespace ics */
} /* namespace spral */
