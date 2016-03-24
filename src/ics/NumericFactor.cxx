#include "NumericFactor.hxx"
#include "WorkspaceManager.hxx"

namespace spral {
namespace ics {

NumericFactor::NumericFactor(const SymbolicFactor &sf, const T aval[])
: sf_(sf), lmem_(sf.get_factor_mem_size())
{
   WorkspaceManager memhandler(sf_.get_max_workspace_size());
   std::vector<int> child_complete(sf_.get_nchunks(), 0);
   /* Iterate over chunks, factorizing */
   for(auto ci = sf_.chunk_begin(); ci != sf_.chunk_end(); ++ci) {
      auto const& chunk = *(*ci);
      if(chunk.get_nchild() == 0) {
         chunk.factor(aval, lmem_.get_ptr(), memhandler);
         /*if(chunk.has_parent()) {
            auto const& parent = chunk.get_parent();
            child_complete[parent.get_idx()]++;
            if(child_complete[parent.get_idx()] != parent.get_nchild())
               printf("Ready to run %d\n", parent.get_idx());
         }*/
      }
   }
}

void NumericFactor::solve(int nrhs, T x[], int ldx) const {
   int ldxperm = sf_.get_n();
   T *xperm = new T[nrhs*ldxperm];
   permute_rhs_to_elim(nrhs, x, ldx, xperm, ldxperm);
   forward_solve(nrhs, xperm, ldxperm);
   backward_solve(nrhs, xperm, ldxperm);
   permute_rhs_from_elim(nrhs, xperm, ldxperm, x, ldx);
   delete[] xperm;
}

void NumericFactor::permute_rhs_to_elim(int nrhs, T const x[], int ldx,
      T xperm[], int ldxperm) const {
   int const* perm = sf_.get_perm();
   for(int i=0; i<sf_.get_n(); i++) {
      int j = perm[i];
      for(int r=0; r<nrhs; r++)
         xperm[r*ldxperm + j] = x[r*ldx + i];
   }
}

void NumericFactor::permute_rhs_from_elim(int nrhs, T const xperm[],
      int ldxperm, T x[], int ldx) const {
   int const* perm = sf_.get_perm();
   for(int i=0; i<sf_.get_n(); i++) {
      int j = perm[i];
      for(int r=0; r<nrhs; r++)
         x[r*ldx + i] = xperm[r*ldxperm + j];
   }
}

void NumericFactor::forward_solve(int nrhs, T x[], int ldx) const {
   WorkspaceManager memhandler(nrhs*sf_.get_n()*sizeof(double));
   /* Iterate over node forward, doing forward solve  */
   for(auto node = sf_.node_begin(); node != sf_.node_end(); ++node) {
      node->forward_solve(
            nrhs, x, ldx, lmem_.get_ptr(), sf_.get_ancestor_iterator(*node),
            sf_.get_ancestor_iterator_root(), memhandler
            );
   }
}

void NumericFactor::backward_solve(int nrhs, T x[], int ldx) const {
   WorkspaceManager memhandler(nrhs*sf_.get_n()*sizeof(double));
   /* Iterate over nodes backwards, doing backward solve */
   for(auto node = sf_.node_rbegin(); node != sf_.node_rend(); ++node) {
      node->backward_solve(
            nrhs, x, ldx, lmem_.get_ptr(), sf_.get_ancestor_iterator(*node),
            sf_.get_ancestor_iterator_root(), memhandler
            );
   }
}

} /* namespace ics */
} /* namespace spral */
