#include "NumericFactor.hxx"
#include "WorkspaceManager.hxx"

namespace spral {
namespace ics {

NumericFactor::NumericFactor(const SymbolicFactor &sf, const T aval[])
: sf_(sf), lmem(sf.get_factor_mem_size())
{
   WorkspaceManager memhandler(sf_.get_max_workspace_size());
   /* Iterate over chunks, factorizing */
   for(auto node = sf_.node_begin(); node != sf_.node_end(); ++node) {
      node->factor(
            aval, lmem.get_ptr(), sf_.get_ancestor_iterator(*node),
            sf_.get_ancestor_iterator_root(), memhandler
            );
   }
}

} /* namespace ics */
} /* namespace spral */
