#include "NumericFactor.hxx"

namespace spral {
namespace ics {

NumericFactor::NumericFactor(const SymbolicFactor &sf, const T aval[])
: sf_(sf), lmem(sf.getFactorMemSize())
{
   /* Iterate over chunks, factorizing */
   /*for(auto chunk = sf.tree_start(); node != sf.tree_end(); ++node) {
      chunk.factor(aval, lmem.ptr);
   }*/
}

} /* namespace ics */
} /* namespace spral */
