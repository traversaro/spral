#pragma once

#include "SymbolicFactor.hxx"

namespace spral {
namespace ics {

class NumericFactor {
public:
   /** Perform a factorization as part of construction */
   NumericFactor(const SymbolicFactor &sf, int n, const int ptr[],
         const int row[], const double val[]);

   /** Perform a solve with factors */
   void solve(int nrhs, double x[]);

private:
   const SymbolicFactor &sf_;
};

} /* namespace ics */
} /* namespace spral */
