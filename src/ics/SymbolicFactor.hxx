#pragma once

#include <cstdlib>

namespace spral {
namespace ics {

class SymbolicFactor {

/* Only routine constructs symbolic factorization from matrix data */
public:
   SymbolicFactor (int n, int ptr[], int row[], int nemin);
   ~SymbolicFactor () {
      if(perm_) delete[] perm_;
   }

   long getFactorMemSize(void) const {
      return factor_mem_size_;
   }

   /* Information */
   const int nemin;
   long nfact, nflop;

protected:
   /* Core data */
   int n_;
   int nnodes_;
   int *perm_;
   long factor_mem_size_;

};

} /* namespace ics */
} /* namespace spral */
