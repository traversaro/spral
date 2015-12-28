#pragma once

#include "AlignedCallocBlock.hxx"
#include "SymbolicFactor.hxx"

namespace spral {
namespace ics {

class NumericFactor {
   typedef double T;
public:
   /** Perform a factorization as part of construction */
   NumericFactor(SymbolicFactor const& sf, T const aval[]);

   /** Perform an in-place solve with factors (single rhs version) */
   void solve(T x[]) const {
      solve(1, x, sf_.get_n()); // Call multiple rhs version with nrhs=1
   }

   /** Perform an in-place solve with factors (multiple rhs version) */
   void solve(int nrhs, T x[], int ldx) const;

   /** Permutes rhs to elimination order */
   void permute_rhs_to_elim(int nrhs, T const x[], int ldx, T xperm[],
         int ldxperm) const;

   /** Permutes rhs from elimination order */
   void permute_rhs_from_elim(int nrhs, T const xperm[], int ldxperm, T x[],
         int ldx) const;

   /** Performs forward solve in-place on rhs in elim order */
   void forward_solve(int nrhs, T x[], int ldx) const;

   /** Performs bacward solve in-place on rhs in elim order */
   void backward_solve(int nrhs, T x[], int ldx) const;

private:
   const SymbolicFactor &sf_;
   AlignedCallocBlock<T> lmem;
};

} /* namespace ics */
} /* namespace spral */
