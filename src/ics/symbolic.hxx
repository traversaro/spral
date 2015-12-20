#pragma once

#include <cstdlib>

namespace spral {
namespace ics {

class symbolic {

/* Only routine constructs symbolic factorization from matrix data */
public:
   symbolic (int n, int ptr[], int row[], int nemin);
   ~symbolic () {
      if(perm_) delete[] perm_;
      if(sptr_) free(sptr_); // Allocated by malloc in C interface fn
      if(sparent_) free(sparent_); // Allocated by malloc in C interface fn
      if(rptr_) free(rptr_); // Allocated by malloc in C interface fn
      if(rlist_) free(rlist_); // Allocated by malloc in C interface fn
   }

   /* Information */
   const int nemin;
   long nfact, nflop;

protected:
   /* Core data */
   int n_;
   int nnodes_;
   int *perm_;
   int *sptr_;
   int *sparent_;
   long *rptr_;
   int *rlist_;

};

} /* namespace ics */
} /* namespace spral */
