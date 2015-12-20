#pragma once

namespace spral {
namespace ics {

class symbolic {

/* Only routine constructs symbolic factorization from matrix data */
public:
   symbolic (int n, int ptr[], int row[], int nemin);
   ~symbolic () {
      if(perm_) delete[] perm_;
   }

   /* Information */
   int nemin;
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
