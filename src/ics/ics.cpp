extern "C" {
#include "spral.h"
}

#include <cstddef>
#include <stdexcept>
using namespace std;

namespace spral { namespace ics {
   class symbolic {
      /* Core data */
      protected:
         int n;
         int nnodes;
         int *perm;
         int *sptr;
         int *sparent;
         long *rptr;
         int *rlist;

      /* Information */
      public:
         int nemin;
         long nfact, nflops;

      /* Only routine constructs symbolic factorization from matrix data */
      public:
         symbolic (int n, int ptr[], int row[], int nemin) :
               n(n), nnodes(0), perm(NULL), sptr(NULL), sparent(NULL),
               rptr(NULL), rlist(NULL), nemin(nemin), nfact(0), nflops(0) {

            /* Perform METIS ordering */
            perm = new int[n];
            int *invp = new int[n];
            int flag = spral_metis_order(n, ptr, row, perm, invp, 0);
            delete[] invp;
            if(flag!=0) throw runtime_error("spral_metis_order() failed");

            /* Find assembly tree */
            flag = spral_core_analyse_basic_analyse(n, ptr, row, perm, &nnodes,
               &sptr, &sparent, &rptr, &rlist, nemin, &nfact, &nflops, 0);
            if(flag!=0) throw runtime_error("spral_core_analyse_basic_analyse() failed");
         }
         ~symbolic () {
         }
   };
}}
