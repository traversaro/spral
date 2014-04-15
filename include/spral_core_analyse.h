#ifndef SPRAL_BASIC_ANALYSE_H
#define SPRAL_BASIC_ANALYSE_H

int spral_core_analyse_basic_analyse(int n, const int ptr[], const int row[],
      int perm[], int *nnodes, int **sptr, int **sparent, long **rptr,
      int **rlist, int nemin, long *nfact, long *nflops, int base);

#endif /* SPRAL_BASIC_ANALYSE_H */
