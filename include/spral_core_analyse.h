#ifndef SPRAL_BASIC_ANALYSE_H
#define SPRAL_BASIC_ANALYSE_H

#ifdef __cplusplus
extern "C" {
#endif

int spral_core_analyse_basic_analyse(int n, int const ptr[], int const row[],
      int perm[], int *nnodes, int **sptr, int **sparent, long **rptr,
      int **rlist, int nemin, long *nfact, long *nflops, int base);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* SPRAL_BASIC_ANALYSE_H */
