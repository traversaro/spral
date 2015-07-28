#ifndef SPRAL_ND_H
#define SPRAL_ND_H

#include <stdbool.h>

struct spral_nd_options {
   int array_base;
   int print_level;
   int unit_diagnostics;
   int unit_error;
   int amd_call;
   int amd_switch1;
   double balance;
   int coarse_partition_method;
   bool find_supervariables;
   int max_improve_cycles;
   int matching;
   double max_reduction;
   double min_reduction;
   int partition_method;
   int refinement_band;
   bool remove_dense_rows;
   int stop_coarsening1;
   int stop_coarsening2;
};

struct spral_nd_inform {
   int flag;
   int dense;
   int nsuper;
   int nzsuper;
   int stat;
};

void spral_nd_default_options(struct spral_nd_options *nd);
void spral_nd_order(int mtx, int n, const int *ptr, const int *row, int *perm, const struct spral_nd_options *options, struct spral_nd_inform *inform);

#endif // SPRAL_ND_H