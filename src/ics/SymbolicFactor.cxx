#include "SymbolicFactor.hxx"

#include <stdexcept>

#include "AssemblyTree.hxx"
#include "Chunk.hxx"
#include "Chunker.hxx"
#include "spral_metis_wrapper.h"

namespace spral {
namespace ics {

/* Constructs symbolic factorization from matrix data */
SymbolicFactor::SymbolicFactor (int n, int ptr[], int row[], int nemin)
: nemin(nemin), n_(n), perm_(nullptr), factor_mem_size_(0),
  max_workspace_size_(0), tree_(n)
{

   /* Perform METIS ordering */
   perm_ = new int[n];
   int *invp = new int[n];
   int flag = spral_metis_order(n, ptr, row, perm_, invp, 0);
   delete[] invp;
   if(flag)
      throw std::runtime_error("spral_metis_order() failed");

   /* Construct AssemblyTree */
   tree_.construct_tree(ptr, row, perm_, nemin);

   /* Construct list of nodes */
   factor_mem_size_ = 0;
   int max_contrib_size = 0;
   for(auto node=tree_.begin(); node!=tree_.end(); ++node) {
      int m = node->get_nrow();
      int n = node->get_ncol();
      nodes_.push_back(Node<double>(*node, factor_mem_size_, m));
      factor_mem_size_ += m*((long) n);
      max_contrib_size = std::max(max_contrib_size, m-n);
   }
   max_workspace_size_ = max_contrib_size*max_contrib_size*sizeof(double) +
      n_*sizeof(int);

   /* Now we know where nodes are in memory, build map */
   int *map = new int[n_];
   for(auto node=nodes_.begin(); node!=nodes_.end(); ++node)
      node->build_contribution_map(
            get_ancestor_iterator(*node), get_ancestor_iterator_root(), map
            );
   delete[] map;

#if 0
   /* Construct list of chunks */
   Chunker chunker(tree_);
   for(auto node=tree_.begin(); node!=tree_.end(); ++node) {
      printf("Node %d (%d x %d) in chunk %d\n", node->idx, node->get_nrow(),
            node->get_ncol(), chunker[*node]);
   }
#endif

#if 0
   /* Construct chunk buckets */
   const int MAXROW=50;
   const int MAXCOL=8;
   int clen[MAXROW+1][MAXCOL+1];
   for(int i=0; i<MAXROW+1; i++)
   for(int j=0; j<MAXCOL+1; j++)
      clen[i][j] = 0;
   for(auto nitr=tree_.leaf_first_begin(); nitr!=tree_.leaf_first_end(); ++nitr) {

      const AssemblyTree::Node node = *nitr;
      int i = (node.get_nrow()-1)/4;
      if(i>MAXROW) i = MAXROW;
      int j = node.get_ncol()-1;
      if(j>MAXCOL) j = MAXCOL;
      if(i==MAXROW || j==MAXCOL)
         printf("Node %d is %d x %d\n", node.idx, node.get_nrow(), node.get_ncol());
      clen[i][j]++;
   }
   printf("Buckets:\n  ");
   for(int i=0; i<MAXCOL; i++) printf(" %4d", i+1);
   printf("   >%d\n", MAXCOL);
   for(int i=0; i<MAXROW+1; i++) {
      printf("%2d", i+1);
      for(int j=0; j<MAXCOL+1; j++)
         printf(" %4d", clen[i][j]);
      printf("\n");
   }
#endif

}

} /* namespace ics */
} /* namespace spral */
