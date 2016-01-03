#include <boost/program_options.hpp>
#include <iostream>
#include <time.h>

#include "spral.h"
#include "src/ics/SymbolicFactor.hxx"
#include "src/ics/NumericFactor.hxx"

float tdiff(const struct timespec &t1, const struct timespec &t2) {
   return t2.tv_sec - t1.tv_sec + 1e-9*(t2.tv_nsec - t1.tv_nsec);
}

struct DriverOptions {
   bool print_matrix;
   bool print_factors;
   bool force_posdef;
   int nemin;

   DriverOptions(int argc, char *const * argv) {
      boost::program_options::options_description desc("Allowed options");
      desc.add_options()
         ("help", "produce help message")
         ("print-matrix", "print input matrix")
         ("print-factors", "print numeric factors")
         ("force-posdef", "randomly generate values for posdef matrix")
         ("nemin", boost::program_options::value<int>(&nemin)->default_value(8),
            "supernode amalgamation parameter")
         ;
      boost::program_options::variables_map vm;
      boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);
      boost::program_options::notify(vm);

      /* Print help and exit if requested */
      if( vm.count("help") ) {
         std::cout << desc << std::endl;
         exit(1);
      }

      /* Set values that require present/not-present */
      print_matrix = ( vm.count("print-matrix") );
      print_factors = ( vm.count("print-factors") );
      force_posdef = ( vm.count("force-posdef") );

      /* Feedback non-obvious settings to user */
      std::cout << "nemin = " << nemin << std::endl;
   }
};


template <typename T>
void spmv(int n, int const* ptr, int const* row, T const* val, T const* x, T *y){
   for(int i=0; i<n; ++i) y[i] = 0.0;
   for(int i=0; i<n; ++i) {
      for(int j=ptr[i]; j<ptr[i+1]; ++j) {
         int k = row[j];
         y[k] += val[j] * x[i];
         if(i==k) continue;
         y[i] += val[j] * x[k];
      }
   }
}

template <typename T>
T matrix_inf_norm(int n, int const* ptr, int const* row, T const* val) {
   /* NB: inf norm is maximum row sum of abs values */
   T *row_sum = new T[n];
   for(int i=0; i<n; ++i) row_sum[i] = 0.0;
   for(int i=0; i<n; ++i) {
      for(int j=ptr[i]; j<ptr[i+1]; ++j) {
         int k = row[j];
         row_sum[k] += fabs(val[j]);
         if(i==k) continue;
         row_sum[i] += fabs(val[j]);
      }
   }

   T best = *std::max_element(row_sum, row_sum+n);
   delete[] row_sum;

   return best;
}

template <typename T>
T max_abs_value(T const* first, T const* end) {
   T bestv = 0.0;
   for(auto v=first; v!=end; ++v)
      bestv = std::max(bestv, fabs(*v));
   return bestv;
}

template <typename T>
T calc_bwd_error(int n, int const* ptr, int const* row, T const* val, T const* soln, T const* rhs) {
   /* Calculate residual */
   T *resid = new T[n];
   spmv(n, ptr, row, val, soln, resid);
   for(int i=0; i<n; i++) resid[i] -= rhs[i];
   /* Calculate scaled norms */
   T resid_inf = max_abs_value(resid, resid+n);
   T A_inf = matrix_inf_norm(n, ptr, row, val);
   T x_inf = max_abs_value(soln, soln+n);
   T b_inf = max_abs_value(rhs, rhs+n);
   delete[] resid;
   return resid_inf / (A_inf*x_inf + b_inf);
}

int main(int argc, char *const * argv) {
   /* Process options */
   DriverOptions options(argc, argv);

   /* Read in a matrix */
   printf("Reading...");
   struct spral_rb_options rb_options;
   spral_rb_default_options(&rb_options);
   if(options.force_posdef) rb_options.values = -3;
   void *read_handle;
   int m, n, *ptr, *row, *col, flag;
   double *val;
   flag = spral_rb_read_i32d("matrix.rb", &read_handle, &m, &n, &ptr, &row, 
         &col, &val, &rb_options, nullptr, nullptr, nullptr);
   if(flag) {
      printf("\nRutherford-Boeing read failed with error %d\n", flag);
      return 1;
   }
   printf("ok\n");

   if(options.print_matrix) {
      printf("Matrix:\n");
      for(int i=0; i<n; i++) {
         printf("%4d [%4d:%4d]:", i, ptr[i], ptr[i+1]-1);
         for(int j=ptr[i]; j<ptr[i+1]; ++j) printf(" %10d", row[j]);
         printf("\n");
         printf("%4d [%4d:%4d]:", i, ptr[i], ptr[i+1]-1);
         for(int j=ptr[i]; j<ptr[i+1]; ++j) printf(" %10.2e", val[j]);
         printf("\n");
      }
   }

   /* Analyse */
   printf("\nAnalyse...");
   struct timespec t1, t2;
   clock_gettime(CLOCK_REALTIME, &t1);
   spral::ics::SymbolicFactor sfact(n, ptr, row, options.nemin);
   clock_gettime(CLOCK_REALTIME, &t2);
   printf("ok\n");
   printf("Analyse took %e\n", tdiff(t1, t2));
   printf("Predicted nfact = %.2le\n", (double) sfact.get_nfact());
   printf("Predicted nflop = %.2le\n", (double) sfact.get_nflop());

   if(options.print_matrix) {
      printf("perm = ");
      int const* perm = sfact.get_perm();
      for(int i=0; i<n; i++) printf(" %d", perm[i]);
      printf("\n");
   }

   /* Factorize */
   printf("\nFactorize...");
   clock_gettime(CLOCK_REALTIME, &t1);
   spral::ics::NumericFactor nfact(sfact, val);
   clock_gettime(CLOCK_REALTIME, &t2);
   printf("ok\n");
   printf("Factorize took %e\n", tdiff(t1, t2));
   if(options.print_factors) nfact.print();

   /* Generate solution and rhs */
   double *x = new double[n];
   for(int i=0; i<n; ++i) x[i] = 0.1*(i+1);
   double *rhs = new double[n];
   spmv(n, ptr, row, val, x, rhs);
   if(options.print_matrix) {
      printf("x = ");
      for(int i=0; i<n; ++i) printf(" %e", x[i]);
      printf("\n");
      printf("rhs = ");
      for(int i=0; i<n; ++i) printf(" %e", rhs[i]);
      printf("\n");
   }

   /* Solve */
   printf("\nSolve...");
   double *soln = new double[n];
   for(int i=0; i<n; i++) soln[i] = rhs[i];
   clock_gettime(CLOCK_REALTIME, &t1);
   nfact.solve(soln);
   clock_gettime(CLOCK_REALTIME, &t2);
   printf("ok\n");
   printf("solve took %e\n", tdiff(t1, t2));
   if(options.print_matrix) {
      printf("soln = ");
      for(int i=0; i<n; ++i) printf(" %e", soln[i]);
      printf("\n");
   }

   /* Errors */
   double fwd_err = 0.0;
   for(int i=0; i<n; i++) fwd_err = std::max(fwd_err, fabs(soln[i]-x[i]));
   printf("\nfwd error = %e\n", fwd_err);
   printf("bwd error = %e\n", calc_bwd_error(n, ptr, row, val, soln, rhs));

   /* Cleanup */
   delete[] soln;
   delete[] rhs;
   delete[] x;
   spral_rb_free_handle(&read_handle);

   return 0; // Success
}
