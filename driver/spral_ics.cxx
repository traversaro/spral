#include <boost/program_options.hpp>
#include <iostream>
#include <time.h>

#include "spral.h"
#include "src/ics/SymbolicFactor.hxx"
#include "src/ics/NumericFactor.hxx"

float tdiff(const struct timespec &t1, const struct timespec &t2) {
   return t2.tv_sec - t1.tv_sec + 1e-9*(t2.tv_nsec - t1.tv_nsec);
}

void process_options(int argc, char *const * argv, bool& print_matrix);

int main(int argc, char *const * argv) {
   /* Process options */
   bool print_matrix;
   process_options(argc, argv, print_matrix);

   /* Read in a matrix */
   printf("Reading...");
   struct spral_rb_options rb_options;
   spral_rb_default_options(&rb_options);
   rb_options.values = 2; /* make up values if necessary */
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

   if(print_matrix) {
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
   spral::ics::SymbolicFactor sfact(n, ptr, row, 8);
   clock_gettime(CLOCK_REALTIME, &t2);
   printf("ok\n");
   printf("Analyse took %e\n", tdiff(t1, t2));
   printf("Predicted nfact = %.2le\n", (double) sfact.get_nfact());
   printf("Predicted nflop = %.2le\n", (double) sfact.get_nflop());

   /* Factorize */
   printf("\nFactorize...");
   clock_gettime(CLOCK_REALTIME, &t1);
   spral::ics::NumericFactor nfact(sfact, val);
   clock_gettime(CLOCK_REALTIME, &t2);
   printf("ok\n");
   printf("Factorize took %e\n", tdiff(t1, t2));

   /* Solve */

   /* Cleanup */
   spral_rb_free_handle(&read_handle);

   return 0; // Success
}

void process_options(int argc, char *const * argv, bool& print_matrix) {
   boost::program_options::options_description desc("Allowed options");
   desc.add_options()
      ("help", "produce help message")
      ("print-matrix", "print input matrix")
      ;
   boost::program_options::variables_map vm;
   boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);
   boost::program_options::notify(vm);

   if( vm.count("help") ) {
      std::cout << desc << std::endl;
      exit(1);
   }

   print_matrix = ( vm.count("print-matrix") );
}
