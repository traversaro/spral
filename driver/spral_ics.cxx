#include <time.h>

#include "spral.h"
#include "src/ics/symbolic.hxx"

float tdiff(const struct timespec &t1, const struct timespec &t2) {
   return t2.tv_sec - t1.tv_sec + 1e-9*(t2.tv_nsec - t1.tv_nsec);
}

int main(void) {
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

   /* Analyse */
   printf("Analyse...");
   struct timespec t1, t2;
   clock_gettime(CLOCK_REALTIME, &t1);
   spral::ics::symbolic afact(n, ptr, row, 8);
   clock_gettime(CLOCK_REALTIME, &t2);
   printf("ok\n");
   printf("Analyse took %e\n", tdiff(t1, t2));
   printf("Predicted nfact = %ld\n", afact.nfact);
   printf("Predicted nflop = %ld\n", afact.nflop);

   /* Factorize */
   /* Solve */

   /* Cleanup */
   spral_rb_free_handle(&read_handle);

   return 0; // Success
}
