#ifndef SPRAL_METIS_WRAPPER_H
#define SPRAL_METIS_WRAPPER_H

/* Order a matrix using METIS */
int spral_metis_order(int n, const int ptr[], const int row[], int perm[], int invp[], int base);

#endif /* SPRAL_METIS_WRAPPER_H */
