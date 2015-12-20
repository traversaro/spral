#ifndef SPRAL_METIS_WRAPPER_H
#define SPRAL_METIS_WRAPPER_H

#ifdef __cplusplus
extern "C" {
#endif

/* Order a matrix using METIS */
int spral_metis_order(int n, const int ptr[], const int row[], int perm[], int invp[], int base);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* SPRAL_METIS_WRAPPER_H */
