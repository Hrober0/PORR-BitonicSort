#ifndef BITONIC_COMMON_H
#define BITONIC_COMMON_H

/* Wspólny interfejs sortowania bitonicznego. */

#include <stdint.h>

/* Próg n dla porównania z qsort. */
#define BITONIC_QSORT_MAX_N (1L << 25)

/* Deterministyczne generowanie tablicy. */
void gen_data(double *a, long n, int seed);

/* Zamiana pary elementów w zadanym kierunku. */
static inline void compare_exchange(double *a, long i, long j, int dir)
{
    dir = (dir != 0);
    if ((a[i] > a[j]) == dir) {
        double tmp = a[i];
        a[i] = a[j];
        a[j] = tmp;
    }
}

/* Sekwencyjne sortowanie bitoniczne. */
void bitonic_sort_seq(double *a, long n);

/* Scalanie bitoniczne w kierunku dir. */
void bitonic_merge_seq(double *a, long n, int dir);

/* Czas monotoniczny w sekundach. */
double timer_now(void);

/* Norma różnicy po dodatkowym scaleniu. */
double check_cauchy(const double *sorted, long n);

/* Odległość wyniku od referencji. */
double check_residual(const double *sorted, const double *original, long n);

/* Suma i XOR bitów tablicy. */
void array_checksum(const double *a, long n, uint64_t *sum_bits, uint64_t *xor_bits);

/* Weryfikacja dużego n przez sumę kontrolną. */
double check_residual_checksum(const double *sorted, long n,
                               uint64_t sum_bits, uint64_t xor_bits);

/* Dopisanie wyniku do CSV. */
int append_result_csv(const char *path, const char *version, long n, int p,
                      int run, int seed, double time_s, double cauchy_norm,
                      double residual, int ok);

#endif /* BITONIC_COMMON_H */
