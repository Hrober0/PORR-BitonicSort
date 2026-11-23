#define _POSIX_C_SOURCE 200809L

#include "bitonic_common.h"

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* Parametry generatora LCG. */
#define LCG_MULTIPLIER 6364136223846793005ULL
#define LCG_INCREMENT  1442695040888963407ULL

/* Deterministyczne generowanie tablicy. */
void gen_data(double *a, long n, int seed)
{
    uint64_t state;
    long i;

    if (a == NULL || n <= 0) {
        return;
    }

    state = (uint64_t)(uint32_t)seed;
    for (i = 0; i < n; i++) {
        state = state * LCG_MULTIPLIER + LCG_INCREMENT;
        a[i] = (double)(state >> 11) * 0x1p-53;
    }
}

/* Alokacja tablicy double. */
static int alloc_doubles(long n, double **out)
{
    if (n <= 0 || (unsigned long)n > ULONG_MAX / sizeof(double)) {
        return -1;
    }
    *out = (double *)malloc((size_t)n * sizeof(double));
    return (*out == NULL) ? -1 : 0;
}

/* Sekwencyjne sortowanie bitoniczne. */
void bitonic_sort_seq(double *a, long n)
{
    long k;
    long j;
    long i;

    if (a == NULL || n <= 1) {
        return;
    }

    for (k = 2; k <= n; k <<= 1) {
        for (j = k >> 1; j > 0; j >>= 1) {
            for (i = 0; i < n; i++) {
                long l = i ^ j;
                if (l > i) {
                    compare_exchange(a, i, l, (i & k) == 0);
                }
            }
        }
    }
}

/* Scalanie bitoniczne w kierunku dir. */
void bitonic_merge_seq(double *a, long n, int dir)
{
    long j;
    long i;

    if (a == NULL || n <= 1) {
        return;
    }

    dir = (dir != 0);
    for (j = n >> 1; j > 0; j >>= 1) {
        for (i = 0; i < n; i++) {
            long l = i ^ j;
            if (l > i) {
                compare_exchange(a, i, l, dir);
            }
        }
    }
}

/* Czas monotoniczny w sekundach. */
double timer_now(void)
{
    struct timespec ts;

    if (clock_gettime(CLOCK_MONOTONIC, &ts) != 0) {
        return 0.0;
    }
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

/* Porównanie wartości double. */
static int same_double(double x, double y)
{
    uint64_t bx;
    uint64_t by;

    if (x == y) {
        return 1;
    }
    bx = 0;
    by = 0;
    memcpy(&bx, &x, sizeof(bx));
    memcpy(&by, &y, sizeof(by));
    return bx == by;
}

/* Odległość euklidesowa tablic. */
static double euclidean_diff(const double *x, const double *y, long n)
{
    double sum_sq;
    long i;

    sum_sq = 0.0;
    for (i = 0; i < n; i++) {
        double d;

        if (same_double(x[i], y[i])) {
            continue;
        }
        d = x[i] - y[i];
        if (!isfinite(d)) {
            return INFINITY;
        }
        sum_sq += d * d;
    }
    if (!isfinite(sum_sq)) {
        return INFINITY;
    }
    return sqrt(sum_sq);
}

/* Liczba sąsiednich inwersji. */
static long count_adjacent_inversions(const double *a, long n)
{
    long count;
    long i;

    count = 0;
    for (i = 0; i + 1 < n; i++) {
        if (!(a[i] <= a[i + 1])) {
            count++;
        }
    }
    return count;
}

/* Komparator dla qsort. */
static int compare_doubles(const void *lhs, const void *rhs)
{
    double x;
    double y;

    x = *(const double *)lhs;
    y = *(const double *)rhs;
    if (x < y) {
        return -1;
    }
    if (x > y) {
        return 1;
    }
    return 0;
}

/* Suma i XOR bitów tablicy. */
void array_checksum(const double *a, long n, uint64_t *sum_bits, uint64_t *xor_bits)
{
    uint64_t sum;
    uint64_t xorb;
    long i;

    sum = 0;
    xorb = 0;
    if (a != NULL && n > 0) {
        for (i = 0; i < n; i++) {
            uint64_t bits;

            bits = 0;
            memcpy(&bits, &a[i], sizeof(bits));
            sum += bits;
            xorb ^= bits;
        }
    }
    if (sum_bits != NULL) {
        *sum_bits = sum;
    }
    if (xor_bits != NULL) {
        *xor_bits = xorb;
    }
}

/* Norma różnicy po dodatkowym scaleniu. */
double check_cauchy(const double *sorted, long n)
{
    double *work;
    double norm;

    if (sorted == NULL || n < 0) {
        return INFINITY;
    }
    if (n <= 1) {
        return 0.0;
    }

    if (alloc_doubles(n, &work) != 0) {
        fprintf(stderr, "check_cauchy: out of memory for n=%ld\n", n);
        return INFINITY;
    }

    memcpy(work, sorted, (size_t)n * sizeof(double));
    bitonic_merge_seq(work, n, 1);
    norm = euclidean_diff(work, sorted, n);
    free(work);
    return norm;
}

/* Weryfikacja dużego n przez sumę kontrolną. */
double check_residual_checksum(const double *sorted, long n,
                               uint64_t sum_bits, uint64_t xor_bits)
{
    uint64_t got_sum;
    uint64_t got_xor;
    long inversions;

    if (sorted == NULL || n < 0) {
        return INFINITY;
    }
    if (n == 0) {
        return 0.0;
    }

    inversions = count_adjacent_inversions(sorted, n);
    array_checksum(sorted, n, &got_sum, &got_xor);

    if (inversions == 0 && got_sum == sum_bits && got_xor == xor_bits) {
        return 0.0;
    }
    if (inversions > 0) {
        return (double)inversions;
    }
    return 1.0;
}

/* Odległość wyniku od referencji. */
double check_residual(const double *sorted, const double *original, long n)
{
    double *ref;
    double norm;
    long inversions;
    uint64_t sum_bits;
    uint64_t xor_bits;

    if (sorted == NULL || original == NULL || n < 0) {
        return INFINITY;
    }
    if (n == 0) {
        return 0.0;
    }

    if (n > BITONIC_QSORT_MAX_N) {
        array_checksum(original, n, &sum_bits, &xor_bits);
        return check_residual_checksum(sorted, n, sum_bits, xor_bits);
    }

    if (alloc_doubles(n, &ref) != 0) {
        fprintf(stderr, "check_residual: out of memory for n=%ld\n", n);
        return INFINITY;
    }
    memcpy(ref, original, (size_t)n * sizeof(double));
    qsort(ref, (size_t)n, sizeof(double), compare_doubles);

    norm = euclidean_diff(sorted, ref, n);
    inversions = count_adjacent_inversions(sorted, n);
    free(ref);

    if (inversions > 0 && norm == 0.0) {
        return 1.0;
    }
    return norm;
}

/* Sprawdzenie, czy plik jest pusty. */
static int file_missing_or_empty(const char *path)
{
    FILE *fp;
    int c;

    fp = fopen(path, "r");
    if (fp == NULL) {
        return 1;
    }
    c = fgetc(fp);
    fclose(fp);
    return c == EOF;
}

/* Znacznik czasu UTC. */
static void format_utc_timestamp(char *buf, size_t buflen)
{
    time_t now;
    struct tm tm_utc;

    now = time(NULL);
    if (gmtime_r(&now, &tm_utc) == NULL ||
        strftime(buf, buflen, "%Y-%m-%dT%H:%M:%S", &tm_utc) == 0) {
        snprintf(buf, buflen, "1970-01-01T00:00:00");
    }
}

/* Dopisanie wyniku do CSV. */
int append_result_csv(const char *path, const char *version, long n, int p,
                      int run, int seed, double time_s, double cauchy_norm,
                      double residual, int ok)
{
    int write_header;
    FILE *fp;
    char timestamp[32];
    int rc;

    if (path == NULL || path[0] == '\0' || version == NULL) {
        fprintf(stderr, "append_result_csv: missing path or version\n");
        return -1;
    }

    write_header = file_missing_or_empty(path);
    fp = fopen(path, "a");
    if (fp == NULL) {
        fprintf(stderr, "append_result_csv: cannot open %s\n", path);
        return -1;
    }

    rc = 0;
    if (write_header) {
        if (fprintf(fp,
                    "version,n,p,run,seed,time_s,cauchy_norm,residual,ok,timestamp\n") < 0) {
            rc = -1;
        }
    }

    format_utc_timestamp(timestamp, sizeof(timestamp));
    if (rc == 0 &&
        fprintf(fp, "%s,%ld,%d,%d,%d,%.6f,%.6e,%.6e,%d,%s\n",
                version, n, p, run, seed, time_s, cauchy_norm, residual, ok,
                timestamp) < 0) {
        rc = -1;
    }

    if (fclose(fp) == EOF) {
        rc = -1;
    }
    if (rc != 0) {
        fprintf(stderr, "append_result_csv: write failed for %s\n", path);
    }
    return rc;
}
