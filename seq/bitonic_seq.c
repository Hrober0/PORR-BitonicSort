#define _POSIX_C_SOURCE 200809L

#include "bitonic_common.h"

#include <errno.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define VERSION_NAME "seq"
#define DEFAULT_CSV  "results/seq.csv"
#define DEFAULT_RUN  1
#define LOG2_N_MAX   31

/* Instrukcja wywołania. */
static void usage(const char *argv0)
{
    fprintf(stderr,
            "usage: %s <log2_n> <p> <seed> [results.csv] [run_id]\n"
            "  log2_n       n = 2^log2_n, integer in 0..%d\n"
            "  p            worker count (sequential version records p=1)\n"
            "  seed         LCG seed (same seed => same input in every version)\n"
            "  results.csv  output CSV (default: %s)\n"
            "  run_id       run number written to the CSV (default: %d)\n",
            argv0, LOG2_N_MAX, DEFAULT_CSV, DEFAULT_RUN);
}

/* Parsowanie liczby całkowitej. */
static int parse_i32(const char *text, int *out)
{
    char *end;
    long value;

    errno = 0;
    value = strtol(text, &end, 10);
    if (errno != 0 || end == text || *end != '\0' ||
        value < (long)INT_MIN || value > (long)INT_MAX) {
        return -1;
    }
    *out = (int)value;
    return 0;
}

/* Pomiar sekwencyjnego sortowania bitonicznego. */
int main(int argc, char **argv)
{
    int log2_n;
    int p_arg;
    int seed;
    int run_id;
    int p_recorded;
    int ok;
    int csv_rc;
    long n;
    const char *csv_path;
    double *a;
    double *original;
    uint64_t sum_bits;
    uint64_t xor_bits;
    double t0;
    double t1;
    double time_s;
    double cauchy_norm;
    double residual;

    if (argc < 4 || argc > 6) {
        usage(argv[0]);
        return 2;
    }
    if (parse_i32(argv[1], &log2_n) != 0 || log2_n < 0 || log2_n > LOG2_N_MAX) {
        fprintf(stderr, "error: log2_n must be an integer in 0..%d\n", LOG2_N_MAX);
        usage(argv[0]);
        return 2;
    }
    if (parse_i32(argv[2], &p_arg) != 0) {
        fprintf(stderr, "error: p must be an integer\n");
        usage(argv[0]);
        return 2;
    }
    if (parse_i32(argv[3], &seed) != 0) {
        fprintf(stderr, "error: seed must be an integer\n");
        usage(argv[0]);
        return 2;
    }

    csv_path = DEFAULT_CSV;
    run_id = DEFAULT_RUN;
    if (argc >= 5) {
        csv_path = argv[4];
        if (csv_path[0] == '\0') {
            fprintf(stderr, "error: results.csv path is empty\n");
            return 2;
        }
    }
    if (argc >= 6 && parse_i32(argv[5], &run_id) != 0) {
        fprintf(stderr, "error: run_id must be an integer\n");
        usage(argv[0]);
        return 2;
    }

    p_recorded = 1;
    if (p_arg != 1) {
        fprintf(stderr,
                "warning: sequential version ignores p=%d and records p=1\n",
                p_arg);
    }

    n = 1L << log2_n;
    if ((unsigned long)n > ULONG_MAX / sizeof(double)) {
        fprintf(stderr, "error: n=%ld does not fit in memory\n", n);
        return 1;
    }

    a = (double *)malloc((size_t)n * sizeof(double));
    if (a == NULL) {
        fprintf(stderr, "error: failed to allocate %ld doubles\n", n);
        return 1;
    }

    gen_data(a, n, seed);

    /* Kopia wejścia albo suma kontrolna, zależnie od n. */
    original = NULL;
    sum_bits = 0;
    xor_bits = 0;
    if (n <= BITONIC_QSORT_MAX_N) {
        original = (double *)malloc((size_t)n * sizeof(double));
        if (original == NULL) {
            fprintf(stderr, "error: failed to allocate the reference copy\n");
            free(a);
            return 1;
        }
        memcpy(original, a, (size_t)n * sizeof(double));
    } else {
        array_checksum(a, n, &sum_bits, &xor_bits);
    }

    fprintf(stderr, "bitonic_seq: sorting n=%ld (p=%d, seed=%d, run=%d)\n",
            n, p_recorded, seed, run_id);
    fflush(stderr);

    t0 = timer_now();
    bitonic_sort_seq(a, n);
    t1 = timer_now();
    time_s = t1 - t0;

    cauchy_norm = check_cauchy(a, n);
    if (original != NULL) {
        residual = check_residual(a, original, n);
        free(original);
    } else {
        residual = check_residual_checksum(a, n, sum_bits, xor_bits);
    }
    free(a);

    ok = (cauchy_norm == 0.0 && residual == 0.0);

    printf("n=%ld\n", n);
    printf("time_s=%.6f\n", time_s);
    printf("cauchy_norm=%.6e\n", cauchy_norm);
    printf("residual=%.6e\n", residual);
    printf("status=%s\n", ok ? "OK" : "FAIL");

    csv_rc = append_result_csv(csv_path, VERSION_NAME, n, p_recorded, run_id,
                               seed, time_s, cauchy_norm, residual, ok);
    if (csv_rc != 0) {
        return 1;
    }
    return ok ? 0 : 1;
}
