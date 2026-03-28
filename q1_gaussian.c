#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include <stdint.h>
#include <float.h>
#include <time.h>

#define NUM_THREADS 8
#define TOTAL_ELEMENTS (1ULL << 40) /* 2^40 */
#define MAX_VAL 1000000000ULL       /* 10^9 */

/* Gaussian parameters: mean=5e8, std=1.5e8 so values spread over [0,10^9] */
#define GAUSS_MEAN 500000000.0
#define GAUSS_STD 150000000.0

/* Per-thread result */
typedef struct
{
    int thread_id;
    uint64_t start;
    uint64_t count;
    uint64_t local_min;
    uint64_t local_max;
    __uint128_t local_sum; /* 128-bit to avoid overflow: 2^40 * 10^9 needs ~70 bits */
} ThreadArg;

/* Box-Muller transform: generate one Gaussian sample */
static inline double gaussian(unsigned int *seed)
{
    double u1, u2;
    do
    {
        u1 = (double)rand_r(seed) / RAND_MAX;
    } while (u1 == 0.0);
    u2 = (double)rand_r(seed) / RAND_MAX;
    return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}

/* Clamp Gaussian sample to [0, MAX_VAL] */
static inline uint64_t gaussian_sample(unsigned int *seed)
{
    double val = GAUSS_MEAN + GAUSS_STD * gaussian(seed);
    if (val < 0.0)
        return 0;
    if (val > (double)MAX_VAL)
        return MAX_VAL;
    return (uint64_t)val;
}

void *thread_func(void *arg)
{
    ThreadArg *t = (ThreadArg *)arg;
    unsigned int seed = (unsigned int)(t->thread_id * 123456789UL + 987654321UL);

    uint64_t lmin = UINT64_MAX;
    uint64_t lmax = 0;
    __uint128_t lsum = 0;

    for (uint64_t i = 0; i < t->count; i++)
    {
        uint64_t v = gaussian_sample(&seed);
        if (v < lmin)
            lmin = v;
        if (v > lmax)
            lmax = v;
        lsum += v;
    }

    t->local_min = lmin;
    t->local_max = lmax;
    t->local_sum = lsum;
    return NULL;
}

int main(void)
{
    pthread_t threads[NUM_THREADS];
    ThreadArg args[NUM_THREADS];

    uint64_t chunk = TOTAL_ELEMENTS / NUM_THREADS;

    printf("BDA Assignment-1 Q1: Min, Max, Mean (Gaussian)\n");
    printf("Total elements : 2^40 = %llu\n", (unsigned long long)TOTAL_ELEMENTS);
    printf("Distribution   : Gaussian(mean=%.0f, std=%.0f) clamped to [0, 10^9]\n",
           GAUSS_MEAN, GAUSS_STD);
    printf("Threads        : %d\n\n", NUM_THREADS);

    clock_t t0 = clock();

    for (int i = 0; i < NUM_THREADS; i++)
    {
        args[i].thread_id = i;
        args[i].start = (uint64_t)i * chunk;
        args[i].count = (i == NUM_THREADS - 1)
                            ? (TOTAL_ELEMENTS - (uint64_t)i * chunk)
                            : chunk;
        pthread_create(&threads[i], NULL, thread_func, &args[i]);
    }

    uint64_t gmin = UINT64_MAX, gmax = 0;
    __uint128_t gsum = 0;

    for (int i = 0; i < NUM_THREADS; i++)
    {
        pthread_join(threads[i], NULL);
        if (args[i].local_min < gmin)
            gmin = args[i].local_min;
        if (args[i].local_max > gmax)
            gmax = args[i].local_max;
        gsum += args[i].local_sum;
    }

    double mean = (double)(gsum) / (double)TOTAL_ELEMENTS;

    clock_t t1 = clock();
    double elapsed = (double)(t1 - t0) / CLOCKS_PER_SEC;

    printf("\n===== RESULTS =====\n");
    printf("Minimum : %llu\n", (unsigned long long)gmin);
    printf("Maximum : %llu\n", (unsigned long long)gmax);
    printf("Mean    : %.6f\n", mean);

    return 0;
}
