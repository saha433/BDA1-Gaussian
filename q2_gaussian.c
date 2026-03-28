#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include <stdint.h>
#include <time.h>

#define NUM_THREADS 8
#define VECTOR_LEN 10000000000ULL /* 10^10 */

/* Gaussian: mean=0, std=0.8 -> most values in (-1,1), round to {-1,0,1} */
#define GAUSS_MEAN 0.0
#define GAUSS_STD 0.8

typedef struct
{
    int thread_id;
    uint64_t count;
    int64_t local_dot; /* partial dot product */
    /* For a true 3-D cross product we only need the first 3 elements */
    int a0, a1, a2; /* first 3 elements of A from this thread (thread 0 only) */
    int b0, b1, b2; /* first 3 elements of B from this thread (thread 0 only) */
} ThreadArg;

/* Box-Muller */
static inline double gaussian_bm(unsigned int *seed)
{
    double u1, u2;
    do
    {
        u1 = (double)rand_r(seed) / RAND_MAX;
    } while (u1 == 0.0);
    u2 = (double)rand_r(seed) / RAND_MAX;
    return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}

/* Map Gaussian sample -> {-1, 0, 1} by rounding to nearest integer */
static inline int gauss_to_discrete(unsigned int *seed)
{
    double v = GAUSS_MEAN + GAUSS_STD * gaussian_bm(seed);
    if (v <= -0.5)
        return -1;
    if (v >= 0.5)
        return 1;
    return 0;
}

void *thread_func(void *arg)
{
    ThreadArg *t = (ThreadArg *)arg;
    unsigned int seedA = (unsigned int)(t->thread_id * 111111111UL + 1);
    unsigned int seedB = (unsigned int)(t->thread_id * 999999999UL + 2);

    int64_t dot = 0;

    for (uint64_t i = 0; i < t->count; i++)
    {
        int a = gauss_to_discrete(&seedA);
        int b = gauss_to_discrete(&seedB);

        /* save first 3 elements for 3-D cross product (thread 0) */
        if (t->thread_id == 0 && i < 3)
        {
            if (i == 0)
            {
                t->a0 = a;
                t->b0 = b;
            }
            if (i == 1)
            {
                t->a1 = a;
                t->b1 = b;
            }
            if (i == 2)
            {
                t->a2 = a;
                t->b2 = b;
            }
        }

        dot += (int64_t)a * b;
    }

    t->local_dot = dot;
    return NULL;
}

int main(void)
{
    pthread_t threads[NUM_THREADS];
    ThreadArg args[NUM_THREADS];

    uint64_t chunk = VECTOR_LEN / NUM_THREADS;

    printf("BDA Assignment-1 Q2: Dot & Cross Products (Gaussian)\n");
    printf("Vector length  : 10^10 = %llu\n", (unsigned long long)VECTOR_LEN);
    printf("Distribution   : Gaussian(mean=%.1f, std=%.1f) rounded to {-1,0,1}\n",
           GAUSS_MEAN, GAUSS_STD);
    printf("Threads        : %d\n\n", NUM_THREADS);

    clock_t t0 = clock();

    for (int i = 0; i < NUM_THREADS; i++)
    {
        args[i].thread_id = i;
        args[i].count = (i == NUM_THREADS - 1)
                            ? (VECTOR_LEN - (uint64_t)i * chunk)
                            : chunk;
        pthread_create(&threads[i], NULL, thread_func, &args[i]);
    }

    int64_t gdot = 0;
    for (int i = 0; i < NUM_THREADS; i++)
    {
        pthread_join(threads[i], NULL);
        gdot += args[i].local_dot;
    }

    /* 3-D Cross product using first 3 elements from thread 0 */
    int a0 = args[0].a0, a1 = args[0].a1, a2 = args[0].a2;
    int b0 = args[0].b0, b1 = args[0].b1, b2 = args[0].b2;
    int cx = a1 * b2 - a2 * b1;
    int cy = a2 * b0 - a0 * b2;
    int cz = a0 * b1 - a1 * b0;

    clock_t t1 = clock();
    double elapsed = (double)(t1 - t0) / CLOCKS_PER_SEC;

    printf("\n===== RESULTS =====\n");
    printf("Dot Product (sum over 10^10 elements) : %lld\n", (long long)gdot);
    printf("\n3-D Cross Product (using first 3 elements of each vector):\n");
    printf("  A = (%d, %d, %d)\n", a0, a1, a2);
    printf("  B = (%d, %d, %d)\n", b0, b1, b2);
    printf("  A x B = (%d, %d, %d)\n", cx, cy, cz);

    return 0;
}
