#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <string.h>
#include <math.h>

#define TOTAL_ELEMENTS 1000000000ULL
#define NUM_THREADS 8

#define RANGE_MIN (-1073741824LL)
#define RANGE_MAX (1073741824LL)
#define RANGE_SPAN 2147483649ULL

#define NUM_BUCKETS 1024

// Gaussian parameters
#define MEAN 0.0
#define STD 300000000.0 // spread (adjustable)

typedef struct
{
    unsigned long long start;
    unsigned long long end;
    unsigned long long bucket_counts[NUM_BUCKETS];
} ThreadArgs;

// Gaussian generator
double gaussian()
{
    double u1 = (double)rand() / RAND_MAX;
    double u2 = (double)rand() / RAND_MAX;

    if (u1 < 1e-9)
        u1 = 1e-9;

    return sqrt(-2.0 * log(u1)) * cos(2 * M_PI * u2);
}

// bucket mapping
static inline int get_bucket(long long val)
{
    unsigned long long shifted = (unsigned long long)(val - RANGE_MIN);
    return (int)(shifted * NUM_BUCKETS / RANGE_SPAN);
}

void *histogram_worker(void *arg)
{
    ThreadArgs *a = (ThreadArgs *)arg;
    memset(a->bucket_counts, 0, sizeof(a->bucket_counts));

    srand(a->start + 123);

    unsigned long long count = a->end - a->start;

    for (unsigned long long i = 0; i < count; i++)
    {

        double z = gaussian();
        double val_d = MEAN + STD * z;

        // clamp
        if (val_d < RANGE_MIN)
            val_d = RANGE_MIN;
        if (val_d > RANGE_MAX)
            val_d = RANGE_MAX;

        long long val = (long long)val_d;

        int b = get_bucket(val);
        a->bucket_counts[b]++;
    }

    return NULL;
}

int main()
{
    printf("BDA Assignment-1 Q4: Median Estimation via Histogram (Gaussian)\n");
    printf("Total elements : 10^9 = %llu\n", TOTAL_ELEMENTS);
    printf("Range          : -2^30 to 2^30\n");
    printf("Distribution   : Gaussian\n");
    printf("Threads        : %d\n\n", NUM_THREADS);

    pthread_t threads[NUM_THREADS];
    ThreadArgs args[NUM_THREADS];

    unsigned long long chunk = TOTAL_ELEMENTS / NUM_THREADS;

    for (int t = 0; t < NUM_THREADS; t++)
    {
        args[t].start = (unsigned long long)t * chunk;
        args[t].end = (t == NUM_THREADS - 1) ? TOTAL_ELEMENTS : (unsigned long long)(t + 1) * chunk;

        pthread_create(&threads[t], NULL, histogram_worker, &args[t]);
    }

    for (int t = 0; t < NUM_THREADS; t++)
        pthread_join(threads[t], NULL);

    // merge histograms
    unsigned long long global_hist[NUM_BUCKETS];
    memset(global_hist, 0, sizeof(global_hist));

    for (int t = 0; t < NUM_THREADS; t++)
        for (int b = 0; b < NUM_BUCKETS; b++)
            global_hist[b] += args[t].bucket_counts[b];

    // find median bucket
    unsigned long long median_rank = TOTAL_ELEMENTS / 2;
    unsigned long long cumsum = 0;
    int median_bucket = 0;

    for (int b = 0; b < NUM_BUCKETS; b++)
    {
        cumsum += global_hist[b];
        if (cumsum >= median_rank)
        {
            median_bucket = b;
            break;
        }
    }

    long long bucket_lo = RANGE_MIN +
                          (long long)((unsigned long long)median_bucket * RANGE_SPAN / NUM_BUCKETS);

    long long bucket_hi = RANGE_MIN +
                          (long long)((unsigned long long)(median_bucket + 1) * RANGE_SPAN / NUM_BUCKETS);

    long long median_estimate = bucket_lo + (bucket_hi - bucket_lo) / 2;

    printf("Median bucket    : %d\n", median_bucket);
    printf("Bucket range     : [%lld, %lld]\n", bucket_lo, bucket_hi);
    printf("Median estimate  : %lld\n", median_estimate);

    // local medians
    printf("\nLocal median estimates:\n");

    long long local_medians[NUM_THREADS];

    for (int t = 0; t < NUM_THREADS; t++)
    {
        unsigned long long half = (args[t].end - args[t].start) / 2;

        unsigned long long cum = 0;
        int mb = 0;

        for (int b = 0; b < NUM_BUCKETS; b++)
        {
            cum += args[t].bucket_counts[b];
            if (cum >= half)
            {
                mb = b;
                break;
            }
        }

        long long lo = RANGE_MIN +
                       (long long)((unsigned long long)mb * RANGE_SPAN / NUM_BUCKETS);

        long long hi = RANGE_MIN +
                       (long long)((unsigned long long)(mb + 1) * RANGE_SPAN / NUM_BUCKETS);

        local_medians[t] = lo + (hi - lo) / 2;

        printf("  Thread %d -> %lld\n", t, local_medians[t]);
    }

    // median of medians
    for (int i = 0; i < NUM_THREADS - 1; i++)
        for (int j = i + 1; j < NUM_THREADS; j++)
            if (local_medians[i] > local_medians[j])
            {
                long long tmp = local_medians[i];
                local_medians[i] = local_medians[j];
                local_medians[j] = tmp;
            }

    long long mom = local_medians[NUM_THREADS / 2];

    printf("\nMedian of Medians : %lld\n", mom);
    printf("Histogram Median  : %lld\n", median_estimate);

    return 0;
}