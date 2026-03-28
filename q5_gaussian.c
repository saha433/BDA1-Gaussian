#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef _WIN32
#include <windows.h>
typedef HANDLE thread_t;
#define THREAD_FUNC DWORD WINAPI
static void thread_create(thread_t *h, LPTHREAD_START_ROUTINE fn, void *arg)
{
    *h = CreateThread(NULL, 0, fn, arg, 0, NULL);
}
static void thread_join(thread_t h)
{
    WaitForSingleObject(h, INFINITE);
    CloseHandle(h);
}
/* Atomic add — Windows intrinsic */
#include <intrin.h>
#define ATOMIC_ADD(ptr, val) _InterlockedExchangeAdd((volatile LONG *)(ptr), (val))
#else
#include <pthread.h>
typedef pthread_t thread_t;
#define THREAD_FUNC void *
static void thread_create(thread_t *h, void *(*fn)(void *), void *arg)
{
    pthread_create(h, NULL, fn, arg);
}
static void thread_join(thread_t h) { pthread_join(h, NULL); }
#define ATOMIC_ADD(ptr, val) __sync_fetch_and_add((ptr), (val))
#endif

#define THREADS 4
#define BINS 2000 /* finer resolution for Gaussian tails */
#define GAUSS_MEAN 0.5
#define GAUSS_SIGMA 0.15
#define FIXED_SEED 42U
#define TOTAL_SECONDS 3600     /* 1 hour */
#define VALUES_PER_SEC 100000L /* 100 000 values / second             */

#define TOTAL_MINUTES (TOTAL_SECONDS / 60)   /* 60 */
#define TEN_MIN_BLOCKS (TOTAL_SECONDS / 600) /* 6  */
#define VALUES_PER_MIN (VALUES_PER_SEC * 60L)
#define VALUES_PER_10M (VALUES_PER_SEC * 600L)
#define TOTAL_VALUES ((long)VALUES_PER_SEC * TOTAL_SECONDS)

static long global_hist[BINS];
static long *minute_hist;  /* [TOTAL_MINUTES][BINS] flat */
static long *ten_min_hist; /* [TEN_MIN_BLOCKS][BINS] flat */

#define MIN_HIST(m, b) minute_hist[(m) * BINS + (b)]
#define TEN_HIST(t, b) ten_min_hist[(t) * BINS + (b)]

static inline unsigned int lcg_next(unsigned int *s)
{
    *s = *s * 1664525u + 1013904223u;
    return *s;
}
static inline double lcg_uniform(unsigned int *s)
{
    return (lcg_next(s) & 0x7FFFFFFFu) / (double)0x7FFFFFFFu;
}

/* ─── Box-Muller: produce one Gaussian sample, clipped to [0,1] ────────── */
static double gauss_sample(unsigned int *s)
{
    /* consume two uniforms; cache second in a static — keep it simple */
    static __thread int have_spare = 0;
    static __thread double spare = 0.0;

    double z;
    if (have_spare)
    {
        z = spare;
        have_spare = 0;
    }
    else
    {
        double u, v, r;
        do
        {
            u = 2.0 * lcg_uniform(s) - 1.0;
            v = 2.0 * lcg_uniform(s) - 1.0;
            r = u * u + v * v;
        } while (r >= 1.0 || r == 0.0);
        double fac = sqrt(-2.0 * log(r) / r);
        spare = v * fac;
        have_spare = 1;
        z = u * fac;
    }

    double val = GAUSS_MEAN + GAUSS_SIGMA * z;
    /* clip to [0,1] */
    if (val < 0.0)
        val = 0.0;
    if (val > 1.0)
        val = 1.0;
    return val;
}

typedef struct
{
    int thread_id;
    int start_sec;
    int end_sec;
} ThreadData;

THREAD_FUNC worker(void *arg)
{
    ThreadData *t = (ThreadData *)arg;
    unsigned int seed = FIXED_SEED + (unsigned int)t->thread_id * 2654435761u;

    for (int sec = t->start_sec; sec < t->end_sec; sec++)
    {
        int minute = sec / 60;
        int ten_min = sec / 600;

        for (long i = 0; i < VALUES_PER_SEC; i++)
        {
            double val = gauss_sample(&seed);
            int bin = (int)(val * BINS);
            if (bin >= BINS)
                bin = BINS - 1;

            ATOMIC_ADD(&global_hist[bin], 1L);
            ATOMIC_ADD(&MIN_HIST(minute, bin), 1L);
            ATOMIC_ADD(&TEN_HIST(ten_min, bin), 1L);
        }
    }
#ifdef _WIN32
    return 0;
#else
    return NULL;
#endif
}

void compute_stats(long *hist, long total,
                   double *mean, double *min_v, double *max_v,
                   double *median, double *p25, double *p75,
                   double *mode, double *stddev)
{
    long cumul = 0;
    double sum = 0.0, sum_sq = 0.0;

    *min_v = -1.0;
    *max_v = -1.0;
    *median = *p25 = *p75 = 0.0;

    long max_cnt = 0;
    int mode_bin = 0;

    for (int i = 0; i < BINS; i++)
    {
        if (!hist[i])
            continue;
        double bc = (i + 0.5) / (double)BINS;
        if (*min_v < 0)
            *min_v = bc;
        *max_v = bc;
        sum += hist[i] * bc;
        sum_sq += hist[i] * bc * bc;
        if (hist[i] > max_cnt)
        {
            max_cnt = hist[i];
            mode_bin = i;
        }
    }

    *mean = sum / (double)total;
    *mode = (mode_bin + 0.5) / (double)BINS;
    *stddev = sqrt(sum_sq / (double)total - (*mean) * (*mean));

    long q25 = (long)(total * 0.25);
    long q50 = (long)(total * 0.50);
    long q75 = (long)(total * 0.75);

    for (int i = 0; i < BINS; i++)
    {
        cumul += hist[i];
        double bc = (i + 0.5) / (double)BINS;
        if (*p25 == 0.0 && cumul >= q25)
            *p25 = bc;
        if (*median == 0.0 && cumul >= q50)
            *median = bc;
        if (*p75 == 0.0 && cumul >= q75)
            *p75 = bc;
    }
}

void print_stats(const char *label,
                 double mean, double min_v, double max_v,
                 double median, double p25, double p75,
                 double mode, double stddev)
{
    printf("%-22s Mean=%7.5f  StdDev=%7.5f  Min=%7.5f  Max=%7.5f  "
           "Median=%7.5f  P25=%7.5f  P75=%7.5f  Mode=%7.5f\n",
           label, mean, stddev, min_v, max_v, median, p25, p75, mode);
}

int main(void)
{

    printf("╔════════════════════════════════════════════════════════════╗\n");
    printf("║   Gaussian Stream Statistics Simulator  (mean=0.5, σ=0.15) ║\n");
    printf("╚════════════════════════════════════════════════════════════╝\n\n");
    printf("  Duration       : %d s  (1 hour)\n", TOTAL_SECONDS);
    printf("  Values/second  : %ld\n", VALUES_PER_SEC);
    printf("  Total values   : %ld\n", TOTAL_VALUES);
    printf("  Minute windows : %d\n", TOTAL_MINUTES);
    printf("  10-min blocks  : %d\n", TEN_MIN_BLOCKS);
    printf("  Worker threads : %d\n", THREADS);
    printf("  Fixed seed     : %u\n\n", FIXED_SEED);

    memset(global_hist, 0, sizeof(global_hist));
    minute_hist = calloc((size_t)TOTAL_MINUTES * BINS, sizeof(long));
    ten_min_hist = calloc((size_t)TEN_MIN_BLOCKS * BINS, sizeof(long));
    if (!minute_hist || !ten_min_hist)
    {
        fprintf(stderr, "Memory allocation failed.\n");
        return 1;
    }

    printf("  Running simulation (Box-Muller Gaussian)...");
    fflush(stdout);

    thread_t threads[THREADS];
    ThreadData tdata[THREADS];
    int chunk = TOTAL_SECONDS / THREADS;

    for (int i = 0; i < THREADS; i++)
    {
        tdata[i].thread_id = i;
        tdata[i].start_sec = i * chunk;
        tdata[i].end_sec = (i == THREADS - 1) ? TOTAL_SECONDS : (i + 1) * chunk;
        thread_create(&threads[i], worker, &tdata[i]);
    }
    for (int i = 0; i < THREADS; i++)
        thread_join(threads[i]);

    printf(" done.\n\n");

    double mean, min_v, max_v, median, p25, p75, mode, stddev;

    compute_stats(global_hist, TOTAL_VALUES,
                  &mean, &min_v, &max_v, &median, &p25, &p75, &mode, &stddev);

    printf("╔══════════════════════════════════════════════════════╗\n");
    printf("║           PART (a): OVERALL STATISTICS               ║\n");
    printf("╚══════════════════════════════════════════════════════╝\n");

    print_stats("Overall:", mean, min_v, max_v, median, p25, p75, mode, stddev);

    printf("\n  Theoretical values for Gaussian (mu=0.5, sigma=0.15) clipped to [0,1]:\n");
    printf("    Mean         ≈ 0.50000   got %.5f\n", mean);
    printf("    Std Dev      ≈ 0.14xxx   got %.5f  (clipping reduces sigma slightly)\n", stddev);
    printf("    Median       ≈ 0.50000   got %.5f\n", median);
    printf("    P25          ≈ 0.39878   got %.5f  (mu - 0.6745*sigma)\n", p25);
    printf("    P75          ≈ 0.60122   got %.5f  (mu + 0.6745*sigma)\n", p75);
    printf("    Mode         ≈ 0.50000   got %.5f  (peak of bell curve)\n", mode);

    printf("\n  Significance of measures:\n");
    printf("    Mean   : Central tendency; equals median for a symmetric Gaussian.\n");
    printf("    Median : Robust to outliers; confirms near-perfect symmetry here.\n");
    printf("    Mode   : Peak of the bell curve; coincides with mean/median.\n");
    printf("    StdDev : Spread of the distribution; ~68%% of data in [mu-sigma, mu+sigma].\n");
    printf("    P25/P75: IQR region; for Gaussian, IQR ≈ 1.3490 * sigma.\n");
    printf("    Min/Max: Empirical bounds; clipping at 0 and 1 creates small\n");
    printf("             pile-up at boundaries (visible in histogram tails).\n");
    printf("╔═══════════════════════════════════════════════════════╗\n");
    printf("║         PART (b): PER-MINUTE ANALYSIS                  ║\n");
    printf("╚════════════════════════════════════════════════════════╝\n");
    printf("%-22s %-9s %-9s %-9s %-9s %-9s %-9s %-9s %-9s\n",
           "Interval", "Mean", "StdDev", "Min", "Max", "Median", "P25", "P75", "Mode");

    FILE *csv_min = fopen("minute_stats.csv", "w");
    fprintf(csv_min, "minute,mean,stddev,min,max,median,p25,p75,mode\n");

    /* Anomaly detection: flag minute if |mean - 0.5| > 3*sigma/sqrt(N) */
    long N_per_min = VALUES_PER_MIN;
    double expected_se = GAUSS_SIGMA / sqrt((double)N_per_min);
    int anomaly_found = 0;

    for (int m = 0; m < TOTAL_MINUTES; m++)
    {
        double mn, mi, ma, med, q1, q3, mo, sd;
        compute_stats(&MIN_HIST(m, 0), VALUES_PER_MIN,
                      &mn, &mi, &ma, &med, &q1, &q3, &mo, &sd);

        char label[32];
        snprintf(label, sizeof(label), "Minute %02d:", m + 1);
        print_stats(label, mn, mi, ma, med, q1, q3, mo, sd);

        /* Anomaly flag */
        if (fabs(mn - 0.5) > 3.0 * expected_se)
        {
            printf("  *** ANOMALY detected in Minute %02d: mean=%.5f deviates "
                   "by >3 SE from 0.5 ***\n",
                   m + 1, mn);
            anomaly_found = 1;
        }

        fprintf(csv_min, "%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n",
                m + 1, mn, sd, mi, ma, med, q1, q3, mo);
    }
    fclose(csv_min);

    if (!anomaly_found)
        printf("\n  No anomalies detected (all minute means within 3 SE of 0.5).\n");

    printf("\n  Trend interpretation:\n");
    printf("    Because values are i.i.d. Gaussian, per-minute means\n");
    printf("    should hover at ~0.5 with standard error ≈ %.7f.\n", expected_se);
    printf("    Any minute deviating by >3 SE is flagged as anomalous.\n");
    printf("    Std-dev per minute ≈ %.5f (mirrors population sigma).\n", GAUSS_SIGMA);
    printf("╔════════════════════════════════════════════════════╗\n");
    printf("║         10-MINUTE BLOCK ANALYSIS                     ║\n");
    printf("╚══════════════════════════════════════════════════════╝\n");
    printf("%-22s %-9s %-9s %-9s %-9s %-9s %-9s %-9s %-9s\n",
           "Interval", "Mean", "StdDev", "Min", "Max", "Median", "P25", "P75", "Mode");

    FILE *csv_10 = fopen("ten_min_stats.csv", "w");
    fprintf(csv_10, "block,mean,stddev,min,max,median,p25,p75,mode\n");

    for (int b = 0; b < TEN_MIN_BLOCKS; b++)
    {
        double mn, mi, ma, med, q1, q3, mo, sd;
        compute_stats(&TEN_HIST(b, 0), VALUES_PER_10M,
                      &mn, &mi, &ma, &med, &q1, &q3, &mo, &sd);

        char label[32];
        snprintf(label, sizeof(label), "Min %02d-%02d:", b * 10 + 1, (b + 1) * 10);
        print_stats(label, mn, mi, ma, med, q1, q3, mo, sd);

        fprintf(csv_10, "%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n",
                b + 1, mn, sd, mi, ma, med, q1, q3, mo);
    }
    fclose(csv_10);

    compute_stats(global_hist, TOTAL_VALUES,
                  &mean, &min_v, &max_v, &median, &p25, &p75, &mode, &stddev);

    double IQR = p75 - p25;
    double lower = p25 - 1.5 * IQR;
    double upper = p75 + 1.5 * IQR;

    long outlier_count = 0;
    for (int i = 0; i < BINS; i++)
    {
        double bc = (i + 0.5) / (double)BINS;
        if (bc < lower || bc > upper)
            outlier_count += global_hist[i];
    }

    printf("╔══════════════════════════════════════════════════════╗\n");
    printf("║       PART (c): IQR & OUTLIER ANALYSIS               ║\n");
    printf("╚══════════════════════════════════════════════════════╝\n");
    printf("  Q1 (P25)     : %.5f\n", p25);
    printf("  Q3 (P75)     : %.5f\n", p75);
    printf("  IQR          : %.5f   (expected ≈ %.5f = 1.3490 × sigma)\n",
           IQR, 1.3490 * GAUSS_SIGMA);
    printf("  Lower fence  : %.5f   (Q1 - 1.5 × IQR)\n", lower);
    printf("  Upper fence  : %.5f   (Q3 + 1.5 × IQR)\n", upper);
    printf("  Outlier range: values < %.5f  OR  > %.5f\n", lower, upper);
    printf("  Outlier count: %ld  (%.4f%% of %ld values)\n",
           outlier_count,
           100.0 * outlier_count / (double)TOTAL_VALUES,
           TOTAL_VALUES);
    printf("\n  For a pure Gaussian, Tukey fences exclude ≈ 0.70%% of the tail.\n");
    printf("  Clipping at [0,1] transfers some tail mass onto the boundary bins,\n");
    printf("  slightly inflating outlier counts near 0 and 1.\n");
    printf("\n  Impact of outliers on central tendency:\n");
    printf("    Mean   : Mildly pulled toward 0 or 1 by clipped values.\n");
    printf("    Median : Virtually unaffected — robust to tails.\n");
    printf("    Mode   : Unaffected (peak is far from tails).\n");
    printf("    StdDev : Slightly reduced because clipping removes extreme tails.\n");

    FILE *csv_g = fopen("global_stats.csv", "w");
    fprintf(csv_g, "interval,mean,stddev,min,max,median,p25,p75,mode\n");
    fprintf(csv_g, "Overall,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n",
            mean, stddev, min_v, max_v, median, p25, p75, mode);
    for (int b = 0; b < TEN_MIN_BLOCKS; b++)
    {
        double mn, mi, ma, med, q1, q3, mo, sd;
        compute_stats(&TEN_HIST(b, 0), VALUES_PER_10M,
                      &mn, &mi, &ma, &med, &q1, &q3, &mo, &sd);
        fprintf(csv_g, "Block_%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n",
                b + 1, mn, sd, mi, ma, med, q1, q3, mo);
    }
    fclose(csv_g);

    /* Config for the plotter */
    FILE *cfg_f = fopen("sim_config.txt", "w");
    fprintf(cfg_f, "total_seconds=%d\n", TOTAL_SECONDS);
    fprintf(cfg_f, "total_minutes=%d\n", TOTAL_MINUTES);
    fprintf(cfg_f, "ten_min_blocks=%d\n", TEN_MIN_BLOCKS);
    fprintf(cfg_f, "values_per_sec=%ld\n", VALUES_PER_SEC);
    fprintf(cfg_f, "total_values=%ld\n", TOTAL_VALUES);
    fprintf(cfg_f, "gauss_mean=%.4f\n", GAUSS_MEAN);
    fprintf(cfg_f, "gauss_sigma=%.4f\n", GAUSS_SIGMA);
    fclose(cfg_f);

    printf("\n  CSV files saved: global_stats.csv  minute_stats.csv  ten_min_stats.csv\n");
    printf("  Run plot_boxplots.py to generate box plots.\n\n");

    free(minute_hist);
    free(ten_min_hist);
    return 0;
}