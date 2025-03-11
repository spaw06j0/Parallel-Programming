#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <pthread.h>

// #define TOSS 100000000
// #define NUM_THREADS 4

pthread_mutex_t mutexSum;

// area of a circle is pi * (1/2)^2
// area of a square is 1
// pi = 4 * (area of circle / area of square)

typedef struct {
    int start;
    int end;
    long long int *numOfhit;
} Arg_t;

void *pi(void *args) {

    Arg_t *arg = (Arg_t *)args;
    int start = arg->start;
    int end = arg->end;
    long long int *numOfhit = arg->numOfhit;
    long long int local_numOfhit = 0;
    
    unsigned int seed = 777;

    for (int i = start; i < end; i++) {
        double x = ((double) rand_r(&seed) / RAND_MAX) * 2.0 - 1.0;
        double y = ((double) rand_r(&seed) / RAND_MAX) * 2.0 - 1.0;
        
        double distance_squared = x * x + y * y;
        if (distance_squared <= 1)
            local_numOfhit++;
    }

    pthread_mutex_lock(&mutexSum);
    *numOfhit += local_numOfhit;
    pthread_mutex_unlock(&mutexSum);

    pthread_exit((void *)0);
}

int main(int argc, char *argv[]) {
    int NUM_THREADS = atoi(argv[1]);
    long long int TOSS = atoll(argv[2]);
    pthread_t threads[NUM_THREADS];
    Arg_t arg[NUM_THREADS];
    int each_thread_toss = TOSS / NUM_THREADS;
    long long int *numOfhit = (long long int *) malloc(sizeof(*numOfhit));
    *numOfhit = 0;
    pthread_mutex_init(&mutexSum, NULL);

    for (int i = 0; i < NUM_THREADS; i++) {
        arg[i].start = i * each_thread_toss;
        arg[i].end = (i + 1) * each_thread_toss;
        arg[i].numOfhit = numOfhit;

        pthread_create(&threads[i], NULL, pi, (void *)&arg[i]);
    }

    for (int i = 0; i < NUM_THREADS; i++) {
        pthread_join(threads[i], NULL);
    }

    pthread_mutex_destroy(&mutexSum);

    // for (int toss = 0; toss < TOSS; toss ++) {
    //     double x = (double)rand() / RAND_MAX * 2 - 1;
    //     double y = (double)rand() / RAND_MAX * 2 - 1;
    //     double distance_squared = x * x + y * y;
    //     if (distance_squared <= 1)
    //         number_in_circle++;
    // }
    double pi_estimate = 4 * (*numOfhit) /(( double ) TOSS);
    printf("pi estimate: %f\n", pi_estimate);
    return 0;
}