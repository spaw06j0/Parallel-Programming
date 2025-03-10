#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define TOSS 100000000
// area of a circle is pi * (1/2)^2
// area of a square is 1
// pi = 4 * (area of circle / area of square)

int main(int argc, char *argv[]) {
    long long int number_in_circle = 0;
    for (int toss = 0; toss < TOSS; toss ++) {
        double x = (double)rand() / RAND_MAX * 2 - 1;
        double y = (double)rand() / RAND_MAX * 2 - 1;
        double distance_squared = x * x + y * y;
        if (distance_squared <= 1)
            number_in_circle++;
    }
    double pi_estimate = 4 * number_in_circle /(( double ) TOSS);
    printf("pi estimate: %f\n", pi_estimate);
    return 0;
}