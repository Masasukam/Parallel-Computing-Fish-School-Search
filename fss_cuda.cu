#include "fss.h"
#include <cmath>
#include <random>
#include <iostream>
#include <cuda.h>
#include <curand_kernel.h>

#define THRUST_IGNORE_DEPRECATED_CPP_DIALECT
#define CUB_IGNORE_DEPRECATED_CPP_DIALECT

#define NUM_THREADS 256

#define MIN(x, y) (x < y ? x : y)
#define SQR(x) ((x) * (x))

// when set to 1, the program will print some info for debugging
#define debug_mode 0

int t;
int blks; // number of blocks
double cfvalnx, cfvalny, cfvald; // numerator and denominator of collective movement

// benchmark function
__device__
double f(double x, double y, int t) {
    // return sin(sqrt(x * x + y * y));
    // return - ((x - 5) * (x - 5)) - ((y - 5) * (y - 5));
    return (-1 + 2.0 * t / nsteps) * (((x - 5) * (x - 5)) + ((y - 5) * (y - 5)));
    // return sin(0.01 * x * x + 0.005 * y * y - 0.05 * x + 2 * sin(0.01 * t));
    // return sin(0.01 * x * x + 0.005 * y * y - 0.05 * x + 2 * sin(0.005 * t));
}

// Function to calculate distance between two fishes
__device__
double distanceSquared(const fish_t &fish1, const fish_t &fish2) {
    return SQR(fish1.x - fish2.x) + SQR(fish1.y - fish2.y);
}

// Todo0: individual fish distances
__global__
void individual_move(fish_t* fish, int nfish, int t, double &cfvalnx, double &cfvalny, double &cfvald) {
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= nfish) return;

    curandState state; // random number generator
    curand_init(clock64(), tid, 0, &state);

    fish[tid].fval = f(fish[tid].x, fish[tid].y, t);
    fish[tid].ax = fish[tid].ay = 0;
    
    double next_fval = f(fish[tid].x, fish[tid].y, t);
    double newvx, newvy;
    for (int j = 0; j < retrynum; ++ j) {
        newvx = curand_uniform(&state) * 2.0f - 1.0f;
        newvy = curand_uniform(&state) * 2.0f - 1.0f;
        if (f(fish[tid].x + newvx * dt, fish[tid].y + newvy * dt, t) > next_fval) {
            printf(".\n");
            next_fval = f(fish[tid].x + fish[tid].vx * dt, fish[tid].y + fish[tid].vy * dt, t);
            fish[tid].vx = newvx;
            fish[tid].vy = newvy;
        }
    }

    cfvalnx += fish[tid].vx * (next_fval - fish[tid].fval);
    cfvalny += fish[tid].vy * (next_fval - fish[tid].fval);
    cfvald += fabs(next_fval - fish[tid].fval);

    // Todo 1: Add weight features
    // fish[tid].weight += (next_fval - fish[tid].fval) / fabs(next_fval - fish[tid].fval); // problematic
    // fish[tid].weight = fmax(1, fish[tid].weight);
    // fish[tid].weight = fmin(Wscale, fish[tid].weight);
}


// Todo2: within group distance
__global__
void collective_move(fish_t* fish, int nfish) {
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= nfish) return;

    double cutoffSquared = SQR(cutoff);
    double totalAdjustmentX = 0.0 - cutoff;
    double totalAdjustmentY = 0.0 - cutoff;
    int count = -1;

    for (int j = 0; j < nfish; ++ j) {
        double distSquared = distanceSquared(fish[tid], fish[j]);
        if (distSquared < cutoffSquared) {
            totalAdjustmentX += cutoff - fish[tid].x + fish[j].x;
            totalAdjustmentY += cutoff - fish[tid].y + fish[j].y;
            count ++ ;
        }
    }

    if (count > 0) {
        fish[tid].ax += totalAdjustmentX / count;
        fish[tid].ay += totalAdjustmentY / count;
    }
}

// Integrate the ODE
__global__
void move_fish(fish_t* fish, int nfish, double size) {
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= nfish) return;

    fish[tid].vx += fish[tid].ax * dt;
    fish[tid].vx = MIN(fish[tid].vx, vmax);
    fish[tid].vy += fish[tid].ay * dt;
    fish[tid].vy = MIN(fish[tid].vy, vmax);

    fish[tid].x += fish[tid].vx * dt;
    fish[tid].y += fish[tid].vy * dt;


    // Bounce from walls
    while (fish[tid].x < 0 || fish[tid].x > size) {
        fish[tid].x = fish[tid].x < 0 ? -fish[tid].x : 2 * size - fish[tid].x;
        fish[tid].vx = -fish[tid].vx;
    }

    while (fish[tid].y < 0 || fish[tid].y > size) {
        fish[tid].y = fish[tid].y < 0 ? -fish[tid].y : 2 * size - fish[tid].y;
        fish[tid].vy = -fish[tid].vy;
    }
}


void init_simulation(fish_t* fish, int nfish, double size) {
	// You can use this space to initialize static, global data objects
    // that you may need. This function will be called once before the
    // algorithm begins. Do not do any fish simulation here
    t = 0;
    blks = (nfish + NUM_THREADS - 1) / NUM_THREADS;

    #if debug_mode & 1
    setvbuf( stdout, NULL, _IONBF, 0 );
    #endif
}

void simulate_one_step(fish_t* fish, int nfish, double size) {
    ++ t;
	cfvalnx = cfvalny = cfvald = 0;

    individual_move<<<blks, NUM_THREADS>>>(fish, nfish, t, cfvalnx, cfvalny, cfvald);

    cudaDeviceSynchronize();

    #if debug_mode & 1
    printf("here 1\n");
    #endif

    collective_move<<<blks, NUM_THREADS>>>(fish, nfish);

    cudaDeviceSynchronize();

    #if debug_mode & 1
    printf("here 2\n");
    #endif

    move_fish<<<blks, NUM_THREADS>>>(fish, nfish, size);

    #if debug_mode & 1
    printf("here 3\n");
    #endif
}