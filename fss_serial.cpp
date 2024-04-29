#include "fss.h"
#include <cmath>
#include <random>
#include <iostream>

#define MIN(x, y) (x < y ? x : y)
#define SQR(x) ((x) * (x))

// random number generator
std::random_device rd;
std::mt19937 gen(1);
std::uniform_real_distribution<float> rand_real(-1.0, 1.0);

int t;
double cfvalnx, cfvalny, cfvald; // numerator and denominator of collective movement

// benchmark function
double f(double x, double y) {
	// return sin(sqrt(x * x + y * y));
    // return ((x - 3.5) * (x-3.5))+ ((y - 3.5) * (y - 3.5));
    // return sin(0.01 * x * x + 0.005 * y * y - 0.05 * x + 2 * sin(0.01 * t));
    return sin(0.01 * x * x + 0.005 * y * y - 0.05 * x + 2 * sin(0.005 * t));
}

// Function to calculate distance between two fishes
double distanceSquared(const fish_t &fish1, const fish_t &fish2) {
    return SQR(fish1.x - fish2.x) + SQR(fish1.y - fish2.y);
}

// Todo0: individual fish distances
void individual_move(fish_t& i) {
    i.fval = f(i.x, i.y);
    i.ax = i.ay = 0;
    
    double next_fval = f(i.x, i.y);
    double newvx, newvy;
    for (int j = 0; j < retrynum; ++ j) {
        newvx = rand_real(gen);
        newvy = rand_real(gen);
        if (f(i.x + newvx * dt, i.y + newvy * dt) > next_fval) {
            next_fval = f(i.x + i.vx * dt, i.y + i.vy * dt);
            i.vx = newvx;
            i.vy = newvy;
        }
    }

    cfvalnx += i.vx * (next_fval - i.fval);
    cfvalny += i.vy * (next_fval - i.fval);
    cfvald += fabs(next_fval - i.fval);

    // Todo 1: Add weight features
    // i.weight += (next_fval - i.fval) / fabs(next_fval - i.fval); // problematic
    // i.weight = fmax(1, i.weight);
    // i.weight = fmin(Wscale, i.weight);
}


// Todo2: within group distance
void collective_move(fish_t& i, fish_t* fishes, int nfish) {
    double cutoffSquared = SQR(cutoff);
    double totalAdjustmentX = 0.0;
    double totalAdjustmentY = 0.0;
    int count = 0;

    for (int j = 0; j < nfish; ++j) {
        if (&i != &fishes[j]) { // Skip self-comparison
            double distSquared = distanceSquared(i, fishes[j]);
            if (distSquared < cutoffSquared) {
                double adjustmentX = i.x - fishes[j].x;
                double adjustmentY = i.y - fishes[j].y;
                totalAdjustmentX += cutoff - adjustmentX;
                totalAdjustmentY += cutoff - adjustmentY;
                count++;
            }
        }
    }

    if (count > 0) {
        i.ax += totalAdjustmentX / count;
        i.ay += totalAdjustmentY / count;
    }
}

// Integrate the ODE
void move_fish(fish_t& i, double size) {
    // Slightly simplified Velocity Verlet integration
    // Conserves energy better than explicit Euler method
    i.vx += i.ax * dt;
    i.vx = MIN(i.vx, vmax);
    i.vy += i.ay * dt;
    i.vy = MIN(i.vy, vmax);

    // make sure that every fish will move
    if (fabs(i.vx) < 1e-3 && fabs(i.vy) < 1e-3) {
        double bestfval = f(i.x, i.y);
        double newvx, newvy;
        for (int j = 0; j < retrynum; ++ j) {
            newvx = rand_real(gen);
            newvy = rand_real(gen);
            if (f(i.x + newvx * dt, i.y + newvy * dt) > bestfval) {
                bestfval = f(i.x + i.vx * dt, i.y + i.vy * dt);
                i.vx = newvx;
                i.vy = newvy;
            }
        }
    }


    i.x += i.vx * dt;
    i.y += i.vy * dt;


    // Bounce from walls
    while (i.x < 0 || i.x > size) {
        i.x = i.x < 0 ? -i.x : 2 * size - i.x;
        i.vx = -i.vx;
    }

    while (i.y < 0 || i.y > size) {
        i.y = i.y < 0 ? -i.y : 2 * size - i.y;
        i.vy = -i.vy;
    }
}


void init_simulation(fish_t* fish, int nfish, double size) {
	// You can use this space to initialize static, global data objects
    // that you may need. This function will be called once before the
    // algorithm begins. Do not do any fish simulation here
    t = 0;
}

void simulate_one_step(fish_t* fish, int nfish, double size) {
    ++ t;
	cfvalnx = cfvalny = cfvald = 0;

    for (int i = 0; i < nfish; ++ i) {
    	individual_move(fish[i]);
    }

    // for (int i = 0; i < nfish; ++ i) {
    // 	collective_move(fish[i]);
    // }
    for (int i = 0; i < nfish; ++i) {
        collective_move(fish[i], fish, nfish);
    }


    // Move fish
    for (int i = 0; i < nfish; ++i) {
        move_fish(fish[i], size);
    }
}