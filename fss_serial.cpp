#include "fss.h"
#include <cmath>
#include <iostream>

#define MIN(x, y) (x < y ? x : y)

int t;
double cfvalnx, cfvalny, cfvald; // numerator and denominator of collective movement

// benchmark function
double f(double x, double y) {
	// return sin(sqrt(x * x + y * y));
    // return ((x - 3.5) * (x-3.5))+ ((y - 3.5) * (y - 3.5));
    return sin(0.01 * x * x + 0.005 * y * y - 0.05 * x + 2 * sin(0.01 * t));
}

void individual_move(fish_t& i) {
    i.fval = f(i.x, i.y);
    double next_fval = f(i.x + i.vx * dt, i.y + i.vy * dt); // units mismatch
    if (next_fval < i.fval) {   // fval at step 0 initialized in main
    	i.ax = - i.vx / dt;  // Setting velocity to 0
    	i.ay = - i.vy / dt;
        i.vx = i.vy = 0;
    }
    else {
    	i.ax = i.ay = 0;
    }

    cfvalnx += i.vx * (next_fval - i.fval);
    cfvalny += i.vy * (next_fval - i.fval);
    cfvald += fabs(next_fval - i.fval);

    // i.weight += (next_fval - i.fval) / fabs(next_fval - i.fval); // problematic
    // i.weight = fmax(1, i.weight);
    // i.weight = fmin(Wscale, i.weight);
}

void collective_move(fish_t& i) {
    if (cfvald < 1e-10) std::cout << "1\n";
	i.ax += cfvalnx / cfvald;
	i.ay += cfvalny / cfvald;
}

// Integrate the ODE
void move_fish(fish_t& i, double size) {
    // Slightly simplified Velocity Verlet integration
    // Conserves energy better than explicit Euler method
    i.vx += i.ax * dt;
    i.vx = MIN(i.vx, vmax);
    i.vy += i.ay * dt;
    i.vy = MIN(i.vy, vmax);
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

    for (int i = 0; i < nfish; ++ i) {
    	collective_move(fish[i]);
    }

    // Move fish
    for (int i = 0; i < nfish; ++i) {
        move_fish(fish[i], size);
    }
}