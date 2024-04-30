#ifndef __FSS_H__
#define __FSS_H__

// Program Constants
#define nsteps   1000
#define savefreq 10
#define density  0.1
// #define Wscale   5000
#define dt       0.02
#define vmax     1.0
#define retrynum 10
#define cutoff   0.1


// Fish Data Structure
typedef struct fish_t {
    double x, y;   // Position
    double vx, vy; // Velocity
    double ax, ay; // Acceleration
    double fval;
    double weight;
} fish_t;

// Benchmark function
double f(double x, double y);

// Simulation routine
void init_simulation(fish_t* fish, int nfish, double size);
void simulate_one_step(fish_t* fish, int nfish, double size);

#endif
