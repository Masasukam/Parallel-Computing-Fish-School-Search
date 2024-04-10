#include "fss.h"
#include <chrono>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>

// =================
// Helper Functions
// =================

// I/O routines
void save(std::ofstream& fsave, fish_t* fish, int nfish, double size) {
    static bool first = true;

    if (first) {
        fsave << nfish << " " << size << "\n";
        first = false;
    }

    for (int i = 0; i < nfish; ++i) {
        fsave << fish[i].x << " " << fish[i].y << "\n";
    }

    fsave << std::endl;
}

// Particle Initialization
void init_fish(fish_t* fish, int nfish, double size, int fish_seed) {
    std::random_device rd;
    std::mt19937 gen(fish_seed ? fish_seed : rd());

    int sx = (int)ceil(sqrt((double)nfish));
    int sy = (nfish + sx - 1) / sx;

    std::vector<int> shuffle(nfish);
    for (int i = 0; i < (int)shuffle.size(); ++i) {
        shuffle[i] = i;
    }

    for (int i = 0; i < nfish; ++i) {
        // Make sure particles are not spatially sorted
        std::uniform_int_distribution<int> rand_int(0, nfish - i - 1);
        int j = rand_int(gen);
        int k = shuffle[j];
        shuffle[j] = shuffle[nfish - i - 1];

        // Distribute particles evenly to ensure proper spacing
        fish[i].x = size * (1. + (k % sx)) / (1 + sx);
        fish[i].y = size * (1. + (k / sx)) / (1 + sy);
        fish[i].fval = f(fish[i].x, fish[i].y);
        fish[i].weight = Wscale / 2;

        // Assign random velocities within a bound
        std::uniform_real_distribution<float> rand_real(-1.0, 1.0);
        fish[i].vx = rand_real(gen);
        fish[i].vy = rand_real(gen);
    }
}

// Command Line Option Processing
int find_arg_idx(int argc, char** argv, const char* option) {
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], option) == 0) {
            return i;
        }
    }
    return -1;
}

int find_int_arg(int argc, char** argv, const char* option, int default_value) {
    int iplace = find_arg_idx(argc, argv, option);

    if (iplace >= 0 && iplace < argc - 1) {
        return std::stoi(argv[iplace + 1]);
    }

    return default_value;
}

char* find_string_option(int argc, char** argv, const char* option, char* default_value) {
    int iplace = find_arg_idx(argc, argv, option);

    if (iplace >= 0 && iplace < argc - 1) {
        return argv[iplace + 1];
    }

    return default_value;
}

// ==============
// Main Function
// ==============

int main(int argc, char** argv) {
    // Parse Args
    if (find_arg_idx(argc, argv, "-h") >= 0) {
        std::cout << "Options:" << std::endl;
        std::cout << "-h: see this help" << std::endl;
        std::cout << "-n <int>: set number of fish" << std::endl;
        std::cout << "-o <filename>: set the output file name" << std::endl;
        std::cout << "-s <int>: set fish initialization seed" << std::endl;
        return 0;
    }

    // Open Output File
    char* savename = find_string_option(argc, argv, "-o", nullptr);
    std::ofstream fsave(savename);

    // Initialize Particles
    int nfish = find_int_arg(argc, argv, "-n", 1000);
    int fish_seed = find_int_arg(argc, argv, "-s", 0);
    double size = sqrt(density * nfish);

    fish_t* fish = new fish_t[nfish];

    init_fish(fish, nfish, size, fish_seed);

    // Algorithm
    auto start_time = std::chrono::steady_clock::now();

    init_simulation(fish, nfish, size);


    for (int step = 0; step < nsteps; ++step) {
        simulate_one_step(fish, nfish, size);

        // Save state if necessary

        if (fsave.good() && (step % savefreq) == 0) {
            save(fsave, fish, nfish, size);
        }
    }

    auto end_time = std::chrono::steady_clock::now();

    std::chrono::duration<double> diff = end_time - start_time;
    double seconds = diff.count();

    // Finalize
    std::cout << "Simulation Time = " << seconds << " seconds for " << nfish << " fish.\n";
    std::cout << "Board Size = " << size << "\n";
    fsave.close();
    delete[] fish;
}
