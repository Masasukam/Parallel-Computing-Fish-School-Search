#include "fss_mpi.h"
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>
#include <mpi.h>

// Define constants and variables for MPI communication
#define MPI_MASTER 0
#define TAG_FISH_DATA 1

// Function prototypes
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

    // int sx = (int)ceil(sqrt((double)nfish));
    // int sy = (nfish + sx - 1) / sx;

    std::vector<int> shuffle(nfish);
    for (int i = 0; i < (int)shuffle.size(); ++i) {
        shuffle[i] = i;
    }

    for (int i = 0; i < nfish; ++i) {
        // Make sure particles are not spatially sorted
        std::uniform_int_distribution<int> rand_int(0, nfish - i - 1);
        int j = rand_int(gen);
        // int k = shuffle[j];
        shuffle[j] = shuffle[nfish - i - 1];

        // Distribute particles evenly to ensure proper spacing
        std::uniform_real_distribution<float> rand_pos(0, 1.0 * size);
        fish[i].x = rand_pos(gen);
        fish[i].y = rand_pos(gen);
        // fish[i].x = size * (1. + (k % sx)) / (1 + sx);
        // fish[i].y = size * (1. + (k / sx)) / (1 + sy);
        // fish[i].weight = Wscale / 2;

        // Assign random velocities within a bound
        std::uniform_real_distribution<float> rand_real(-1.0, 1.0);
        fish[i].vx = rand_real(gen);
        fish[i].vy = rand_real(gen);

        // printf("fish %d: vx: %.2f   vy: %.2f\n", i, fish[i].vx, fish[i].vy); // for debugging
    }

    for (int i = 0; i < nfish; ++i) {
        fish[i].id = i + 1;
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

MPI_Datatype FISH;

int main(int argc, char** argv) {

    // Parse Args
    if (find_arg_idx(argc, argv, "-h") >= 0) {
        //if (rank == MPI_MASTER) {
        std::cout << "Options:" << std::endl;
        std::cout << "-h: see this help" << std::endl;
        std::cout << "-n <int>: set number of fish" << std::endl;
        std::cout << "-o <filename>: set the output file name" << std::endl;
        std::cout << "-s <int>: set fish initialization seed" << std::endl;
        //}
        return 0;
    }

    // Open Output File
    char* savename = find_string_option(argc, argv, "-o", nullptr);
    std::ofstream fsave(savename);

    int num_procs, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Create MPI Particle Type
    const int nitems = 9;
    int blocklengths[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1};
    MPI_Datatype types[9] = {MPI_UINT64_T, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
                             MPI_DOUBLE,   MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
    MPI_Aint offsets[9];
    offsets[0] = offsetof(fish_t, id);
    offsets[1] = offsetof(fish_t, x);
    offsets[2] = offsetof(fish_t, y);
    offsets[3] = offsetof(fish_t, vx);
    offsets[4] = offsetof(fish_t, vy);
    offsets[5] = offsetof(fish_t, ax);
    offsets[6] = offsetof(fish_t, ay);
    offsets[7] = offsetof(fish_t, fval);
    offsets[8] = offsetof(fish_t, weight);
    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &FISH);
    MPI_Type_commit(&FISH);

    // Initialize Particles
    int nfish = find_int_arg(argc, argv, "-n", 1000);
    int fish_seed = find_int_arg(argc, argv, "-s", 0);
    double size = sqrt(density * nfish);

    fish_t* fish = new fish_t[nfish];

    if (rank == 0) {
        init_fish(fish, nfish, size, fish_seed);
    }

    MPI_Bcast(fish, nfish, FISH, 0, MPI_COMM_WORLD);

    // Algorithm
    auto start_time = std::chrono::steady_clock::now();

    init_simulation(fish, nfish, size, rank, num_procs);

    for (int step = 0; step < nsteps; ++step) {
        simulate_one_step(fish, nfish, size, rank, num_procs);

        // Save state if necessary
        if (fsave.good() && (step % savefreq) == 0) {
            gather_for_save(fish, nfish, size, rank, num_procs);
            if (rank == MPI_MASTER) {
                save(fsave, fish, nfish, size);
            }
        }
    }

    auto end_time = std::chrono::steady_clock::now();
    std::chrono::duration<double> diff = end_time - start_time;
    double seconds = diff.count();

    // Finalize
    if (rank == 0) {
        std::cout << "Simulation Time = " << seconds << " seconds for " << nfish << " fish." << std::endl;
        std::cout << "Board Size = " << size << std::endl;
    }
    if (fsave) {
        fsave.close();
    }
    delete[] fish;
    MPI_Finalize();
    return 0;
}