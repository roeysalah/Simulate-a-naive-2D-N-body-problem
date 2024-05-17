#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <omp.h>

#define NUM_THREADS 8
#define N 1000 // stars
#define num_steps 10000 
#define timestep 2e15 // rate 
#define G 6.674e-11 // Gravitational constant
#define m_sun 2e30 // mass of sun
#define v_avg 200000 // Average speed
#define domain 9e17 //  
#define PI 3.14159265358979323846 // for random theta

// this struct was made for represent one star in our galaxy
typedef struct {
    double x, y;     // position
    double vx, vy;   // velocity
    double fx, fy;   // force
} Body;

Body galaxy[N]; // this array represent the whole galaxy (all stars)

// this function gets the array of all the stars positions and saves the values from the array into the text file
void snapshots(int step_snap) {
    char filename[20];
    sprintf(filename, "snapshot_t%d.txt", step_snap);
    FILE* file = fopen(filename, "w");
    if (file == NULL) // try to create a file
    {
        printf("Error opening file for writing: %s\n", filename);
        return;
    }
    // fill the file with the stars positions
    int i;
    for (i = 0; i < N; i++) {
        fprintf(file, "%lf %lf\n", galaxy[i].x, galaxy[i].y); // x and y coordinate
    }
    fclose(file); // close this file
}

// this function generates a random double between min and max
double randomDouble(double min, double max) {
    return min + ((double)rand() / RAND_MAX) * (max - min); // generate rand number using rand() function
}

// this function generate initial random positions and velocities, in a square domain.
// velocities in the range 0.5v<v<1.5v with uniform direction distributions
void initialize() {
    int i;
    omp_set_num_threads(NUM_THREADS);
    #pragma omp parallel for private(i) shared(galaxy)
    for (i = 0; i < N; i++) 
    {
        galaxy[i].x = randomDouble(0.0, 1) * domain; // choose a random x position
        galaxy[i].y = randomDouble(0.0, 1) * domain; // choose a random y position
        double v = randomDouble(0.5, 1.5); // range of velocity: 0.5v<v<1.5v
        double theta = randomDouble(0.0, 2.0 * PI);
        galaxy[i].vx = v * cos(theta); // choose a random velocity in X axis 
        galaxy[i].vy = v * sin(theta); // choose a random velocity in Y axis 
    }
}

// this function compute force on every star (name: body) 
void compute_forces() {
    double dx, dy, distance, force;
    int i, j;
    omp_set_num_threads(NUM_THREADS);
    #pragma omp parallel for private(i, j, dx, dy, distance, force) shared(galaxy)
    for (i = 0; i < N; i++) {
        galaxy[i].fx = 0.0;
        galaxy[i].fy = 0.0;
        for (j = 0; j < N; j++) // check which stars influence (in terms of force) the specific i star
        {
            if (j != i) // check that is another star
            {
                dx = galaxy[j].x - galaxy[i].x;
                dy = galaxy[j].y - galaxy[i].y;
                distance = sqrt(dx * dx + dy * dy);
                force = (G * m_sun * m_sun) / (distance * distance);
                galaxy[i].fx += force * (dx / distance);
                galaxy[i].fy += force * (dy / distance);
            }
        }
    }
}

// this function compute and update the positions and the velocities of all the stars
void update_positions_velocities() {
    int i;
    omp_set_num_threads(NUM_THREADS);
    #pragma omp parallel for private(i) shared(galaxy)
    for (i = 0; i < N; i++) {
        // update velocity & position
        galaxy[i].x += galaxy[i].x * timestep; // new velocity in X axis
        galaxy[i].y += galaxy[i].y * timestep; // new velocity in Y axis
        galaxy[i].x += galaxy[i].vx * timestep; // new position in X axis
        galaxy[i].y += galaxy[i].vy * timestep; // new position in Y axis

        // check if a body exits the domain we will choose new position (in the domain)
        if (galaxy[i].x < 0.0) galaxy[i].x = randomDouble(0.0, 1) * domain;
        if (galaxy[i].x > domain) galaxy[i].x = randomDouble(0.0, 1) * domain;
        if (galaxy[i].y < 0.0) galaxy[i].y = randomDouble(0.0, 1) * domain;
        if (galaxy[i].y > domain) galaxy[i].y = randomDouble(0.0, 1) * domain;
    }
}

int main(int argc, char** argv) {
    srand(time(NULL)); // seed for random numbers
    double start_time, end_time;
    Body *bodies;
    initialize(); // init all stars (whole galaxy)
    snapshots(0); // take a snapshot of the initalize galaxy


	// start simulation
    start_time = omp_get_wtime(); // start timer
    int step;
    for (step = 0; step < num_steps; step++) {
        compute_forces(); // compute forces
        update_positions_velocities(); // update positions and velocities
        if (step == num_steps / 2) // if we pass half of the total iterations 
        {
            snapshots(num_steps / 2); // take a snapshot in t= Tmax/2 of the galaxy
        }
    }

    end_time = omp_get_wtime(); // stop timer
    snapshots(num_steps); // take a snapshot in t= Tmax of the galaxy
    double elapsed_time = end_time - start_time; // calculate the total time of the simulation
    printf("Done Simulation!!! Elapsed time: %f seconds\n", elapsed_time); // print the total time of the simulation

    return 0;
}

