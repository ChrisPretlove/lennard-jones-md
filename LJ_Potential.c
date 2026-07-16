#include "stdio.h"
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>
#include <math.h>

// Physical constraints of the simulation enviroment
#define BOX_WIDTH 10.0e-9 // Physical width of the 2D simulation box (10 nanometers)
#define BOX_HEIGHT 10.0e-9 // Physical height of the 2D simulation box (10 nanometers)
#define NUM_PARTICLES 200 // Total number of Argon atoms in the system

// Structure representing the 2D phase-space of a single Argon atom
struct atom
{
    double x, y; // Position coordinates (m)
    double vx, vy; // Velocity componants (m/s)
    double fx, fy; // Accumulator for forces acting on the atom (N)
};

int main()
{
    // Seed the random number generator for initial position and velocity spawning using the current real life time
    srand((unsigned int)time(NULL));

    // Open output files to record physical trajectories and system energy metrics
    FILE *data_ptr = fopen("LJ_data.csv", "w");
    FILE *energy_ptr = fopen("LJ_energy.csv", "w");
    if (data_ptr == NULL) 
    {
        printf("Error opening data file!\n");
        return 1;
    }
    fprintf(data_ptr, "time,body_id,x,y\n");

    if (energy_ptr == NULL) 
    {
        printf("Error opening energy file!\n");
        return 1;
    }
    fprintf(energy_ptr, "time,Kinetic energy,Potential energy,Total energy\n");

    // PHYSICAL CONSTANTS & SIMULATION PARAMETERS (Argon System)
    struct atom particles[NUM_PARTICLES];
    double epsilon = 1.65e-21; // Well-depth energy of Argon (Joules)
    double sigma = 3.4e-10; // Collision diameter of Argon (meters)
    double sigma6 = pow(sigma,6);
    double sigma12 = sigma6 * sigma6;
    double dt = 2.0e-15; // Integration time step (seconds, 2 femtoseconds)
    double r_cutoff = 2.5 * sigma; // Force cutoff distance (meters) to optimise performance
    
    // Shift potential at cutoff to prevent a step change (discontinuity) in energy at r = r_cutoff
    double cutoff_potential = 4 * epsilon * (sigma12 / pow(r_cutoff,12) - sigma6 / pow(r_cutoff,6));
    double m = 6.63e-26; // Mass of an Argon atom (kg)

    // INITIALISATION: Spawning Particles and Assigning Velocities
    // Initialize the first particle's velocity components randomly
    particles[0].vx = (-100.0 + ((double)rand() / RAND_MAX) * 200.0);
    particles[0].vy = (-100.0 + ((double)rand() / RAND_MAX) * 200.0);
    particles[0].vx = -1 + (double)rand() / RAND_MAX * (1 - -1);
    particles[0].vy = -1 + (double)rand() / RAND_MAX * (1 - -1);
    
    // Spawn the remaining particles, ensuring no physical overlapping (clashing)
    for (int i = 1; i < NUM_PARTICLES; i++)
    {
        bool clash = true;
        while (clash == true)
        {
            clash = false;
            // Generate random spatial coordinates within the boundaries of the box
            particles[i].x = ((double)rand()/(double)(RAND_MAX)) * BOX_WIDTH;
            particles[i].y = ((double)rand()/(double)(RAND_MAX)) * BOX_HEIGHT;
            
            // Check if this new coordinate overlaps with any existing particles
            for (int j = 0; j < i; j++)
            {
                double dx = particles[i].x - particles[j].x;
                double dy = particles[i].y - particles[j].y;
                double r = sqrt(dx*dx + dy*dy);

                // If particles are overlapping (closer than 95% of their atomic diameter), trigger a clash
                if (r < 0.95 * sigma) 
                {
                    clash = true;
                    break;
                }
            }
        }
        // Assign random initial velocities to mimic thermal motion
        particles[i].vx = (-100.0 + ((double)rand() / RAND_MAX) * 200.0);
        particles[i].vy = (-100.0 + ((double)rand() / RAND_MAX) * 200.0);
    }

    //INITIAL FORCE CALCULATION (Pre-loop setup)
    // Initialize all force vectors to zero
    for (int _ = 0; _ < NUM_PARTICLES; _ ++)
    {
        particles[_].fx = 0;
        particles[_].fy = 0; 
    }

    // Compute initial pairwise Lennard-Jones forces prior to starting integration
    for (int i = 0; i < NUM_PARTICLES; i++)
    {
        for (int j = i + 1; j < NUM_PARTICLES; j++)
        {
            double dx = particles[j].x - particles[i].x;
            double dy = particles[j].y - particles[i].y;
            double r_2 = dx * dx + dy * dy;

            // Apply neighbor list cutoff to ignore interactions beyond the cutoff radius
            if (r_2 < r_cutoff*r_cutoff)
            {
                double r_6 = r_2 * r_2 * r_2;
                double r_12 = r_6 * r_6;

                // Analytical force magnitude from the LJ potential gradient: F(r) = -dV/dr
                double F = (24 * epsilon / r_2) * ((2*sigma12 / r_12) - (sigma6/r_6));

                // Apply Newton's third law to assign equal and opposite pairwise force components
                particles[i].fx -= F*dx;
                particles[i].fy -= F*dy;
                particles[j].fx += F*dx;
                particles[j].fy += F*dy;
            }
            else
            {
                continue;
            }
        }
    }
    
    // MAIN MOLECULAR DYNAMICS INTEGRATION LOOP (Velocity Verlet)
    double t_max = 50.0e-12; // Run simulation for 50 picoseconds
    for (double t = 0.0; t < t_max; t += dt)
    {
        double PE = 0.0; // Potential Energy accumulator
        double KE = 0.0; // Kinetic Energy accumulator

        // Record physical trajectories to data file for Python visualisation
        for (int i = 0; i < NUM_PARTICLES; i++)
        {
            fprintf(data_ptr, "%e,%d,%e,%e\n", t, i, particles[i].x, particles[i].y);
        }

        // 1. First half-step velocity update
        // v(t + dt/2) = v(t) + 0.5 * a(t) * dt
        for (int i = 0; i < NUM_PARTICLES; i++)
        {
            particles[i].vx += 0.5 * particles[i].fx / m * dt;
            particles[i].vy += 0.5 * particles[i].fy / m * dt;
        }

        // 2. Full-step position update and wall reflections
        // r(t + dt) = r(t) + v(t + dt/2) * dt
        for (int i = 0; i < NUM_PARTICLES; i++)
        {
            particles[i].x += particles[i].vx * dt;
            particles[i].y += particles[i].vy * dt;

            // Reflective Boundary Conditions: Maintain physical boundaries with elastic wall collisions
            if (particles[i].x < 0)
            {
                particles[i].x = -particles[i].x;
                particles[i].vx = -particles[i].vx;
            }
            else if (particles[i].x > BOX_WIDTH)
            {
                particles[i].x = 2 * BOX_WIDTH - particles[i].x;
                particles[i].vx = -particles[i].vx;
            }

            if (particles[i].y < 0)
            {
                particles[i].y = -particles[i].y;
                particles[i].vy = -particles[i].vy;
            }
            else if (particles[i].y > BOX_HEIGHT)
            {
                particles[i].y = 2 * BOX_HEIGHT - particles[i].y;
                particles[i].vy = -particles[i].vy;
            }
        }

        // 3. Clear forces
        // Prepare accumulators to compute fresh forces at the new spatial positions
        for (int i = 0; i < NUM_PARTICLES; i++)
        {
            particles[i].fx = 0.0;
            particles[i].fy = 0.0;
        }

        // 4. Compute new pair forces (at the new positions)
        for (int i = 0; i < NUM_PARTICLES; i++)
        {
            for (int j = i + 1; j < NUM_PARTICLES; j++)
            {
                double dx = particles[j].x - particles[i].x;
                double dy = particles[j].y - particles[i].y;
                double r2 = dx*dx + dy*dy;

                if (r2 < r_cutoff*r_cutoff)
                {
                    double r6 = r2*r2*r2;
                    double r12 = r6*r6;

                    // Compute force vector magnitude
                    double F = (24.0 * epsilon / r2) * ((2.0 * sigma12 / r12) - (sigma6 / r6));
                    
                    // Accumulate Potential Energy, subtracting cutoff potential to maintain smooth energy metrics
                    PE += 4 * epsilon * (sigma12 / r12 - sigma6 / r6) - cutoff_potential;

                    // Distribute force components using Newton's third law
                    particles[i].fx -= F * dx;
                    particles[i].fy -= F * dy;
                    particles[j].fx += F * dx;
                    particles[j].fy += F * dy;
                }
                else
                {
                    continue;
                }
            }
        }

        // 5. Second half-step velocty update (with new forces)
        // v(t + dt) = v(t + dt/2) + 0.5 * a(t + dt) * dt
        for (int i = 0; i < NUM_PARTICLES; i++)
        {
            particles[i].vx += 0.5 * particles[i].fx / m * dt;
            particles[i].vy += 0.5 * particles[i].fy / m * dt;
            
            // Accumulate kinetic energy: KE = 0.5 * m * v^2
            KE += 0.5 * m * (particles[i].vx * particles[i].vx + particles[i].vy * particles[i].vy);
        }

        // Log simulation time, kinetic, potential, and total energy (conserved quantity)
        fprintf(energy_ptr, "%e,%e,%e,%e\n", t, KE, PE, KE+PE);
    }

    // Clean up file connections
    fclose(data_ptr);
    fclose(energy_ptr);
    printf("Simulation complete. Atomic path data saved to LJ_data.csv and energy saved to LJ_energy.csv\n");
    return 0;
}