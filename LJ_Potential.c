#include "stdio.h"
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>
#include <math.h>

const float BOX_WIDTH = 5.0e-9;
const float BOX_HEIGHT = 5.0e-9;
const int NUM_PARTICLES = 80;

struct atom
{
    double x, y; //meters
    double vx, vy; //meters per second
    double fx, fy; //newtons (meters squared per second squared)
};

int main()
{
    srand((unsigned int)time(NULL));

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

    struct atom particles[NUM_PARTICLES];
    double epsilon = 1.65e-21;
    double sigma = 3.4e-10;
    double sigma6 = pow(sigma,6);
    double sigma12 = sigma6 * sigma6;
    double dt = 2.0e-15;
    double r_cutoff = 2.5 * sigma;
    double cutoff_potential = 4 * epsilon * (sigma12 / pow(r_cutoff,12) - sigma6 / pow(r_cutoff,6));
    double m = 6.63e-26;

    particles[0].vx = (-100.0 + ((double)rand() / RAND_MAX) * 200.0);
    particles[0].vy = (-100.0 + ((double)rand() / RAND_MAX) * 200.0);
    particles[0].vx = -1 + (double)rand() / RAND_MAX * (1 - -1);
    particles[0].vy = -1 + (double)rand() / RAND_MAX * (1 - -1);
    for (int i = 1; i < NUM_PARTICLES; i++)
    {
        bool clash = true;
        while (clash == true)
        {
            clash = false;
            particles[i].x = ((double)rand()/(double)(RAND_MAX)) * BOX_WIDTH;
            particles[i].y = ((double)rand()/(double)(RAND_MAX)) * BOX_HEIGHT;
            for (int j = 0; j < i; j++)
            {
                double dx = particles[i].x - particles[j].x;
                double dy = particles[i].y - particles[j].y;
                double r = sqrt(dx*dx + dy*dy);

                if (r < 0.95 * sigma) 
                {
                    clash = true;
                    break;
                }
            }
        }
        particles[i].vx = (-100.0 + ((double)rand() / RAND_MAX) * 200.0);
        particles[i].vy = (-100.0 + ((double)rand() / RAND_MAX) * 200.0);
    }

    for (int _ = 0; _ < NUM_PARTICLES; _ ++)
    {
        particles[_].fx = 0;
        particles[_].fy = 0; 
    }
    for (int i = 0; i < NUM_PARTICLES; i++)
    {
        for (int j = i + 1; j < NUM_PARTICLES; j++)
        {
            double dx = particles[j].x - particles[i].x;
            double dy = particles[j].y - particles[i].y;
            double r_2 = dx * dx + dy * dy;
            if (r_2 < r_cutoff*r_cutoff)
            {
                double r_6 = r_2 * r_2 * r_2;
                double r_12 = r_6 * r_6;

                double F = (24 * epsilon / r_2) * ((2*sigma12 / r_12) - (sigma6/r_6));

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
    
    double t_max = 50.0e-12;
    for (double t = 0.0; t < t_max; t += dt)
    {
        double PE = 0.0;
        double KE = 0.0;

        // Save positions
        for (int i = 0; i < NUM_PARTICLES; i++)
        {
            fprintf(data_ptr, "%e,%d,%e,%e\n", t, i, particles[i].x, particles[i].y);
        }

        // 1. FIRST HALF-STEP VELOCITY UPDATE
        for (int i = 0; i < NUM_PARTICLES; i++)
        {
            particles[i].vx += 0.5 * particles[i].fx / m * dt;
            particles[i].vy += 0.5 * particles[i].fy / m * dt;
        }

        // 2. FULL-STEP POSITION UPDATE & WALL REFLECTIONS
        for (int i = 0; i < NUM_PARTICLES; i++)
        {
            particles[i].x += particles[i].vx * dt;
            particles[i].y += particles[i].vy * dt;

            // Reflect off walls
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

        // 3. CLEAR FORCES
        for (int i = 0; i < NUM_PARTICLES; i++)
        {
            particles[i].fx = 0.0;
            particles[i].fy = 0.0;
        }

        // 4. COMPUTE NEW PAIR FORCES (At the new positions)
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

                    double F = (24.0 * epsilon / r2) * ((2.0 * sigma12 / r12) - (sigma6 / r6));
                    PE += 4 * epsilon * (sigma12 / r12 - sigma6 / r6) - cutoff_potential;

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

        // 5. SECOND HALF-STEP VELOCITY UPDATE (With new forces)
        for (int i = 0; i < NUM_PARTICLES; i++)
        {
            particles[i].vx += 0.5 * particles[i].fx / m * dt;
            particles[i].vy += 0.5 * particles[i].fy / m * dt;
            KE += 0.5 * m * (particles[i].vx * particles[i].vx + particles[i].vy * particles[i].vy);
        }
        fprintf(energy_ptr, "%e,%e,%e,%e\n", t, KE, PE, KE+PE);
    }
    fclose(data_ptr);
    fclose(energy_ptr);
    printf("Simulation complete. Data saved to LJ_data.csv and energy saved to LJ_energy.csv\n");
    return 0;
}