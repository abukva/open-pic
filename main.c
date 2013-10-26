#include <stdio.h>
#include <stdlib.h>
#include "lib.h"
#include <time.h>
#include <math.h>
#include "mpi.h"

#define GRID_X 3000 /* x dimension of a grid */
#define GRID_Y 160 /* y dimension of a grid */
#define print_every 120 /* print data every n steps */

int main(int argc, char *argv[])
{
	clock_t begin, end;
	double time_spend;
	int rank, size;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if(rank==0)
	{		
		begin=clock();
	}	

	double t,dt;
	int i,j, frame_switch;
	frame_switch=0; /* on(1) off(0) switch for frame moving */
	double MAX_X, MAX_Y;
	MAX_X=600; /* maximum x coordinate */
	MAX_Y=800; /* maximum y coordinate */
	int horizontal=2; /* number of particles on verical in one cell */
	int vertical=2; /* number of particles on horizontal in one cell */
	double density=0.0016;
	int number_of_particles;
	int per_process, rest;

	t=1200;
	dt=0.1;

	double x_scale, y_scale;
	int per_cell=horizontal*vertical;
	x_scale=((double)GRID_X)/MAX_X;
	y_scale=((double)GRID_Y)/MAX_Y;

	/*laser parameters*/
	double A0, omega0, tfwhm, w0, xc, yc;
	xc=480;
	yc=MAX_Y/2;
	A0=4;
	omega0=1;
	w0=50;
	tfwhm=100;

	/* init of particles and grids */
	if(rank==0)
		printf("Initialization of grid...\n");
	grid *grid_all=init_grid(GRID_X, GRID_Y, x_scale, y_scale, A0, omega0, tfwhm, w0, xc, yc);
	
	int start_particle_x, end_particle_x;
	per_process=GRID_X/size;
	rest=GRID_X-(size)*per_process;
	if(rank==0)
	{
		start_particle_x=0;
		end_particle_x=per_process+rest;
	}
	else
	{
		start_particle_x=(rank)*per_process+rest;
		end_particle_x=(rank+1)*per_process+rest;
	}
	number_of_particles=(end_particle_x-start_particle_x)*GRID_Y*per_cell;
	if(rank==0)
		printf("Initialization of particles...\n");
	particle *particles=init_particles(horizontal, vertical, GRID_X, GRID_Y, MAX_X, MAX_Y, &number_of_particles, per_cell, start_particle_x, end_particle_x, density);
	/* Printing options */
	int ex=1;
	int ey=1;
	int ez=1;
	int bx=1;
	int by=1;
	int bz=1;
	int charge=1;
	if(rank==0)
		printf("Simulation is starting...\n");
	/* Start simulation */

	simulation(particles,t,dt,&number_of_particles,grid_all,print_every, x_scale, y_scale, GRID_X, GRID_Y, frame_switch, horizontal, vertical, ex, ey, ez, bx, by, bz, charge, density);

	if(rank==0)
	{
		end=clock();
		time_spend=(double)(end-begin)/CLOCKS_PER_SEC;
		printf("\nTime of execution: %.1f sec, %.1f min %.2f h\n", time_spend, time_spend/60, time_spend/3600.0);
	}

	MPI_Finalize();
	return 0;
}
