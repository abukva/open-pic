#include <stdio.h>
#include <stdlib.h>
#include "lib.h"
#include <time.h>
#include <math.h>
#include "mpi.h"
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

int main(int argc, char *argv[])
{
	if(argc!=2)
	{
		printf("Wrong number of arguments");
		return 0;
	}
	clock_t begin, end;
	double time_spend;
	int rank, size;
	int per_process, rest;
	int number_of_particles;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if(rank==0)
	{		
		begin=clock();
	}	

	/*all variables*/	

	int GRID_X, GRID_Y; /* x and y dimensions of a grid */
	double MAX_X, MAX_Y; /*maximum x and y coordinates*/
	double t,dt; /* time and time step */
	int print_every; /* print data every n steps */	
	int frame_switch; /* on(1) off(0) switch for frame moving */
	int horizontal; /* number of particles on verical in one cell */
	int vertical; /* number of particles on horizontal in one cell */
	double density;
	double A0, omega0, tfwhm, w0, xc, yc; /*laser parameters*/
	int ex, ey, ez, bx, by, bz, jx, jy, jz, charge, particles_print; /* printing parameters */

	/* Read data */
	char buf[1000];
	if(rank==0)
		printf("Reading the data...\n");

	FILE* data_file = fopen(argv[1], "r");
	fscanf(data_file, "%d%d", &GRID_X, &GRID_Y); fgets(buf, 1000, data_file);
	fscanf(data_file, "%lf%lf", &MAX_X, &MAX_Y); fgets(buf, 1000, data_file);
	fscanf(data_file, "%lf%lf", &t, &dt); fgets(buf, 1000, data_file);
	fscanf(data_file, "%d", &print_every); fgets(buf, 1000, data_file);
	fscanf(data_file, "%d", &frame_switch); fgets(buf, 1000, data_file);
	fscanf(data_file, "%d%d", &horizontal, &vertical); fgets(buf, 1000, data_file);
	fscanf(data_file, "%lf", &density); fgets(buf, 1000, data_file);
	fscanf(data_file, "%lf", &A0); fgets(buf, 1000, data_file);
	fscanf(data_file, "%lf", &omega0); fgets(buf, 1000, data_file);
	fscanf(data_file, "%lf", &tfwhm); fgets(buf, 1000, data_file);
	fscanf(data_file, "%lf", &w0); fgets(buf, 1000, data_file);
	fscanf(data_file, "%lf", &xc); fgets(buf, 1000, data_file);
	fscanf(data_file, "%d%d%d%d%d%d%d%d%d%d%d", &ex, &ey, &ez, &bx, &by, &bz, &jx, &jy, &jz, &charge, &particles_print); fgets(buf, 1000, data_file);
	fclose(data_file);

	

	if(rank==0)
	{
		char folder_name[100]="";
		strcat(folder_name, argv[1]);
		strcat(folder_name, ".out");
		mkdir(folder_name, 0777);
		chdir(folder_name);
		mkdir("out", 0777);
		mkdir("out/ex", 0777);
		mkdir("out/ey", 0777);
		mkdir("out/ez", 0777);
		mkdir("out/bx", 0777);
		mkdir("out/by", 0777);
		mkdir("out/bz", 0777);
		mkdir("out/jx", 0777);
		mkdir("out/jy", 0777);
		mkdir("out/jz", 0777);
		mkdir("out/charge", 0777);
		mkdir("out/particles", 0777);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	if(rank!=0)
	{
		char folder_name[100]="";
		strcat(folder_name, argv[1]);
		strcat(folder_name, ".out");
		chdir(folder_name);
	}	

	double x_scale, y_scale;
	int per_cell=horizontal*vertical;
	x_scale=((double)GRID_X)/MAX_X;
	y_scale=((double)GRID_Y)/MAX_Y;
	
	yc=MAX_Y/2; /* y position of laser center */

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
	
	if(rank==0)
		printf("Simulation is starting...\n");

	/* Start simulation */

	simulation(particles,t,dt,&number_of_particles,grid_all,print_every, x_scale, y_scale, GRID_X, GRID_Y, frame_switch, horizontal, vertical, ex, ey, ez, bx, by, bz, charge, density, jx, jy, jz, particles_print);

	if(rank==0)
	{
		end=clock();
		time_spend=(double)(end-begin)/CLOCKS_PER_SEC;
		printf("\nTime of execution: %.1f sec, %.1f min %.2f h\n", time_spend, time_spend/60, time_spend/3600.0);
	}

	MPI_Finalize();
	return 0;
}
