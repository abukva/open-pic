#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "lib.h"
#include "mpi.h"

int gcd(int a, int b)
{
	if(b==0)
		return a;
	else
		return gcd(b, a%b);
}

void move_frame(particle *particle, grid *grid_all, int* number_of_particles, int horizontal, int vertical, int GRID_X, int GRID_Y, double MAX_X, double MAX_Y, double density)
{
	int j, rank, size;
	int  move=moving_frame_get_move(grid_all->grid_bz);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	/* setting the new most right raw as frash new fields with all zeros */

	moving_frame_move(grid_all->grid_ex);
	moving_frame_move(grid_all->grid_ey);
	moving_frame_move(grid_all->grid_ez);
	moving_frame_move(grid_all->grid_bx);
	moving_frame_move(grid_all->grid_by);
	moving_frame_move(grid_all->grid_bz);
	moving_frame_move(grid_all->grid_jx);
	moving_frame_move(grid_all->grid_jy);
	moving_frame_move(grid_all->grid_jz);

	/*initialize new particles at the most right*/
	
	int per_process=GRID_Y/size;
	int rest=GRID_Y-(size)*per_process;
	if(rank==0)
	{
		per_process+=rest;
	}
	int new_number_of_particles=per_process*horizontal*vertical;	
	double scale_x=((double)(MAX_X)/(GRID_X));
	double scale_y=((double)(MAX_Y)/(GRID_Y));
	double start_x=(scale_x/(horizontal*2));
	start_x+=((GRID_X)+(move))*(scale_x);
	int start_particle_y;
	start_particle_y=per_process*(rank)+rest;
	if(rank==0)
		start_particle_y=0;
	double start_y=scale_y/(vertical*2)+(start_particle_y*scale_y);
	int per_cell=horizontal*vertical;
	for(j=(*number_of_particles);j<(new_number_of_particles+(*number_of_particles));j++)
	{
		particle[j].x=start_x;
		particle[j].y=start_y;

		start_x+=scale_x/(horizontal);

		if(start_x>=((GRID_X)+(move)+1)*(scale_x))
		{
			start_x=((GRID_X)+(move))*(scale_x)+(scale_x/(horizontal*2));
			start_y+=scale_y/(vertical);
		}

		/*set initial value for momentum */
		
		particle[j].py=0.0;
		particle[j].pz=0.0;
		particle[j].px=0.0;
		particle[j].q=-1.0;

		get_weight(particle, scale_x, scale_y, per_cell, &j, &new_number_of_particles, density);
	}
	*number_of_particles+=new_number_of_particles;
}
