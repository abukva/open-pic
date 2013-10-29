#include <stdio.h>
#include <stdlib.h>
#include "lib.h"

void field_value_print(grid* grid_all, int counter, int GRID_X, int GRID_Y, int ex, int ey, int ez, int bx, int by, int bz, int charge)
{
	/* EX */
	if(ex)
	{
		moving_frame_print(grid_all->grid_ex, "out/ex/ex_%d.data", counter);
	}
	/* EY */
	if(ey)
	{
		moving_frame_print(grid_all->grid_ey, "out/ey/ey_%d.data", counter);
	}
	/* EZ */
	if(ez)
	{
		moving_frame_print(grid_all->grid_ez, "out/ez/ez_%d.data", counter);		
	}
	/* BX */
	if(bx)
	{
		moving_frame_print(grid_all->grid_bx, "out/bx/bx_%d.data", counter);
	}
	/* BY */
	if(by)
	{
		moving_frame_print(grid_all->grid_by, "out/by/by_%d.data", counter);
	}
	/* BZ */
	if(bz)
	{
		moving_frame_print(grid_all->grid_bz, "out/bz/bz_%d.data", counter);
	}
	/* CHARGE */
	if(charge)
	{
		char file_name[30];
		int i, j;
		sprintf(file_name, "out/charge/charge_%d.data", counter);
		FILE *fp;
		fp=fopen(file_name, "w");
		for(j=0;j<GRID_Y;j++)
		{
			for(i=0;i<GRID_X;i++)
			{
				fprintf(fp, "%3.8f ", grid_all->charge->elements[i*GRID_Y+j]);
			}
			fprintf(fp, "\n");
		}
		fclose(fp);
	}
	

}
void deposit_charge(particle* particle, grid* grid_all, int GRID_X, int GRID_Y, int number_of_particles, double scale_x, double scale_y)
{
	int i, j, x, y;

	for(i=0;i<GRID_X;i++)
	{
		for(j=0;j<GRID_Y;j++)
		{
			grid_all->charge->elements[i*GRID_Y+j]=0;
		}
	}
	for(i=0;i<number_of_particles;i++)
	{
		x=(int)(particle[i].x*scale_x);
		y=(int)(particle[i].y*scale_y);
		x=x-(grid_all->grid_bz->move);
		(grid_all->charge->elements[x*GRID_Y+y])+=particle[i].q_weight;
	}
}
void particle_value_print(particle* particle, int counter, int number_of_particles, int rank)
{
	int i;
	char file_name[30];
	sprintf(file_name, "out/particles/particles_%d_%d.data", counter, rank);
	FILE *fp;
	fp=fopen(file_name, "a+");
	for(i=0;i<number_of_particles;i++)
	{
		if(particle[i].px>=10.0)
		{
			fprintf(fp, "%.8f %.8f %.8f %.8f %.8f %.8f\n", particle[i].x, particle[i].y, particle[i].px, particle[i].py, particle[i].pz, particle[i].q_weight);
		}
	}
	fclose(fp);
}
