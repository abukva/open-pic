#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "lib.h"
#include "mpi.h"

void simulation(particle* particle, double t, double dt, int* number_of_particles, grid *grid_all, int print_every, double scale_x, double scale_y, int GRID_X, int GRID_Y, int frame_switch, int horizontal, int vertical, int ex, int ey, int ez, int bx, int by, int bz, int charge, double density, int jx, int jy, int jz, int particles) 
{
	int number_of_itter_v=number_of_itter(t,dt);
	int i, m, n;
	double MAX_X=GRID_X/scale_x;
	double MAX_Y=GRID_Y/scale_y;
	int rank;
	int m_x, n_t;
	m_x=100000.0/scale_x;
	n_t=dt*100000;
	int gcd_v=gcd(m_x, n_t);
	m_x=m_x/gcd_v;
	n_t=n_t/gcd_v;
	int to_move=0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	double* recive=(double*)malloc((GRID_X*GRID_Y)*sizeof(double));
	for(i=0;i<=number_of_itter_v;i++) {
		if(i%print_every==0)
		{
			if(rank==0)
			{
				double percent=((double)i/number_of_itter_v)*100;
				printf("Simulation: %.2f%% complete.\n", percent);
			}
			if(charge)
			{
				deposit_charge(particle, grid_all, GRID_X, GRID_Y, *number_of_particles, scale_x, scale_y);

				MPI_Reduce(grid_all->charge->elements, recive, ((GRID_X)*GRID_Y), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
				if(rank==0)
				{
					for(m=0;m<GRID_X;m++)
					{
						for(n=0;n<GRID_Y;n++)
						{
							grid_all->charge->elements[(m)*(GRID_Y)+n]=recive[(m)*(GRID_Y)+n];
						}
					}
				}
			}

			push_one_set(particle, dt, number_of_particles, grid_all, scale_x, scale_y, GRID_X, GRID_Y);

			MPI_Barrier(MPI_COMM_WORLD);

			MPI_Allreduce(grid_all->grid_jx->elements, recive, ((GRID_X)*GRID_Y), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			for(m=0;m<GRID_X;m++)
			{
				for(n=0;n<GRID_Y;n++)
				{
					grid_all->grid_jx->elements[(m)*(GRID_Y)+n]=recive[(m)*(GRID_Y)+n];
				}
			}

			MPI_Allreduce(grid_all->grid_jy->elements, recive, ((GRID_X)*GRID_Y), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			for(m=0;m<GRID_X;m++)
			{
				for(n=0;n<GRID_Y;n++)
				{
					grid_all->grid_jy->elements[(m)*(GRID_Y)+n]=recive[(m)*(GRID_Y)+n];
				}
			}

			MPI_Allreduce(grid_all->grid_jz->elements, recive, ((GRID_X)*GRID_Y), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			for(m=0;m<GRID_X;m++)
			{
				for(n=0;n<GRID_Y;n++)
				{
					grid_all->grid_jz->elements[(m)*(GRID_Y)+n]=recive[(m)*(GRID_Y)+n];
				}
			}
			
			MPI_Barrier(MPI_COMM_WORLD);
			
			if(rank==0)
			{
				if(jx)
					moving_frame_print(grid_all->grid_jx, "out/jx/jx_%d.data", i);
				if(jy)
					moving_frame_print(grid_all->grid_jy, "out/jy/jy_%d.data", i);
				if(jz)
					moving_frame_print(grid_all->grid_jz, "out/jz/jz_%d.data", i);
			}
			
			solve_fields(grid_all, scale_x, scale_y, dt, GRID_X, GRID_Y);

			if(rank==0)
			{
				field_value_print(grid_all, i, GRID_X, GRID_Y, ex, ey, ez, bx, by, bz, charge);
			}
			if(particles)
				particle_value_print(particle, i, *number_of_particles, rank);

			MPI_Barrier(MPI_COMM_WORLD);

			if(frame_switch)
			{
				to_move+=n_t;
				if(to_move>m_x)
				{
					to_move-=m_x;
					move_frame(particle, grid_all, number_of_particles, horizontal, vertical, GRID_X, GRID_Y, MAX_X, MAX_Y, density);
				}
			}
		}
		else
		{
			
			push_one_set(particle, dt, number_of_particles, grid_all, scale_x, scale_y, GRID_X, GRID_Y);

			MPI_Barrier(MPI_COMM_WORLD);

			MPI_Allreduce(grid_all->grid_jx->elements, recive, ((GRID_X)*GRID_Y), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			for(m=0;m<GRID_X;m++)
			{
				for(n=0;n<GRID_Y;n++)
				{
					grid_all->grid_jx->elements[(m)*(GRID_Y)+n]=recive[(m)*(GRID_Y)+n];
				}
			}

			MPI_Allreduce(grid_all->grid_jy->elements, recive, ((GRID_X)*GRID_Y), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			for(m=0;m<GRID_X;m++)
			{
				for(n=0;n<GRID_Y;n++)
				{
					grid_all->grid_jy->elements[(m)*(GRID_Y)+n]=recive[(m)*(GRID_Y)+n];
				}
			}

			MPI_Allreduce(grid_all->grid_jz->elements, recive, ((GRID_X)*GRID_Y), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			for(m=0;m<GRID_X;m++)
			{
				for(n=0;n<GRID_Y;n++)
				{
					grid_all->grid_jz->elements[(m)*(GRID_Y)+n]=recive[(m)*(GRID_Y)+n];
				}
			}
			
			MPI_Barrier(MPI_COMM_WORLD);			
			
			solve_fields(grid_all, scale_x, scale_y, dt, GRID_X, GRID_Y);

			if(frame_switch)
			{
				to_move+=n_t;
				if(to_move>m_x)
				{
					to_move-=m_x;
					move_frame(particle, grid_all, number_of_particles, horizontal, vertical, GRID_X, GRID_Y, MAX_X, MAX_Y, density);
				}
			}
		}
		
	}
	free(recive);
}


