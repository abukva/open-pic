#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "lib.h"

#define LIMIT 0.00000000001

double abs_double(double value)
{
	if(value<0)
		return -value;
	else
		return value;
}

void solve_fields(grid *grid_all, double scale_x, double scale_y, double dt, int GRID_X, int GRID_Y)
{
	double dx=1.0/scale_x;
	double dy=1.0/scale_y;

	int i, j, i_help;
	double set=0;
	int move=moving_frame_get_move(grid_all->grid_bz);

	moving_frame_get_guard_cells(grid_all->grid_ex);
	moving_frame_get_guard_cells(grid_all->grid_ey);
	moving_frame_get_guard_cells(grid_all->grid_ez);
	moving_frame_get_guard_cells(grid_all->grid_bx);
	moving_frame_get_guard_cells(grid_all->grid_by);
	moving_frame_get_guard_cells(grid_all->grid_bz);

	for (i=move; i<move+GRID_X; i++)
	{
		for(j=0; j<GRID_Y; j++)
		{
			set=moving_frame_get(grid_all->grid_bx, i, j)-(moving_frame_get(grid_all->grid_ez, i, j)-moving_frame_get(grid_all->grid_ez, i, j-1))/dy*dt/2.0;
			moving_frame_set(grid_all->grid_bx, i, j, set);
			set=moving_frame_get(grid_all->grid_by, i, j)+(moving_frame_get(grid_all->grid_ez, i, j)-moving_frame_get(grid_all->grid_ez, i-1, j))/dx*dt/2.0;
			moving_frame_set(grid_all->grid_by, i, j, set);
			set=moving_frame_get(grid_all->grid_bz, i, j)-((moving_frame_get(grid_all->grid_ey, i, j)-moving_frame_get(grid_all->grid_ey, i-1, j))/dx-(moving_frame_get(grid_all->grid_ex, i, j)-moving_frame_get(grid_all->grid_ex, i, j-1))/dy)*dt/2.0;
			moving_frame_set(grid_all->grid_bz, i, j, set);
		}
	}
	for (i=move; i<move+GRID_X; i++)
	{
		for(j=0; j<GRID_Y; j++)
		{
			set=moving_frame_get(grid_all->grid_ex,i,j)+((moving_frame_get(grid_all->grid_bz, i, j+1)-moving_frame_get(grid_all->grid_bz, i, j))/dy - moving_frame_get(grid_all->grid_jx, i, j))*dt;
			moving_frame_set(grid_all->grid_ex, i, j, set);
			set=moving_frame_get(grid_all->grid_ey, i, j)+(-(moving_frame_get(grid_all->grid_bz, i+1, j)-moving_frame_get(grid_all->grid_bz, i, j))/dx - moving_frame_get(grid_all->grid_jy, i, j))*dt;
			moving_frame_set(grid_all->grid_ey, i, j, set);
			set=moving_frame_get(grid_all->grid_ez, i, j)+(-(moving_frame_get(grid_all->grid_bx, i, j+1)-moving_frame_get(grid_all->grid_bx, i, j))/dy + (moving_frame_get(grid_all->grid_by, i+1, j)-moving_frame_get(grid_all->grid_by, i, j))/dx - moving_frame_get(grid_all->grid_jz, i, j))*dt;
			moving_frame_set(grid_all->grid_ez, i, j, set);
		}
	}
	for (i=move; i<move+GRID_X; i++)
	{
		for(j=0; j<GRID_Y; j++)
		{
			set=moving_frame_get(grid_all->grid_bx, i, j)-(moving_frame_get(grid_all->grid_ez, i, j)-moving_frame_get(grid_all->grid_ez, i, j-1))/dy*dt/2.0;
			moving_frame_set(grid_all->grid_bx, i, j, set);
			set=moving_frame_get(grid_all->grid_by, i, j)+(moving_frame_get(grid_all->grid_ez, i, j)-moving_frame_get(grid_all->grid_ez, i-1, j))/dx*dt/2.0;
			moving_frame_set(grid_all->grid_by, i, j, set);
			set=moving_frame_get(grid_all->grid_bz, i, j)-((moving_frame_get(grid_all->grid_ey, i, j)-moving_frame_get(grid_all->grid_ey, i-1, j))/dx-(moving_frame_get(grid_all->grid_ex, i, j)-moving_frame_get(grid_all->grid_ex, i, j-1))/dy)*dt/2.0;
			moving_frame_set(grid_all->grid_bz, i, j, set);
		}
	}

	moving_frame_set_guard_cells(grid_all->grid_ex);
	moving_frame_set_guard_cells(grid_all->grid_ey);
	moving_frame_set_guard_cells(grid_all->grid_ez);
	moving_frame_set_guard_cells(grid_all->grid_bx);
	moving_frame_set_guard_cells(grid_all->grid_by);
	moving_frame_set_guard_cells(grid_all->grid_bz);

	moving_frame_zero(grid_all->grid_jx);
	moving_frame_zero(grid_all->grid_jy);
	moving_frame_zero(grid_all->grid_jz);
}

void get_coordinates(double x_start, double x_end, double y_start, double y_end, int* x_start_coordinate, int* x_end_coordinate, int* y_start_coordinate, int* y_end_coordinate, double scale_x, double scale_y)
{
	double delta1, delta2;

	/* x-coodinate */
	
	(*x_start_coordinate)=(int)(x_start*scale_x);
	(*x_end_coordinate)=(int)(x_end*scale_x);

/*	delta1=(x_start-((*x_start_coordinate)/scale_x))*scale_x;
	delta2=(x_end-((*x_end_coordinate)/scale_x))*scale_x;

	if(abs_double(0.5-delta1)<LIMIT && abs_double(0.5-delta2)<LIMIT)
	{
		if(x_start<x_end)
			(*x_start_coordinate)=(*x_end_coordinate);
		else
			(*x_end_coordinate)=(*x_start_coordinate);
	}
	else if(abs_double(0.5-delta1)<LIMIT)
	{
		if(x_start<x_end)
		{
			if(delta2>0.5)
				(*x_end_coordinate)++;
			(*x_start_coordinate)=(*x_end_coordinate);
		}
		else
		{		
			if(delta2>0.5)
				(*x_end_coordinate)++;
		}
	}
	else if(abs_double(0.5-delta2)<LIMIT)
	{
		if(x_start>x_end)
		{
			if(delta1>0.5)
				(*x_start_coordinate)++;
			(*x_end_coordinate)=(*x_start_coordinate);
		}
		else
		{
			if(delta1>0.5)
				(*x_start_coordinate)++;
		}
	}
	else
	{
		if(delta1>0.5)
			(*x_start_coordinate)++;
		if(delta2>0.5)
			(*x_end_coordinate)++;
	}
*/
	/* y-coordinate */

	(*y_start_coordinate)=(int)(y_start*scale_y);
	(*y_end_coordinate)=(int)(y_end*scale_y);

/*	delta1=(y_start-((*y_start_coordinate)/scale_y))*scale_y;
	delta2=(y_end-((*y_end_coordinate)/scale_y))*scale_y;

	if(abs_double(0.5-delta1)<LIMIT && abs_double(0.5-delta2)<LIMIT)
	{
		if(y_start<y_end)
			(*y_start_coordinate)=(*y_end_coordinate);
		else
			(*y_end_coordinate)=(*y_start_coordinate);
	}
	else if(abs_double(0.5-delta1)<LIMIT)
	{
		if(y_start<y_end)
		{
			if(delta2>0.5)
				(*y_end_coordinate)++;
			(*y_start_coordinate)=(*y_end_coordinate);
		}
		else
		{
			if(delta2>0.5)
				(*y_end_coordinate)++;
		}
	}
	else if(abs_double(0.5-delta2)<LIMIT)
	{
		if(y_start>y_end)
		{
			if(delta1>0.5)
				(*y_start_coordinate)++;
			(*y_end_coordinate)=(*y_start_coordinate);
		}
		else
		{
			if(delta1>0.5)
				(*y_start_coordinate)++;
		}
	}
	else
	{
		if(delta1>0.5)
			(*y_start_coordinate)++;
		if(delta2>0.5)
			(*y_end_coordinate)++;
	}*/
}

void get_splits(split* split_data, int *splits, double x_start, double y_start, double x_end, double y_end, double scale_x, double scale_y)
{
	double delta_x;
	double delta_y;
	int x_start_coordinate, x_end_coordinate, y_start_coordinate, y_end_coordinate;

	get_coordinates(x_start, x_end, y_start, y_end, &x_start_coordinate, &x_end_coordinate, &y_start_coordinate, &y_end_coordinate, scale_x, scale_y);

	(*splits)=abs(x_start_coordinate-x_end_coordinate)+abs(y_start_coordinate-y_end_coordinate);

	if(*splits==0)
	{
		split_data->x0[0]=x_start*scale_x-x_start_coordinate-0.5;
		split_data->x1[0]=x_end*scale_x-x_end_coordinate-0.5;
		split_data->y0[0]=y_start*scale_y-y_start_coordinate-0.5;
		split_data->y1[0]=y_end*scale_y-y_end_coordinate-0.5;
		split_data->i_cell[0]=x_start_coordinate;
		split_data->j_cell[0]=y_start_coordinate;
	}

	if(*splits==1)
	{
		split_data->x0[0]=x_start*scale_x-x_start_coordinate-0.5;
		split_data->y0[0]=y_start*scale_y-y_start_coordinate-0.5;
		split_data->i_cell[0]=x_start_coordinate;
		split_data->j_cell[0]=y_start_coordinate;
		if(x_start_coordinate-x_end_coordinate!=0)
		{
			if(x_end_coordinate>x_start_coordinate)
				split_data->x1[0]=0.5;
			else
				split_data->x1[0]=-0.5;

			delta_x=(split_data->x1[0]-split_data->x0[0]);
			split_data->y1[0]=split_data->y0[0]+((delta_x*((y_end-y_start))/((x_end-x_start))));
			split_data->x0[1]=-split_data->x1[0];
			split_data->x1[1]=x_end*scale_x-x_end_coordinate-0.5;
			split_data->y0[1]=split_data->y1[0];
			split_data->y1[1]=y_end*scale_y-y_end_coordinate-0.5;
			if(x_end_coordinate>x_start_coordinate)
				split_data->i_cell[1]=split_data->i_cell[0]+1;
			else
				split_data->i_cell[1]=split_data->i_cell[0]-1;
			split_data->j_cell[1]=split_data->j_cell[0];
		}
		else
		{
			if(y_end_coordinate>y_start_coordinate)
				split_data->y1[0]=0.5;
			else
				split_data->y1[0]=-0.5;
			delta_y=(split_data->y1[0]-split_data->y0[0]);
			split_data->x1[0]=split_data->x0[0]+(delta_y*((x_end-x_start))/((y_end-y_start)));
			split_data->y0[1]=-split_data->y1[0];
			split_data->y1[1]=y_end*scale_y-y_end_coordinate-0.5;
			split_data->x0[1]=split_data->x1[0];
			split_data->x1[1]=x_end*scale_x-x_end_coordinate-0.5;
			split_data->i_cell[1]=split_data->i_cell[0];
			if(y_end_coordinate>y_start_coordinate)
				split_data->j_cell[1]=split_data->j_cell[0]+1;
			else
				split_data->j_cell[1]=split_data->j_cell[0]-1;
		}
	}

	if(*splits==2)
	{
		split_data->x0[0]=x_start*scale_x-x_start_coordinate-0.5;
		split_data->y0[0]=y_start*scale_y-y_start_coordinate-0.5;
		split_data->i_cell[0]=x_start_coordinate;
		split_data->j_cell[0]=y_start_coordinate;
		double x_start_help, y_start_help;
		x_start_help=split_data->x0[0];
		y_start_help=split_data->y0[0];
		if(split_data->x0[0]<0)
			x_start_help=-split_data->x0[0];
		if(split_data->y0[0]<0)
			y_start_help=-split_data->y0[0];
		double x_positive_delta, y_positive_delta;
		if(x_end-x_start<0)
			x_positive_delta=(x_start-x_end);
		else
			x_positive_delta=(x_end-x_start);
		if(y_end-y_start<0)
			y_positive_delta=(y_start-y_end);
		else
			y_positive_delta=(y_end-y_start);

		if(((0.5-y_start_help)/(0.5-x_start_help))>((y_positive_delta)/(x_positive_delta)))
		{
			if(x_end_coordinate>x_start_coordinate)
				split_data->x1[0]=0.5;
			else
				split_data->x1[0]=-0.5;

			delta_x=(split_data->x1[0]-split_data->x0[0]);
			split_data->y1[0]=split_data->y0[0]+(delta_x*((y_end-y_start))/((x_end-x_start)));
			split_data->x0[1]=-split_data->x1[0];
			split_data->y0[1]=split_data->y1[0];
			if(y_end_coordinate>y_start_coordinate)
				split_data->y1[1]=0.5;
			else
				split_data->y1[1]=-0.5;
			delta_y=(split_data->y1[1]-split_data->y0[1]);
			split_data->x1[1]=split_data->x0[1]+(delta_y*((x_end-x_start))/((y_end-y_start)));
			if(x_end_coordinate>x_start_coordinate)
				split_data->i_cell[1]=split_data->i_cell[0]+1;
			else
				split_data->i_cell[1]=split_data->i_cell[0]-1;
			split_data->j_cell[1]=split_data->j_cell[0];
			split_data->x0[2]=split_data->x1[1];
			split_data->y0[2]=-split_data->y1[1];
			split_data->x1[2]=x_end*scale_x-x_end_coordinate-0.5;
			split_data->y1[2]=y_end*scale_y-y_end_coordinate-0.5;
			split_data->i_cell[2]=split_data->i_cell[1];
			if(y_end_coordinate>y_start_coordinate)
				split_data->j_cell[2]=split_data->j_cell[1]+1;
			else
				split_data->j_cell[2]=split_data->j_cell[1]-1;
		}
		else
		{
			if(y_end_coordinate>y_start_coordinate)
				split_data->y1[0]=0.5;
			else
				split_data->y1[0]=-0.5;
			delta_y=(split_data->y1[0]-split_data->y0[0]);
			split_data->x1[0]=split_data->x0[0]+(delta_y*((x_end-x_start))/((y_end-y_start)));
			split_data->y0[1]=-split_data->y1[0];
			split_data->x0[1]=split_data->x1[0];
			if(x_end_coordinate-x_start_coordinate)
				split_data->x1[1]=0.5;
			else
				split_data->x1[1]=-0.5;
			delta_x=(split_data->x1[1]-split_data->x0[1]);
			split_data->y1[1]=split_data->y0[1]+(delta_x*((y_end-y_start))/((x_end-x_start)));
			split_data->i_cell[1]=split_data->i_cell[0];
			if(y_end_coordinate>y_start_coordinate)
				split_data->j_cell[1]=split_data->j_cell[0]+1;
			else
				split_data->j_cell[1]=split_data->j_cell[0]-1;
			split_data->x0[2]=-split_data->x1[1];
			split_data->y0[2]=split_data->y1[1];
			split_data->x1[2]=x_end*scale_x-x_end_coordinate-0.5;
			split_data->y1[2]=y_end*scale_y-y_end_coordinate-0.5;
			if(x_end_coordinate>x_start_coordinate)
				split_data->i_cell[2]=split_data->i_cell[1]+1;
			else
				split_data->i_cell[2]=split_data->i_cell[1]-1;
			split_data->j_cell[2]=split_data->j_cell[1];			
		}
	}
}

void current_deposition (grid *grid_all, double x_start, double y_start, double x_end, double y_end, double q, double dt, double scale_x, double scale_y, int GRID_X, int GRID_Y, double pz, double q_weight)
{
	split* split_data=(split*)malloc(sizeof(split));
	int splits, i;
	splits=0;
	get_splits(split_data, &splits, x_start, y_start, x_end, y_end, scale_x, scale_y);
	double jnorm_x=1.0/(2.0*dt*scale_x);
	double jnorm_y=1.0/(2.0*dt*scale_y);
	double x0, x1, y0, y1, qnx, qny, qvz, S0x_0, S0x_1, S1x_0, S1x_1, S0y_0, S0y_1, S1y_0, S1y_1, wl1, wl2, wp1_0, wp1_1, wp2_0, wp2_1;
	int i_cell, j_cell;
	for(i=0;i<=splits;i++)
	{
		x0=split_data->x0[i];
		x1=split_data->x1[i];
		y0=split_data->y0[i];
		y1=split_data->y1[i];

		i_cell=split_data->i_cell[i];
		j_cell=split_data->j_cell[i];

		qnx=q*q_weight*jnorm_x;
		qny=q*q_weight*jnorm_y;
		qvz=(q*q_weight*pz)/3.0;

		S0x_0=0.5-x0;
		S0x_1=0.5+x0;

		S1x_0=0.5-x1;
		S1x_1=0.5+x1;

		S0y_0=0.5-y0;
		S0y_1=0.5+y0;

		S1y_0=0.5-y1;
		S1y_1=0.5+y1;

		wl1=qnx*(-x0+x1);
		wl2=qny*(-y0+y1);

		wp1_0=S0y_0+S1y_0;
		wp1_1=S0y_1+S1y_1;

		wp2_0=S0x_0+S1x_0;
		wp2_1=S0x_1+S1x_1;		

		/*
		moving_frame_set(grid_all->grid_jx, i_cell, j_cell, moving_frame_get(grid_all->grid_jx, i_cell, j_cell)+wl1*wp1_0);
		moving_frame_set(grid_all->grid_jx, i_cell, j_cell+1, moving_frame_get(grid_all->grid_jx, i_cell, j_cell+1)+wl1*wp1_1);
		
		moving_frame_set(grid_all->grid_jy, i_cell, j_cell, moving_frame_get(grid_all->grid_jy, i_cell, j_cell)+wl2*wp2_0);
		moving_frame_set(grid_all->grid_jy, i_cell+1, j_cell, moving_frame_get(grid_all->grid_jy, i_cell+1, j_cell)+wl2*wp2_1);

		moving_frame_set(grid_all->grid_jz, i_cell, j_cell, moving_frame_get(grid_all->grid_jz, i_cell, j_cell)+qvz*(S0x_0*S0y_0+S1x_0*S1y_0+S0x_0*S1y_0+S1x_0*S0y_0)/2);
		moving_frame_set(grid_all->grid_jz, i_cell+1, j_cell, moving_frame_get(grid_all->grid_jz, i_cell+1, j_cell)+qvz*(S0x_1*S0y_0+S1x_1*S1y_0+S0x_1*S1y_0+S1x_1*S0y_0)/2);
		moving_frame_set(grid_all->grid_jz, i_cell, j_cell+1, moving_frame_get(grid_all->grid_jz, i_cell, j_cell+1)+qvz*(S0x_0*S0y_1+S1x_0*S1y_1+S0x_0*S1y_1+S1x_0*S0y_1)/2);
		moving_frame_set(grid_all->grid_jz, i_cell+1, j_cell+1, moving_frame_get(grid_all->grid_jz, i_cell+1, j_cell+1)+qvz*(S0x_1*S0y_1+S1x_1*S1y_1+S0x_1*S1y_1+S1x_1*S0y_1)/2);

		*/		
	
		//ove nezakomentarisane sam izveo sa marijom a ove gore su iz rikardovog koda direktno prepisane

		

		moving_frame_set(grid_all->grid_jx, i_cell, j_cell-1, moving_frame_get(grid_all->grid_jx, i_cell, j_cell-1)+wl1*wp1_0);
		moving_frame_set(grid_all->grid_jx, i_cell, j_cell, moving_frame_get(grid_all->grid_jx, i_cell, j_cell)+wl1*wp1_1);
		
		moving_frame_set(grid_all->grid_jy, i_cell-1, j_cell, moving_frame_get(grid_all->grid_jy, i_cell-1, j_cell)+wl2*wp2_0);
		moving_frame_set(grid_all->grid_jy, i_cell, j_cell, moving_frame_get(grid_all->grid_jy, i_cell, j_cell)+wl2*wp2_1);

		moving_frame_set(grid_all->grid_jz, i_cell, j_cell, moving_frame_get(grid_all->grid_jz, i_cell, j_cell)+qvz*(S0x_0*S0y_0+S1x_0*S1y_0+S0x_0*S1y_0+S1x_0*S0y_0)/2.0);
		moving_frame_set(grid_all->grid_jz, i_cell+1, j_cell, moving_frame_get(grid_all->grid_jz, i_cell+1, j_cell)+qvz*(S0x_1*S0y_0+S1x_1*S1y_0+S0x_1*S1y_0+S1x_1*S0y_0)/2.0);
		moving_frame_set(grid_all->grid_jz, i_cell, j_cell+1, moving_frame_get(grid_all->grid_jz, i_cell, j_cell+1)+qvz*(S0x_0*S0y_1+S1x_0*S1y_1+S0x_0*S1y_1+S1x_0*S0y_1)/2.0);
		moving_frame_set(grid_all->grid_jz, i_cell+1, j_cell+1, moving_frame_get(grid_all->grid_jz, i_cell+1, j_cell+1)+qvz*(S0x_1*S0y_1+S1x_1*S1y_1+S0x_1*S1y_1+S1x_1*S0y_1)/2.0);

		//da proverim indekse za ovaj jz, treba sve da ide -1

		
	}
	free(split_data);
}
