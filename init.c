#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "lib.h"
#include "mpi.h"
 
void get_weight(particle *particle, double scale_x, double scale_y, int per_cell, int *i, int* number_of_particles, double density)
{
	double x=((int)(particle[*i].x*(1/scale_x)))*(scale_x);
	x+=scale_x/2;
	double y=((int)(particle[*i].y*(1/scale_y)))*(scale_y);
	y+=scale_y/2;
	
	/*if(x>100 && x<120)
	{
		particle[*i].q_weight=density*(x/10.0 - 10.0)/(per_cell);
	}
	else if(x>120 && x<140)
	{
		particle[*i].q_weight=density*(-x/20.0 + 8.0)/per_cell;
	}
	else if(x>140)
	{
		particle[*i].q_weight=(density*1.0)/per_cell;
	}
	else
	{
		(*i)--;
		(*number_of_particles)--;
	}*/
	if(x>20 && x<60)
	{
		if(y>9 && y<11)
		{
			particle[*i].q_weight=1.0/per_cell;
		}
		else
		{
			(*i)--;
			(*number_of_particles)--;
		}
	}	
	else
	{
		(*i)--;
		(*number_of_particles)--;
	}
	
	
	
	
	
	
}

/*initialize n partcles as structure with x and y components */

particle *init_particles(int horizontal, int vertical, int GRID_X, int GRID_Y, double MAX_X, double MAX_Y, int* number_of_particles, int per_cell, int start_particle_x, int end_particle_x, double density) 
{
	int help=2*(GRID_X*GRID_Y*horizontal*vertical);	
	particle *particle_s=(particle*)malloc(help*sizeof(particle));
	int i;
	double scale_x=((double)(MAX_X)/(GRID_X));
	double scale_y=((double)MAX_Y/GRID_Y);
	double start_x=(scale_x/(horizontal*2))+(start_particle_x*(scale_x));
	double start_y=scale_y/(vertical*2);
	
	for(i=0;i<*number_of_particles;i++) 
	{
		particle_s[i].x=start_x;

		particle_s[i].y=start_y;

		start_x+=scale_x/(horizontal);

		if(start_x>=(end_particle_x*(scale_x)))
		{
			start_x=(scale_x/(horizontal*2))+(start_particle_x*(scale_x));
			start_y+=scale_y/(vertical);
		}

		/*set initial value for momentum */
		
		particle_s[i].px=1000.0;
		particle_s[i].py=0.0;
		particle_s[i].pz=0.0;

		/* set value for "weight" */		

		particle_s[i].q=-1.0;

		get_weight(particle_s,scale_x,scale_y, per_cell, &i, number_of_particles, density);
	}

	return particle_s;
}

/* function that return the initial value for field */

double initial_field_value(double x, double y, double A0, double omega0, double tfwhm, double w0, double xc, double yc)
{
	double r=y-yc;
	double z=x-xc;
	
	if(z>3*w0)
		return 0;
	else if(z<(-3*w0))
		return 0;
		
	double E0=A0*omega0;
	double k=omega0;
	double zr=(omega0*w0*w0)/2;
	double w=w0*sqrt(1+(z/zr)*(z/zr));
	double R=z*(1+(zr/z)*(zr/z));
	double gouy=atan(z/zr);
	double result=(E0*w0/w)*exp(-(r*r)/(w*w))*cos(-k*z-k*(r*r)/(2*R)+gouy);
	double a=(x-xc)/tfwhm;
	if(a<0)
		a=-a;
	a=1-a;
	double b=a;
	if(a<0)
		b=0;
	double envelope=10*(b*b*b)-15*(b*b*b*b)+6*(b*b*b*b*b);
	result*=envelope;
	return result;
}

/* initialize structure 'grid' with magnetic and electric fields fields */

grid *init_grid(int x, int y, double x_scale, double y_scale, double A0, double omega0, double tfwhm, double w0, double xc, double yc) 
{
	grid *grid_s=(grid*)malloc(sizeof(grid));
	int i,j;
	double dx=1.0/x_scale;
	double dy=1.0/y_scale;

	/* initialize each grid as a 2D matrix, or array of pointers to pointers of a arrays of doubles */

	grid_s->grid_ex=moving_frame_create(x,y);
	grid_s->grid_ey=moving_frame_create(x,y);
	grid_s->grid_ez=moving_frame_create(x,y);
	grid_s->grid_bx=moving_frame_create(x,y);
	grid_s->grid_by=moving_frame_create(x,y);
	grid_s->grid_bz=moving_frame_create(x,y);
	grid_s->grid_jx=moving_frame_create(x,y);
	grid_s->grid_jy=moving_frame_create(x,y);
	grid_s->grid_jz=moving_frame_create(x,y);
	grid_s->charge=moving_frame_create(x,y);

	/* settings initial values for the fields */

	double k=0;
	double n=0;

	for(i=0;i<x;i++)
	{		
		for(j=0;j<y;j++)
		{
				/*moving_frame_set(grid_s->grid_ex, i, j, 7*k+5*(n+dy/2));
				moving_frame_set(grid_s->grid_ey, i, j, 7*(k+dx/2)+5*n);
				moving_frame_set(grid_s->grid_ez, i, j, 7*(k+dx/2)+5*(n+dy/2));
				moving_frame_set(grid_s->grid_bx, i, j, 7*(k+dx/2)+5*n);
				moving_frame_set(grid_s->grid_by, i, j, 7*k+5*(n+dy/2));
				moving_frame_set(grid_s->grid_bz, i, j, 7*k+5*n);*/
				//moving_frame_set(grid_s->grid_ey, i, j, initial_field_value(k+dx/2, n, A0, omega0, tfwhm, w0, xc, yc));
				//moving_frame_set(grid_s->grid_bz, i, j, initial_field_value(k, n, A0, omega0, tfwhm, w0, xc, yc));
				n+=dy;
		}
		k+=dx;
		n=0;
	}
	return grid_s;
}
