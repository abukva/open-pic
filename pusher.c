#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "lib.h"

int number_of_itter(double t, double dt) {
	return (int)t/dt;
}

void boundary_conditions_particles(particle *particle, int GRID_X, int GRID_Y, int* number_of_particles, int *i, double scale_x, double scale_y, grid* grid_all, double x_start, double y_start, double dt)
{
	int move=moving_frame_get_move(grid_all->grid_bz);
	double MAX_X=(GRID_X/scale_x);
	MAX_X+=(move+1)/scale_x;
	double MIN_X=(move+1)/scale_x;
	double MAX_Y=(GRID_Y/scale_y);

	if((particle[*i].x>=MIN_X && particle[*i].x<=MAX_X) && (particle[*i].y<=MAX_Y && particle[*i].y>=0))
	{
		current_deposition (grid_all, x_start, y_start, particle[*i].x, particle[*i].y, particle[*i].q, dt, scale_x, scale_y, GRID_X, GRID_Y, particle[*i].pz, particle[*i].q_weight);
	}
	else
	{
		if(particle[*i].x>MAX_X)
		{
			particle[*i].x=particle[(*number_of_particles)-1].x;
			particle[*i].y=particle[(*number_of_particles)-1].y;
			particle[*i].px=particle[(*number_of_particles)-1].px;
			particle[*i].py=particle[(*number_of_particles)-1].py;
			particle[*i].pz=particle[(*number_of_particles)-1].pz;
			particle[*i].q=particle[(*number_of_particles)-1].q;
			particle[*i].q_weight=particle[(*number_of_particles)-1].q_weight;
			(*number_of_particles)--;
			(*i)--;
		}
		else if(particle[*i].x<MIN_X)
		{
			particle[*i].x=particle[(*number_of_particles)-1].x;
			particle[*i].y=particle[(*number_of_particles)-1].y;
			particle[*i].px=particle[(*number_of_particles)-1].px;
			particle[*i].py=particle[(*number_of_particles)-1].py;
			particle[*i].pz=particle[(*number_of_particles)-1].pz;
			particle[*i].q=particle[(*number_of_particles)-1].q;
			particle[*i].q_weight=particle[(*number_of_particles)-1].q_weight;
			(*number_of_particles)--;
			(*i)--;
		}
		else if(particle[*i].y>MAX_Y)
		{
			//particle[*i].y=particle[*i].y-MAX_Y;
			particle[*i].x=particle[(*number_of_particles)-1].x;
			particle[*i].y=particle[(*number_of_particles)-1].y;
			particle[*i].px=particle[(*number_of_particles)-1].px;
			particle[*i].py=particle[(*number_of_particles)-1].py;
			particle[*i].pz=particle[(*number_of_particles)-1].pz;
			particle[*i].q=particle[(*number_of_particles)-1].q;
			particle[*i].q_weight=particle[(*number_of_particles)-1].q_weight;
			(*number_of_particles)--;
			(*i)--;
				
		}
		else if(particle[*i].y<0)
		{
			//particle[*i].y=MAX_Y+particle[*i].y;
			particle[*i].x=particle[(*number_of_particles)-1].x;
			particle[*i].y=particle[(*number_of_particles)-1].y;
			particle[*i].px=particle[(*number_of_particles)-1].px;
			particle[*i].py=particle[(*number_of_particles)-1].py;
			particle[*i].pz=particle[(*number_of_particles)-1].pz;
			particle[*i].q=particle[(*number_of_particles)-1].q;
			particle[*i].q_weight=particle[(*number_of_particles)-1].q_weight;
			(*number_of_particles)--;
			(*i)--;
		}
	}
}

void push_one_set(particle *particle, double dt, int* number_of_particles, grid *grid_all, double scale_x, double scale_y, int GRID_X, int GRID_Y) {
	int i;
	double gamma; /* gamma factor */
	double l,l_x,l_y,l_z,s_x,s_y,s_z; /* intermiediate B intervals */
	double px_prime,py_prime,pz_prime,px_help,py_help,pz_help; /*helper values */
	double x_start, y_start;
	for(i=0;i<*number_of_particles;i++) 
	{
		/*return magnetic and electric field at particles's current position */

		field* curr_field=particle_field(particle[i].x, particle[i].y, grid_all, scale_x, scale_y, GRID_X, GRID_Y);

		/*first set of equations for all three cooridnates (half tranistion of magnetic and electric field) */

		px_prime=particle[i].px+particle[i].q*(dt/2.0)*curr_field->ex; 
		py_prime=particle[i].py+particle[i].q*(dt/2.0)*curr_field->ey; 
		pz_prime=particle[i].pz+particle[i].q*(dt/2.0)*curr_field->ez;

		/*calculating rel. gamma factor */

		gamma=sqrt(1.0+(px_prime*px_prime+py_prime*py_prime+pz_prime*pz_prime));

		/* magnetic rotation */

		l_x=particle[i].q*(dt/2.0)*(curr_field->bx)/gamma;
		l_y=particle[i].q*(dt/2.0)*(curr_field->by)/gamma;
		l_z=particle[i].q*(dt/2.0)*(curr_field->bz)/gamma;

		particle[i].px=px_prime+(py_prime*l_z-pz_prime*l_y);
		particle[i].py=py_prime+(pz_prime*l_x-px_prime*l_z);
		particle[i].pz=pz_prime+(px_prime*l_y-py_prime*l_x);

		l=l_x*l_x+l_y*l_y+l_z*l_z;

		/* full magnetic rotation */

		s_x=2.0*l_x/(1.0+l);
		s_y=2.0*l_y/(1.0+l);
		s_z=2.0*l_z/(1.0+l);

		px_help=particle[i].px;
		py_help=particle[i].py;
		pz_help=particle[i].pz;
		
		particle[i].px=px_prime+(py_help*s_z-pz_help*s_y);
		particle[i].py=py_prime+(pz_help*s_x-px_help*s_z);
		particle[i].pz=pz_prime+(px_help*s_y-py_help*s_x);

		/*adding the other half of electric field */

		particle[i].px=particle[i].px+particle[i].q*(dt/2.0)*curr_field->ex;
		particle[i].py=particle[i].py+particle[i].q*(dt/2.0)*curr_field->ey;
		particle[i].pz=particle[i].pz+particle[i].q*(dt/2.0)*curr_field->ez;

		/*updating rel. gamma factor */

		gamma=sqrt(1.0+(particle[i].px*particle[i].px+particle[i].py*particle[i].py+particle[i].pz*particle[i].pz));

		/*updating coordinates */

		x_start=particle[i].x;
		y_start=particle[i].y;

		particle[i].x+=(particle[i].px/gamma)*dt;
		particle[i].y+=(particle[i].py/gamma)*dt;

		/* boundary conditions for particle */

		boundary_conditions_particles(particle, GRID_X, GRID_Y, number_of_particles, &i, scale_x, scale_y, grid_all, x_start, y_start, dt);	
		
		free(curr_field);
	}
}
