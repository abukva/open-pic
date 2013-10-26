#include <stdio.h>
#include <stdlib.h>
#include "lib.h"

field *particle_field(double x, double y, grid *grid_all, double scale_x, double scale_y, int GRID_X, int GRID_Y) 
{
	/*field for the current particle */

	field *curr=(field*)malloc(sizeof(field));
	
	double delta_x, delta_y, e_x, e_y, e_z, b_x, b_y, b_z, dx_help, dy_help;
	int x_abs, y_abs;

	x_abs=(int)(x*scale_x);
	y_abs=(int)(y*scale_y);

	delta_x=x*scale_x-x_abs-0.5;
	delta_y=y*scale_y-y_abs-0.5;

	/*interpolating fields at the particle location */

/*****************************************/

	if(delta_y>0)
	{
		dy_help=-(0.5-delta_y);
		e_x=((moving_frame_get(grid_all->grid_ex, x_abs, y_abs))*(0.5-delta_x)*(0.5-dy_help)+(moving_frame_get(grid_all->grid_ex, x_abs+1, y_abs))*(0.5+delta_x)*(0.5-dy_help)+(moving_frame_get(grid_all->grid_ex, x_abs, y_abs+1))*(0.5-delta_x)*(0.5+dy_help)+(moving_frame_get(grid_all->grid_ex, x_abs+1, y_abs+1))*(0.5+delta_x)*(0.5+dy_help));
	}
	else
	{
		dy_help=(0.5+delta_y);
		e_x=(moving_frame_get(grid_all->grid_ex, x_abs, y_abs-1)*(0.5-delta_x)*(0.5-dy_help)+moving_frame_get(grid_all->grid_ex, x_abs+1, y_abs-1)*(0.5+delta_x)*(0.5-dy_help)+moving_frame_get(grid_all->grid_ex, x_abs, y_abs)*(0.5-delta_x)*(0.5+dy_help)+moving_frame_get(grid_all->grid_ex, x_abs+1, y_abs)*(0.5+delta_x)*(0.5+dy_help));
	}	

/*****************************************/

	if(delta_x>0)
	{
		dx_help=-(0.5-delta_x);
		e_y=((moving_frame_get(grid_all->grid_ey, x_abs, y_abs))*(0.5-dx_help)*(0.5-delta_y)+(moving_frame_get(grid_all->grid_ey, x_abs+1, y_abs))*(0.5+dx_help)*(0.5-delta_y)+(moving_frame_get(grid_all->grid_ey, x_abs, y_abs+1))*(0.5-dx_help)*(0.5+delta_y)+(moving_frame_get(grid_all->grid_ey, x_abs+1, y_abs+1))*(0.5+dx_help)*(0.5+delta_y));
	}
	else
	{
		dx_help=(0.5+delta_x);
		e_y=((moving_frame_get(grid_all->grid_ey, x_abs-1, y_abs))*(0.5-dx_help)*(0.5-delta_y)+(moving_frame_get(grid_all->grid_ey, x_abs, y_abs))*(0.5+dx_help)*(0.5-delta_y)+(moving_frame_get(grid_all->grid_ey, x_abs-1, y_abs+1))*(0.5-dx_help)*(0.5+delta_y)+(moving_frame_get(grid_all->grid_ey, x_abs, y_abs+1))*(0.5+dx_help)*(0.5+delta_y));
	}	

/*****************************************/

	if(delta_x>0)
	{
		if(delta_y>0)
		{
			dx_help=-(0.5-delta_x);
			dy_help=-(0.5-delta_y);
			e_z=((moving_frame_get(grid_all->grid_ez, x_abs, y_abs))*(0.5-dx_help)*(0.5-dy_help)+(moving_frame_get(grid_all->grid_ez, x_abs+1, y_abs))*(0.5+dx_help)*(0.5-dy_help)+(moving_frame_get(grid_all->grid_ez, x_abs, y_abs+1))*(0.5-dx_help)*(0.5+dy_help)+(moving_frame_get(grid_all->grid_ez, x_abs+1, y_abs+1))*(0.5+dx_help)*(0.5+dy_help));
		}
		else
		{
			dx_help=-(0.5-delta_x);
			dy_help=0.5+delta_y;
			e_z=((moving_frame_get(grid_all->grid_ez, x_abs, y_abs-1))*(0.5-dx_help)*(0.5-dy_help)+(moving_frame_get(grid_all->grid_ez, x_abs+1, y_abs-1))*(0.5+dx_help)*(0.5-dy_help)+(moving_frame_get(grid_all->grid_ez, x_abs, y_abs))*(0.5-dx_help)*(0.5+dy_help)+(moving_frame_get(grid_all->grid_ez, x_abs+1, y_abs))*(0.5+dx_help)*(0.5+dy_help));
		}
	}
	else
	{
		if(delta_y>0)
		{
			dx_help=0.5+delta_x;
			dy_help=-(0.5-delta_y);
			e_z=((moving_frame_get(grid_all->grid_ez, x_abs-1, y_abs))*(0.5-dx_help)*(0.5-dy_help)+(moving_frame_get(grid_all->grid_ez, x_abs, y_abs))*(0.5+dx_help)*(0.5-dy_help)+(moving_frame_get(grid_all->grid_ez, x_abs-1, y_abs+1))*(0.5-dx_help)*(0.5+dy_help)+(moving_frame_get(grid_all->grid_ez, x_abs, y_abs+1))*(0.5+dx_help)*(0.5+dy_help));
		}
		else
		{
			dx_help=0.5+delta_x;
			dy_help=0.5+delta_y;
			e_z=((moving_frame_get(grid_all->grid_ez, x_abs-1, y_abs-1))*(0.5-dx_help)*(0.5-dy_help)+(moving_frame_get(grid_all->grid_ez, x_abs, y_abs-1))*(0.5+dx_help)*(0.5-dy_help)+(moving_frame_get(grid_all->grid_ez, x_abs-1, y_abs))*(0.5-dx_help)*(0.5+dy_help)+(moving_frame_get(grid_all->grid_ez, x_abs, y_abs))*(0.5+dx_help)*(0.5+dy_help));
		}
	}

/*****************************************/

	if(delta_x>0)
	{
		dx_help=-(0.5-delta_x);
		b_x=((moving_frame_get(grid_all->grid_bx, x_abs, y_abs))*(0.5-dx_help)*(0.5-delta_y)+(moving_frame_get(grid_all->grid_bx, x_abs+1, y_abs))*(0.5+dx_help)*(0.5-delta_y)+(moving_frame_get(grid_all->grid_bx, x_abs, y_abs+1))*(0.5-dx_help)*(0.5+delta_y)+(moving_frame_get(grid_all->grid_bx, x_abs+1, y_abs+1))*(0.5+dx_help)*(0.5+delta_y));
	}
	else
	{
		dx_help=(0.5+delta_x);
		b_x=((moving_frame_get(grid_all->grid_bx, x_abs-1, y_abs))*(0.5-dx_help)*(0.5-delta_y)+(moving_frame_get(grid_all->grid_bx, x_abs, y_abs))*(0.5+dx_help)*(0.5-delta_y)+(moving_frame_get(grid_all->grid_bx, x_abs-1, y_abs+1))*(0.5-dx_help)*(0.5+delta_y)+(moving_frame_get(grid_all->grid_bx, x_abs, y_abs+1))*(0.5+dx_help)*(0.5+delta_y));
	}	

/*****************************************/

	if(delta_y>0)
	{
		dy_help=-(0.5-delta_y);
		b_y=((moving_frame_get(grid_all->grid_by, x_abs, y_abs))*(0.5-delta_x)*(0.5-dy_help)+(moving_frame_get(grid_all->grid_by, x_abs+1, y_abs))*(0.5+delta_x)*(0.5-dy_help)+(moving_frame_get(grid_all->grid_by, x_abs, y_abs+1))*(0.5-delta_x)*(0.5+dy_help)+(moving_frame_get(grid_all->grid_by, x_abs+1, y_abs+1))*(0.5+delta_x)*(0.5+dy_help));
	}
	else
	{
		dy_help=(0.5+delta_y);
		b_y=((moving_frame_get(grid_all->grid_by, x_abs, y_abs-1))*(0.5-delta_x)*(0.5-dy_help)+(moving_frame_get(grid_all->grid_by, x_abs+1, y_abs-1))*(0.5+delta_x)*(0.5-dy_help)+(moving_frame_get(grid_all->grid_by, x_abs, y_abs))*(0.5-delta_x)*(0.5+dy_help)+(moving_frame_get(grid_all->grid_by, x_abs+1, y_abs))*(0.5+delta_x)*(0.5+dy_help));
	}	

/*****************************************/

	b_z=(moving_frame_get(grid_all->grid_bz, x_abs, y_abs)*(0.5-delta_x)*(0.5-delta_y)+moving_frame_get(grid_all->grid_bz, x_abs+1, y_abs)*(0.5+delta_x)*(0.5-delta_y)+moving_frame_get(grid_all->grid_bz, x_abs, y_abs+1)*(0.5-delta_x)*(0.5+delta_y)+moving_frame_get(grid_all->grid_bz, x_abs+1, y_abs+1)*(0.5+delta_x)*(0.5+delta_y));

/*****************************************/	

	curr->ex=e_x;
	curr->ey=e_y;
	curr->ez=e_z;
	curr->bx=b_x;
	curr->by=b_y;
	curr->bz=b_z;

	return curr;
}
