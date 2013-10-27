#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "lib.h"

void moving_frame_zero(moving_frame* moving_frame)
{
	int i, j;
	for(i=0;i<moving_frame->GRID_X;i++)
	{		
		for(j=0;j<moving_frame->GRID_Y;j++)
		{
				moving_frame->elements[i*moving_frame->GRID_Y+j]=0;
		}
	}
	for(i=0;i<2;i++)
	{
		for(j=0;j<moving_frame->GRID_Y;j++)
		{
			moving_frame->guard_cells[i*moving_frame->GRID_Y+j]=0;
		}
	}
}

moving_frame* moving_frame_create(int GRID_X, int GRID_Y)
{
	int i,j;
	moving_frame* new_moving_frame=(moving_frame*)malloc(sizeof(moving_frame));
	new_moving_frame->GRID_X=GRID_X;
	new_moving_frame->GRID_Y=GRID_Y;
	new_moving_frame->elements=(double*)malloc((GRID_X*GRID_Y)*sizeof(double));
	new_moving_frame->guard_cells=(double*)malloc(4*GRID_Y*sizeof(double));
	new_moving_frame->move=0;
	moving_frame_zero(new_moving_frame);
	return new_moving_frame;
}

double moving_frame_get(moving_frame* moving_frame, int lab_x, int lab_y)
{
	if(lab_x<moving_frame->move || lab_x>=(moving_frame->move+moving_frame->GRID_X))
	{
		int y_position=lab_y;
		if(lab_y>moving_frame->GRID_Y-1)
			y_position-=moving_frame->GRID_Y;
		if(lab_y<0)
			y_position+=moving_frame->GRID_Y;
		if(lab_x<moving_frame->move)
			return moving_frame->guard_cells[lab_y];
		else
			return moving_frame->guard_cells[moving_frame->GRID_Y+lab_y];
	}
	else
	{
		int x_position=lab_x%moving_frame->GRID_X;
		int y_position=lab_y;
		if(lab_y>=moving_frame->GRID_Y)
			y_position-=moving_frame->GRID_Y;
		if(lab_y<0)
			y_position+=moving_frame->GRID_Y;		
		return moving_frame->elements[x_position*(moving_frame->GRID_Y)+y_position];
	}
}

void moving_frame_set(moving_frame* moving_frame, int lab_x, int lab_y, double value)
{
	if(lab_x<moving_frame->move || lab_x>=(moving_frame->move+moving_frame->GRID_X))
	{
		
	}
	else
	{
		int x_position=lab_x%(moving_frame->GRID_X);
		int y_position=lab_y;
		if(lab_y>=moving_frame->GRID_Y)
			y_position-=moving_frame->GRID_Y;
		if(lab_y<0)
			y_position+=moving_frame->GRID_Y;
		moving_frame->elements[x_position*(moving_frame->GRID_Y)+y_position]=value;
	}
}

int moving_frame_get_move(moving_frame* moving_frame)
{
	return moving_frame->move;
}

void moving_frame_move(moving_frame* moving_frame)
{
	int j;
	
	for(j=0;j<moving_frame->GRID_Y;j++)
	{
		moving_frame->elements[(((moving_frame->move)%(moving_frame->GRID_X)))*(moving_frame->GRID_Y)+j]=0;
	}
	(moving_frame->move)++;
}

void moving_frame_print(moving_frame* moving_frame, char* pattern, int counter)
{
	char filename[30];
	sprintf(filename, pattern, counter);
	FILE *fp;
	fp=fopen(filename, "w");
	int i, j;
	for(j=0;j<moving_frame->GRID_Y;j++)
	{
		for(i=moving_frame->move;i<(moving_frame->GRID_X+moving_frame->move);i++)
		{
			fprintf(fp, "%10.8f ", moving_frame_get(moving_frame, i, j));
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
}

void moving_frame_print_charge(moving_frame* moving_frame, char* pattern, int counter)
{
	char filename[30];
	sprintf(filename, pattern, counter);
	FILE *fp;
	fp=fopen(filename, "w");
	int i, j;
	for(j=0;j<moving_frame->GRID_Y;j++)
	{
		for(i=0;i<moving_frame->GRID_X;i++)
		{
			fprintf(fp, "%10.8f ", moving_frame->elements[i*(moving_frame->GRID_Y)+j]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
}

void moving_frame_set_charge(moving_frame* moving_frame, int lab_x, int lab_y, double value)
{
	int x_position=lab_x-(moving_frame->move);
	int y_position=lab_y;
	moving_frame->elements[x_position*(moving_frame->GRID_Y)+y_position]=value;
}

double moving_frame_get_charge(moving_frame* moving_frame, int lab_x, int lab_y)
{
	int x_position=lab_x-moving_frame->move;
	int y_position=lab_y;		
	return moving_frame->elements[x_position*(moving_frame->GRID_Y)+y_position];
}

void moving_frame_get_guard_cells(moving_frame* moving_frame)
{
	int j;
	for(j=0;j<moving_frame->GRID_Y;j++)
	{
		moving_frame->guard_cells[2*moving_frame->GRID_Y+j]=moving_frame->elements[(moving_frame->move % moving_frame->GRID_X)*moving_frame->GRID_Y+j];
		moving_frame->guard_cells[3*moving_frame->GRID_Y+j]=moving_frame->elements[((moving_frame->GRID_X-1+moving_frame->move) % moving_frame->GRID_X)*moving_frame->GRID_Y+j];
	}
}

void moving_frame_set_guard_cells(moving_frame* moving_frame)
{
	int j;
	for(j=0;j<moving_frame->GRID_Y;j++)
	{
		moving_frame->guard_cells[j]=moving_frame->guard_cells[2*moving_frame->GRID_Y+j];
		moving_frame->guard_cells[moving_frame->GRID_Y+j]=moving_frame->guard_cells[3*moving_frame->GRID_Y+j];
	}
}
