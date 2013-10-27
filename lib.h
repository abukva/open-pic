typedef struct particle_s {
	double x,y,px,py,pz;
	double q, q_weight;
} particle; /*strucutre of a particle */

typedef struct field_s {
	double ex,ey,ez,bx,by,bz; /* components of a electric and magnetic field */
} field;

typedef struct moving_frame_s {
	double *elements;
	double *guard_cells;
	int GRID_X;
	int GRID_Y;
	int move;
} moving_frame;

typedef struct grid_s {
	moving_frame *grid_ex, *grid_ey, *grid_ez, *grid_bx, *grid_by, *grid_bz, *grid_jx, *grid_jy, *grid_jz, *charge;
} grid; /* structure of fields at grid */

typedef struct split_s {
	double x0[3];
	double x1[3];
	double y0[3];
	double y1[3];
	int i_cell[3];
	int j_cell[3];
} split;

particle *init_particles(int, int, int, int, double, double, int*, int, int, int, double); /* initializations of particles positions at beggining */

field *particle_field(double, double, grid *, double, double, int, int); /* return field at particle position */

grid *init_grid(int , int, double, double, double, double, double, double, double, double); /* initialization of grid */

int number_of_itter(double,double); /*return number of iterations */

void push_one_set(particle *, double, int*, grid *, double, double, int, int); /*push all particles once and calculate the currents*/

void simulation(particle *, double , double , int*, grid * ,int, double, double, int, int, int, int, int, int, int, int, int, int, int, int, double, int, int, int, int); /* push all particles number of itter */

void solve_fields(grid *, double, double, double, int, int); /* solve fields */

double initial_field_value(double , double, double, double, double, double, double, double); /* function for the initial value of the field */

void current_deposition (grid *, double, double, double, double, double, double, double, double, int, int, double, double); /* function for depositin currents */

void boundary_conditions_particles(particle *particle, int, int, int*, int*, double, double, grid*, double, double, double); /* boundary conditions for particles */

void get_weight(particle*, double, double, int, int*, int*, double);

void move_frame(particle *, grid *, int*, int, int, int, int, double, double, double);

void field_value_print(grid*, int, int, int, int, int, int, int, int, int, int);

void deposit_charge(particle*, grid*, int, int, int, double, double);

void get_splits(split*, int*, double, double, double, double, double, double);

void get_coordinates(double, double, double, double, int*, int*, int*, int*, double, double);

int gcd(int, int);

moving_frame* moving_frame_create(int, int);

double moving_frame_get(moving_frame*, int, int);

void moving_frame_set(moving_frame*, int, int, double);

int moving_frame_get_move(moving_frame*);

void moving_frame_move(moving_frame*);

void moving_frame_zero(moving_frame*);

void moving_frame_print(moving_frame*, char*, int);

void moving_frame_print_charge(moving_frame*, char*, int);

void moving_frame_set_charge(moving_frame*, int, int, double);

double moving_frame_get_charge(moving_frame*, int, int);

void moving_frame_get_guard_cells(moving_frame*);

void moving_frame_set_guard_cells(moving_frame*);

void particle_value_print(particle*, int, int, int);

double abs_double(double);
