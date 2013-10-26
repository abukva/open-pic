export CC=mpicc
all:

	$(CC) -c diagnostic.c field.c frame.c init.c main.c moving_frame.c pusher.c simulation.c solver.c
	$(CC) diagnostic.o field.o frame.o init.o main.o moving_frame.o pusher.o simulation.o solver.o -o simulacija

