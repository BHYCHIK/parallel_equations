/*
 ============================================================================
 Name        : rotation.c
 Author      : Ivan Remen
 Version     :
 Copyright   : Your copyright notice
 Description : Hello MPI World in C 
 ============================================================================
 */
#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include <math.h>

#define N 12

static void print_equation(int matr[N][N + 1]) {
	int i = 0, j = 0;
	for (i = 0; i < N; ++i) {
		for (j = 0; j < N + 1; ++j) {
			fprintf(stdout, "%d ", matr[i][j]);
		}
		fprintf(stdout, "\n");
	}
}

static void fill_equation(int matr[N][N + 1]) {
	int i = 0, j = 0;
	for (i = 0; i < N; ++i) {
		for (j = 0; j < N + 1; ++j) {
			fscanf(stdin, "%d", &matr[i][j]);
		}
	}
}

static void operate_equations(int *first_equation, int matr[][N+1], int eq_num, int n) {
	int i, j;
	for (i = 0; i < eq_num; ++i) {
		//if (n == 0) printf("%d)", n);
		for (j = 0; j < N + 1; ++j) {
			//if (n == 0) printf("%dx%d=%d ", i,j, matr[i][j]);
			++matr[i][j];
		}
		//if (n == 0) printf("\n");
	}
}

static void solve_recursively_server(int proc_id, int proc_num, int matr[][N+1], int iter) {
	//MPI_Barrier( MPI_COMM_WORLD);
	if (iter == N + 1) return;
	MPI_Bcast(matr, N + 1, MPI_INT, 0, MPI_COMM_WORLD);
	int *first_equation = ((int *)matr);

	int equations_per_node = ceil((N - iter) / ( (proc_num - 0) * 1.0));

	int i = 0;
	int equation_sended = 0;
	for (i = 1; i < proc_num; ++i) {
		int equations_left = N - iter - equation_sended;
		int to_send = ((equations_left <  equations_per_node) ?
				equations_left : equations_per_node);
		if (to_send == 0) break;
		MPI_Send(first_equation + (equation_sended + 1) * (N + 1),
				(N + 1) * to_send,
				MPI_INT, i, 0, MPI_COMM_WORLD);
		//printf("%d) %d Iteration %d\n", i, to_send, iter);
		equation_sended += to_send;
	}
	int equations_left = (N - iter) - equation_sended;
	//printf("LEFT: %d\n\n", equations_left);

	operate_equations(first_equation, first_equation + (N - iter - equations_left + 1) * (N + 1),
			equations_left, proc_id);

	int equation_received = 0;
	for (i = 1; i < proc_num; ++i) {
		int equations_left = N - iter - equation_received;
		int to_recv = ((equations_left <  equations_per_node) ?
				equations_left : equations_per_node);
		if (to_recv == 0) break;

		MPI_Status status;
		MPI_Recv(first_equation + (equation_received + 1) * (N + 1),
				(N + 1) * to_recv,
				MPI_INT, i, 0, MPI_COMM_WORLD, &status);
		//printf("%d) %d Iteration %d\n", i, to_send, iter);
		equation_received += to_recv;
	}

	solve_recursively_server(proc_id, proc_num, first_equation + N + 1, iter + 1);
}

static void solve_recursively_client(int proc_id, int proc_num, int matr[][N+1], int iter) {
	//MPI_Barrier( MPI_COMM_WORLD);
	if (iter == N + 1) return;

	static int first_equation[N + 1];
	MPI_Bcast(first_equation, N + 1, MPI_INT, 0, MPI_COMM_WORLD);

	int equations_per_node = ceil((N - iter) / ( (proc_num - 0) * 1.0));
	int total_equations = N - iter;
	int start_of_equations = (proc_id - 1) * equations_per_node;
	if (total_equations < start_of_equations) {
		solve_recursively_client(proc_id, proc_num, matr, iter + 1);
		return;
	}

	int need_to_solve_equations =
			(total_equations < start_of_equations + equations_per_node) ?
					total_equations - start_of_equations : equations_per_node;

	if (need_to_solve_equations == 0) {
		solve_recursively_client(proc_id, proc_num, matr, iter + 1);
		return;
	}

	//printf("Iteration %d. I am %d. Solving from %d. Need to solve %d\n", iter, proc_id, start_of_equations,
			//need_to_solve_equations);

	MPI_Status status;
	MPI_Recv(matr, need_to_solve_equations * (N+1), MPI_INT, 0, 0, MPI_COMM_WORLD, &status);

	operate_equations(first_equation, matr, need_to_solve_equations, proc_id);

	MPI_Send(matr, (N + 1) * need_to_solve_equations, MPI_INT, 0, 0, MPI_COMM_WORLD);

	solve_recursively_client(proc_id, proc_num, matr, iter + 1);
}

static void solve_recursively(int proc_id, int proc_num, int matr[N][N+1], int iter) {
	if (proc_id == 0) solve_recursively_server(proc_id, proc_num, matr, iter);
	else solve_recursively_client(proc_id, proc_num, matr, iter);
}

static void solve_equation(int proc_id, int proc_num) {
	static int matr[N][N+1];
	if (proc_id == 0) {
		fill_equation(matr);
		print_equation(matr);
	}

	solve_recursively(proc_id, proc_num, matr, 1);
	if (proc_id == 0) {
		printf("\n\n\n\n\n\n\n");
		print_equation(matr);
	}
}

int main(int argc, char* argv[]){
	int  proc_id; /* rank of process */
	int  proc_num;       /* number of processes */
	
	/* start up MPI */
	
	MPI_Init(&argc, &argv);
	
	/* find out process rank */
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
	
	/* find out number of processes */
	MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
	
	solve_equation(proc_id, proc_num);
	
	/* shut down MPI */
	MPI_Finalize(); 
	
	
	return 0;
}
