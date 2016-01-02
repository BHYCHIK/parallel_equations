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

#define N 2
#define PROC_NUM 2

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

static void solve_recursively_server(int proc_id, int proc_num, int matr[N][N+1], int iter) {
	if (iter == N + 1) return;
	MPI_Bcast(matr, N + 1, MPI_INT, 0, MPI_COMM_WORLD);
	int *first_equation = ((int *)matr);

	solve_recursively_server(proc_id, proc_num, first_equation + N + 1, iter + 1);
}

static void solve_recursively_client(int proc_id, int proc_num, int matr[N][N+1], int iter) {
	if (iter == N + 1) return;

	static int first_equation[N + 1];
	MPI_Bcast(first_equation, N + 1, MPI_INT, 0, MPI_COMM_WORLD);
	int i = 0;
	for (i = 0; i < N + 1; ++i)
		printf("%d ", first_equation[i]);
	printf("\n");
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
