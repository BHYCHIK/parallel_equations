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
#include <stdlib.h>

#define N 3

static void fill_equation(double matr[N + 1][N]) {
	int i = 0, j = 0;
	for (i = 0; i < N; ++i) {
		for (j = 0; j < N + 1; ++j) {
			fscanf(stdin, "%lf", &matr[j][i]);
		}
	}
}

static double sqr(double a) {
	return a * a;
}

static void transform_matrix(int proc_id, double matr[][N], int columns_per_proc) {
	int i, j, k;
	for (i = 0; i < N - 1; ++i) {
		for (j = i + 1; j < N; ++j) {

			double c = 0.0, s = 0.0;

			int worker_with_i = i / columns_per_proc;
			double buf[2];

			buf[0] = matr[i % columns_per_proc][i] / sqrt(sqr(matr[i % columns_per_proc][i]) + matr[i % columns_per_proc][j]);
			buf[1]= matr[i % columns_per_proc][j] / sqrt(sqr(matr[i % columns_per_proc][i]) + matr[i % columns_per_proc][j]);

			MPI_Bcast(&buf, 2, MPI_DOUBLE, worker_with_i, MPI_COMM_WORLD);

			c = buf[0];
			s = buf[1];

			static double first_equation_copy[N + 1];
			static double second_equation_copy[N + 1];

			for (k = 0; k < columns_per_proc; ++k) {
				first_equation_copy[k] = matr[k][i];
				second_equation_copy[k] = matr[k][j];

				matr[k][i] = first_equation_copy[k] * c + second_equation_copy[k] * s;
				matr[k][j] = second_equation_copy[k] * c - first_equation_copy[k] * s;
			}
		}
	}
}

static void solve(int proc_id, int proc_num, double matr[N + 1][N]) {
	int i, j;
	int columns_per_proc = ceil((N + 1.0) / proc_num);

	void *work_buffer = malloc(columns_per_proc * N * sizeof(double));

	MPI_Scatter(matr, columns_per_proc * N, MPI_DOUBLE, work_buffer,
			columns_per_proc * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	transform_matrix(proc_id, work_buffer, columns_per_proc);

	MPI_Gather(work_buffer, columns_per_proc * N, MPI_DOUBLE, matr, columns_per_proc * N,
			MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if (proc_id == 0) {

		for(i = 0; i < N; ++i) {
			for (j = 0; j < N + 1; ++j) {
				printf("%lf ", matr[j][i]);
			}
			printf("\n");
		}

	}

	free(work_buffer);
}

static void solve_equation(int proc_id, int proc_num) {
	static double matr[N + 1][N];
	if (proc_id == 0) {
		fill_equation(matr);
	}

	solve(proc_id, proc_num, matr);
}

int main(int argc, char* argv[]){
	int  proc_id; /* rank of process */
	int  proc_num;       /* number of processes */

	/* start up MPI */

	MPI_Init(&argc, &argv);

	double start = MPI_Wtime();


	/* find out process rank */
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);

	/* find out number of processes */
	MPI_Comm_size(MPI_COMM_WORLD, &proc_num);

	solve_equation(proc_id, proc_num);

	double stop = MPI_Wtime();

	if (proc_id == 0) printf("Used %lf seconds", stop - start);

	/* shut down MPI */
	MPI_Finalize(); 

	return 0;
}
