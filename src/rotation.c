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
#include <assert.h>

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

static void transform_matrix2(int proc_id, double matr[][N], int columns_per_proc) {
	int i, j, k;
	for (i = 0; i < N - 1; ++i) {

		int worker_with_i = i / columns_per_proc;
		static double column_i[N];

		if (worker_with_i == proc_id) {
			//printf("YYY%dYYY\n",proc_id);
			for (j = 0; j < N; ++j) {
				column_i[j] = matr[i % columns_per_proc][j];
				//printf("XXX%lfXXX\n", matr[i % columns_per_proc][j]);
			}
		}
		//printf("worker_with_id %d\n", worker_with_i);
		MPI_Bcast(column_i, N, MPI_DOUBLE, worker_with_i, MPI_COMM_WORLD);
		//printf("i=%d proc=%d ", i, proc_id);
		for (j = 0; j < N; ++j) {
			//printf("%lf ", column_i[j]);
		}
		//printf("\n");

		for (j = i + 1; j < N; ++j) {

			double c = 0.0, s = 0.0;

			double buf[2];

			//printf("CALCING: %lf %lf %lf\n", column_i[i], column_i[i], column_i[j]);
			c = column_i[i] / sqrt(sqr(column_i[i]) + sqr(column_i[j]));
			s = column_i[j] / sqrt(sqr(column_i[i]) + sqr(column_i[j]));

			static double first_equation_copy[N + 1];
			static double second_equation_copy[N + 1];

			//printf("BEFORE2: %lf %lf\n", column_i[j], matr[i][j]);
			//printf("BEFORE1: %lf %lf\n", column_i[i], matr[i][i]);


			for (k = 0; k < columns_per_proc; ++k) {
				first_equation_copy[k] = matr[k][i];
				second_equation_copy[k] = matr[k][j];

				matr[k][i] = first_equation_copy[k] * c + second_equation_copy[k] * s;
				matr[k][j] = second_equation_copy[k] * c - first_equation_copy[k] * s;
			}

			double tmp_column_i_i = column_i[i];
			double tmp_column_i_j = column_i[j];
			printf("c=%lf s=%lf\n", c, s);
			column_i[i] = tmp_column_i_i * c + tmp_column_i_j * s;
			column_i[j] = tmp_column_i_j * c - tmp_column_i_i * s;

			if (column_i[i] != matr[i][i]) {
				//printf("FUCKOF1: %d %d %lf %lf\n", i, j, column_i[i], matr[i][i]);
			}
			if (column_i[j] != matr[i][j]) {
				//printf("FUCKOF2: %d %d %lf %lf\n", i, j, column_i[j], matr[i][j]);
			}
			//assert(column_i[i] == matr[i][i]);
			//assert(column_i[j] == matr[i][j]);
		}
	}

}

static void solve(int proc_id, int proc_num, double matr[N + 1][N]) {
	int i, j;
	int columns_per_proc = ceil((N + 1.0) / proc_num);
	printf("columns_per_proc=%d\n", columns_per_proc);

	void *work_buffer = malloc(columns_per_proc * N * sizeof(double));

	MPI_Scatter(matr, columns_per_proc * N, MPI_DOUBLE, work_buffer,
			columns_per_proc * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	transform_matrix2(proc_id, work_buffer, columns_per_proc);

	MPI_Gather(work_buffer, columns_per_proc * N, MPI_DOUBLE, matr, columns_per_proc * N,
			MPI_DOUBLE, 0, MPI_COMM_WORLD);
	free(work_buffer);

	if (proc_id == 0) {

		for(i = 0; i < N; ++i) {
			for (j = 0; j < N + 1; ++j) {
				printf("%lf ", matr[j][i]);
			}
			printf("\n");
		}

	}

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

	if (proc_id == 0) printf("Used %lf seconds\n", stop - start);

	/* shut down MPI */
	MPI_Finalize(); 

	return 0;
}
