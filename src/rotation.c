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

#define N 3

static void print_equation(double matr[N][N + 1]) {
	int i = 0, j = 0;
	for (i = 0; i < N; ++i) {
		for (j = 0; j < N + 1; ++j) {
			fprintf(stdout, "%lf ", matr[i][j]);
		}
		fprintf(stdout, "\n");
	}
}

static void fill_equation(double matr[N][N + 1]) {
	int i = 0, j = 0;
	for (i = 0; i < N; ++i) {
		for (j = 0; j < N + 1; ++j) {
			fscanf(stdin, "%lf", &matr[i][j]);
		}
	}
}

static double sqr(double a) {
	return a * a;
}

static find_solution(double matr[][N+1], double x_solution[N + 1]) {
	int i, j;
	for (i = N - 1; i >= 0; --i) {
		double tmp = matr[i][N];
		for (j = N - 1; j > i; --j) {
			tmp -= x_solution[j] * matr[i][j];
		}
		x_solution[i] = tmp / matr[i][i];
	}
}

static void operate(double *first_equation, double *second_equation, double c, double s, int columns_per_proc) {
	/*
	int i;
	for (i = 0; i < columns_per_proc; ++i)
		second_equation[i]++;
	*/


	static double first_equation_copy[N + 1];
	static double second_equation_copy[N + 1];

	int i = 0;
	for (i = 0; i < columns_per_proc; ++i) {
		first_equation_copy[i]= first_equation[i];
		second_equation_copy[i]= second_equation[i];

		first_equation[i] = first_equation_copy[i] * c + second_equation_copy[i] * s;
		second_equation[i] = second_equation_copy[i] * c - first_equation_copy[i] * s;
	}

}

static void solve(int proc_id, int proc_num, double matr[N][N+1]) {
	int i, j;
	int columns_per_proc = ceil((N + 1.0) / proc_num);
	//printf("columns_per_proc: %d\n", columns_per_proc);
	for (i = 0; i < N - 1; ++i) {
		for (j = i + 1; j < N; ++j) {

			double c = 0;
			double s = 0;
			{
				double buf[2];
				buf[0] = matr[i][i] / sqrt(sqr(matr[i][i]) + matr[j][i]);
				buf[1] = matr[j][i] / sqrt(sqr(matr[i][i]) + matr[j][i]);

				MPI_Bcast(&buf, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				c = buf[0];
				s = buf[1];

				//if (proc_id == 0) printf("%d %d %lf %lf\n", i, j, c, s);
			}

			double *mem = NULL;

			static double first_equation[N + 1];
			mem = ((double *)matr) + i * (N + 1);
			MPI_Scatter(mem, columns_per_proc, MPI_DOUBLE, first_equation, columns_per_proc,
					MPI_DOUBLE, 0, MPI_COMM_WORLD);

			static double second_equation[N + 1];
			mem = ((double *)matr) + j * (N + 1);
			MPI_Scatter(mem, columns_per_proc, MPI_DOUBLE, second_equation, columns_per_proc,
						MPI_DOUBLE, 0, MPI_COMM_WORLD);

			operate(first_equation, second_equation, c, s, columns_per_proc);

			MPI_Gather(second_equation, columns_per_proc, MPI_DOUBLE, mem,
					columns_per_proc, MPI_DOUBLE, 0, MPI_COMM_WORLD);

			mem = ((double *)matr) + i * (N + 1);
			MPI_Gather(first_equation, columns_per_proc, MPI_DOUBLE, mem,
								columns_per_proc, MPI_DOUBLE, 0, MPI_COMM_WORLD);

			/*if (proc_id == 0) {
				printf("\n\n\n\n\n\n\n");
				print_equation(matr);
			}*/
		}
	}

	if (proc_id != 0) return;
	static double x_solution[N + 1];
	find_solution(matr, x_solution);
	printf("\n\n\n\n\n\n\n");
	print_equation(matr);
	printf("\n\n\n\n\n\n\n");
	for(i = 0; i < N; ++i) {
		printf("x_solution[%i] = %lf\n", i, x_solution[i]);
	}
}

static void solve_equation(int proc_id, int proc_num) {
	static double matr[N][N+1];
	if (proc_id == 0) {
		fill_equation(matr);
		print_equation(matr);
	}

	solve(proc_id, proc_num, matr);
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
