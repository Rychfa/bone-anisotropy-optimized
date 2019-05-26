#include "tsc_x86.h"
#include <stdlib.h>
#include <stdio.h>

#define NUM_RUNS 5
#define CALIBRATE
#define CYCLES_REQUIRED 1e10

void mm1(double *A, double *B, double *C, int n)
{
	int i,j,k;

	double (*a)[n][n] = (double (*)[n][n]) A;
	double (*b)[n][n] = (double (*)[n][n]) B;
	double (*c)[n][n] = (double (*)[n][n]) C;

	for (i=0; i<n; i++) {
		for (j=0; j<n; j++) {
			for (k=0; k<n; k++) {
				c[i][j][k] = a[i][j][k]*b[i][j][k];
			}
		}
	}
}

/* 
 */
void mm2(double *A, double *B, double *C, int n)
{
	int i,j,k;

	double (*a)[n][n] = (double (*)[n][n]) A;
	double (*b)[n][n] = (double (*)[n][n]) B;
	double (*c)[n][n] = (double (*)[n][n]) C;

	for (i=0; i<n; i++) {
		for (j=0; j<n; j++) {
			for (k=0; k<n; k+=4) {
				double a0 = a[i][j][k],
						   a1 = a[i][j][k+1],
						   a2 = a[i][j][k+2],
						   a3 = a[i][j][k+3];
				double b0 = b[i][j][k],
							 b1 = b[i][j][k+1],
							 b2 = b[i][j][k+2],
							 b3 = b[i][j][k+3];

				double c0 = a0*b0,
							 c1 = a1*b1,
							 c2 = a2*b2,
							 c3 = a3*b3;

				c[i][j][k]   = c0;
				c[i][j][k+1] = c1;
				c[i][j][k+2] = c2;
				c[i][j][k+3] = c3;
			}
		}
	}
}


void fill_vector(double * x, int n) {
    for(int i=0; i < n; i++) {
        x[i] = 1; //(double) rand() / RAND_MAX;
    }
}

double rdtsc(double *A, double *B, double *C, int n, 
	void (*func)(double *, double *, double *, int)) 
{
    myInt64 i, num_runs;
    myInt64 cycles;
    myInt64 start;
    num_runs = NUM_RUNS;

    /* 
     * The CPUID instruction serializes the pipeline.
     * Using it, we can create execution barriers around the code we want to time.
     * The calibrate section is used to make the computation large enough so as to 
     * avoid measurements bias due to the timing overhead.
     */
#ifdef CALIBRATE
    while(1) {
        start = start_tsc();
        for (i = 0; i < num_runs; ++i) {
        	func(A, B, C, n);
        }
        cycles = stop_tsc(start);

        if(cycles >= CYCLES_REQUIRED) break;

        num_runs *= 2;
    }
#endif

    start = start_tsc();
    for (i = 0; i < num_runs; ++i) {
        func(A, B, C, n);
    }

    cycles = stop_tsc(start)/num_runs;
    return (double) cycles;
}

int main(int argc, char **argv)
{
	int n = 8;
	if (argc > 1) {
		n = atoi(argv[1]);
	}

	/* initialize data */
	double *A = malloc(sizeof(double)*n*n*n);
	double *B = malloc(sizeof(double)*n*n*n);
	double *C = malloc(sizeof(double)*n*n*n);

	fill_vector(A, n*n*n);
	fill_vector(B, n*n*n);
	fill_vector(C, n*n*n);

	/* perform timings */
	double cycles = rdtsc(A, B, C, n, mm1);
	printf("%f, %f, %f\n", A[0], B[0], C[0]);
	printf("mm1: %g, %g\n", cycles, (n*n*n)/cycles);

	cycles = rdtsc(A, B, C, n, mm2);
	printf("%f, %f, %f\n", A[0], B[0], C[0]);
	printf("mm2: %g, %g\n", cycles, (n*n*n)/cycles);

	free(A);
	free(B);
	free(C);
}