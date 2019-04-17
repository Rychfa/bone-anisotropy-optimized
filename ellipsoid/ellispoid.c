#include "tsc_x86.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define NUM_RUNS 1
// #define CALIBRATE
#define CYCLES_REQUIRED 1e8

/* some parameters */
static const double EPSILON = 1e-9;
static const double ALPHA   = 0.25;
static const double BETA    = 0.5;


/**
 *
 */ 
static double _quadratic_form(const double p[3], const double Q[3][3])
{
	double Qp[3] = {0};

	for (int i=0; i<3; i++){
		for (int j=0; j<3; j++) {
			Qp[i] += Q[i][j] * p[j];
		}
	}

	double ret = 0;
	for (int i=0; i<3; i++) {
		ret += p[i] * Qp[i];
	}

	return ret;
}

/**
 *
 */
static double _cost(const double (*p)[3], int n, const double Q[3][3])
{
	double ret = 0;

	for (int i=0; i<n; i++){
		ret += (_quadratic_form(p[i], Q) - 1)*(_quadratic_form(p[i], Q) - 1);
	}

	return ret;
}

/**
 * Gradient descent with backtracking line search [Boyd Chapter 9] 
 * applied to the ellisoid fitting problem:
 * Q = argmin \sum_i (p_i^T Q p_i - 1)^2
 */
void fit_ellipsoid(const double (*p)[3], int n, double Q[3][3])
{
	double Qk[3][3] = {{1,0,0}, {0,1,0}, {0,0,1}}; 

	while(1) {
		/* determine gradient */

		double grad[3][3] = {{0}};

		for (int i=0; i<n; i++) {
			double residual_i = _quadratic_form(p[i], Qk) - 1;
			double jacobian_i[3][3] = {{0}};
		
			/* jacobian_i =  outer product */
			for (int j=0; j<3; j++) {
				for (int k=0; k<3; k++) {
					jacobian_i[j][k] = p[i][j]*p[i][k];
				}
			}

			/* grad += 2*residual_i*jacobian_i */
			for (int j=0; j<3; j++) {
				for (int k=0; k<3; k++) {
					grad[j][k] += 2*residual_i*jacobian_i[j][k];
				}
			}
		}

		/* take step direction to be negative gradient */
		
		double step[3][3];
		for (int i=0; i<3; i++) {
			for (int j=0; j<3; j++) {
				step[i][j] = -grad[i][j];
			}
		}

		/* backtracking line search */
		
		double t = 1;
		double Qk_plus_tstep[3][3] = {{0}};
		for (int i=0; i<3; i++) {
			for (int j=0; j<3; j++) {
				Qk_plus_tstep[i][j] = Qk[i][j] + t*step[i][j];
			}
		}
		double trace_gradstep = 0;
		for (int i=0; i<3; i++){
			for (int j=0; j<3; j++) {
				trace_gradstep += grad[i][j]*step[j][i];
			}
		}
		while (_cost(p, n, Qk_plus_tstep) > _cost(p, n, Qk) + ALPHA*t*trace_gradstep) {
			t = t*BETA;

			for (int i=0; i<3; i++) {
				for (int j=0; j<3; j++) {
					Qk_plus_tstep[i][j] = Qk[i][j] + t*step[i][j];
				}
			}
		}

		/* stopping critera check */

		double normFro_tstep = 0;
		for (int i=0; i<3; i++) {
			for (int j=0; j<3; j++) {
				normFro_tstep += (t*step[i][j])*(t*step[i][j]);
			}
		}

		if (normFro_tstep < EPSILON) break;

		memcpy(Qk, Qk_plus_tstep, sizeof(Qk));
	} /* main while loop */

	memcpy(Q, Qk, sizeof(Qk));
}

#if 0
void fill_vector(double * x, int n) {
    for(int i=0; i < n; i++) {
        x[i] = (double) rand();
    }
}

double rdtsc(double *A, double *B, double *C, long n, 
	void (*func)(double *, double *, double *, long)) 
{
    int i, num_runs;
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
    while(num_runs < (1 << 14)) {
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
#endif

int main(int argc, char **argv)
{
	/* initialize data */
	double Ipoints[][3] = 
		{{-0.28687072, -0.5844695 , -0.29465669},
        { 0.29664863,  0.62228291,  0.36936989},
        { 0.01220042,  0.04451405,  0.28935352},
        { 0.03887869, -0.06267655, -0.02044771},
        { 0.16627267,  0.62938256, -0.08679698},
        { 0.2123087 ,  0.55534194, -0.09466649},
        { 0.26719668,  0.47538131,  0.58769478},
        { 0.28394338,  0.54458556,  0.51808755},
        {-0.01285678,  0.39130018, -0.55237714},
        { 0.08337157,  0.20009325, -0.22645592}};

  /* test correctness. Q should be according to CVXPY
  	[[239.24816358 -81.82197492 -47.58347817]
 		[-81.82197492  29.89943344  16.36709216]
 		[-47.58347817  16.36709216  10.84916192]]
 	*/
  double Q[3][3];
  int n = sizeof(Ipoints)/sizeof(double)/3;
  printf("%d\n", n);
  fit_ellipsoid(Ipoints, n, Q);
  for (int i=0; i<3; i++) {
  	for (int j=0; j<3; j++) {
  		printf("%f, ", Q[i][j]);
  	}
  	printf("\n");
  }

	/* perform timings */
	// double cycles = rdtsc(A, B, C, n, mm1);
	// printf("mm1: %g, %g\n", cycles, (2*n*n*n)/cycles);

	// cycles = rdtsc(A, B, C, n, mm2);
	// printf("mm2: %g, %g\n", cycles, (2*n*n*n)/cycles);
}