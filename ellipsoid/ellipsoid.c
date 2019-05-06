#include <string.h>
#include "ellipsoid.h"

/* some parameters */
static const double EPSILON = 1e-9;
static const double ALPHA   = 0.25;
static const double BETA    = 0.5;

/**
 * returns p^T Q p 
 * flop count = 12adds + 12mults
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
 * returns \sum_i (p_i^T Q p_i - 1)^2
 * flop count = n*(3 adds + 1mult  + 2*C(_quadratic_form))
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
 * flop count PER iter_step = 
 * 				n*( C(_quadratic_form)+1add + 9mults + 9*(2mults+1adds) ) + 
 *        9mults + 
 *        9*(1add + 1mult) + 9*(1mult+1add) + 
 *					iters_bt*(C(_cost)*2 + 1add + 2mult + 1mult + 9*(1add + 1mult)) +
 *        9*(1add + 3 mults)
 *	where iters_bt = number of back tracking iterations
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

/*
 * convenience function, wraps fit_ellipsoid above, given the mils (lengths) along
 * each of the DIRECTIONS defined above
 *
 */
void fit_ellipsoid_mils(const double *mils, double Q[3][3])
{
	/* construct the points */
	double p[NUM_DIRECTIONS][3];
	for (int i=0; i<NUM_DIRECTIONS; i++) {
		for (int j=0; j<3; j++) {
			p[i][j] = mils[i] * DIRECTIONS_NORMALIZED[i][j];
		}
	}

	/* call the main routine */
	fit_ellipsoid(p, NUM_DIRECTIONS, Q);
}