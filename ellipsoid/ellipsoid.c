#include <string.h>
#include <math.h>
#include "ellipsoid.h"
#ifdef DEBUG
	#include <stdio.h>
#endif

/* some parameters */
static const double EPSILON = 1e-3;
static const double ALPHA   = 0.25;
static const double BETA    = 0.5;

/*****************************************************************************
 *  Private data
 ****************************************************************************/
#ifdef DEBUG
	static unsigned long ellipsoid_flop_count  = 0;
	static unsigned long ellipsoid_outer_iters = 0;
#endif

/*****************************************************************************
 *  Debug tools
 ****************************************************************************/

#ifdef DEBUG
	void fit_ellipsoid_debug_init(void)
	{
		ellipsoid_flop_count  = 0;
		ellipsoid_outer_iters = 0;
	}
	long fit_ellipsoid_debug_deinit(void)
	{
		printf("[ellipsoid] flop count = %ld\n", ellipsoid_flop_count);
		printf("[ellipsoid] iters      = %ld\n", ellipsoid_outer_iters);
		return ellipsoid_flop_count;
	}
#endif

/*****************************************************************************
 *  Initial implementation
 ****************************************************************************/

/**
 * returns p^T Q p 
 * flop count = 8adds + 12mults
 *					  = 20
 */ 
static double _quadratic_form(const double p[3], const double (*Q)[3][3])
{
	double Qp[3] = {0};

	for (int i=0; i<3; i++) {
		for (int j=0; j<3; j++) {
			Qp[i] += (*Q)[i][j] * p[j];
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
 * flop count = n*(C(_quadratic_Form)+1sub + 1add+1mult)
 *            = n*(20+3) = 23n
 */
static double _cost(const double (*p)[3], int n, const double (*Q)[3][3])
{
	double ret = 0;

	for (int i=0; i<n; i++){
		double residual_i = _quadratic_form(p[i], Q) - 1;
		ret += residual_i*residual_i;
	}
	return ret;
}

/*
 * ||A||_F
 * flop count = 9*(1mult + 1add) + 1sqrt
 */
static double _normFro(const double A[3][3])
{
	double ret = 0;
	for (int i=0; i<3; i++) {
		for (int j=0; j<3; j++) {
			ret += (A[i][j])*(A[i][j]);
		}
	}
	return sqrt(ret);
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
void fit_ellipsoid(const double (*p)[3], int n, double (*Q)[3][3])
{
	(*Q)[0][0] = 1; (*Q)[0][1] = 0; (*Q)[0][2] = 0;
	(*Q)[1][0] = 0; (*Q)[1][1] = 1; (*Q)[1][2] = 0;
	(*Q)[2][0] = 0; (*Q)[2][1] = 0; (*Q)[2][2] = 1;

	while(1) {
		/* determine gradient */

		double grad[3][3] = {{0}};

		for (int i=0; i<n; i++) {
			double residual_i = _quadratic_form(p[i], Q) - 1;
			double jacobian_i[3][3] = {{0}};
		
			/* jacobian_i =  outer product */
			for (int j=0; j<3; j++) {
				for (int k=0; k<3; k++) {
					jacobian_i[j][k] = p[i][j]*p[i][k];
				} 
			} 

			/* grad += 2*residual_i*jacobian_i */
			// coeff = 2*residual_i
			for (int j=0; j<3; j++) {
				for (int k=0; k<3; k++) {
					grad[j][k] += 2*residual_i*jacobian_i[j][k];
					// grad[j][k] += coeff*jacobian_i[j][k];
				}
			}
		}

#ifdef DEBUG
		ellipsoid_flop_count += n*(20+1 + 9+1+9*2 +2); 
		ellipsoid_flop_count += 9*(1+1) + 1; 
		ellipsoid_outer_iters++;
#endif

		if (_normFro(grad) < EPSILON) break;

		/* take step direction to be negative gradient */
		
		double step[3][3];
		for (int i=0; i<3; i++) {
			for (int j=0; j<3; j++) {
				step[i][j] = -grad[i][j];
			}
		}
#ifdef DEBUG
		// ellipsoid_flop_count += 9;
#endif

		/* backtracking line search */
		
		double t = 100;
		double Qk_plus_tstep[3][3] = {{0}};
		for (int i=0; i<3; i++) {
			for (int j=0; j<3; j++) {
				Qk_plus_tstep[i][j] = (*Q)[i][j] + t*step[i][j];
			}
		}
#ifdef DEBUG
		ellipsoid_flop_count += 9*(1+1);
#endif

		double trace_gradstep = 0;
		for (int i=0; i<3; i++) {
			for (int j=0; j<3; j++) {
				trace_gradstep += grad[i][j]*step[j][i];
			}
		}
#ifdef DEBUG
		// ellipsoid_flop_count += 9*(1+1); 
		ellipsoid_flop_count += 23*n;// + 1+2; 
#endif
		while (_cost(p, n, (double (*)[3][3]) Qk_plus_tstep) > 
																_cost(p, n, Q) + ALPHA*t*trace_gradstep) {
			t = t*BETA;
			
			for (int i=0; i<3; i++) {
				for (int j=0; j<3; j++) {
					Qk_plus_tstep[i][j] = (*Q)[i][j] + t*step[i][j];
				}
			}
#ifdef DEBUG
			/* not counting the second cost, it is loop invariant code */
			ellipsoid_flop_count += 1 + 9*(1+1) + 23*n + 1+2;
#endif
		} /* back tracking while loop */
		memcpy(Q, Qk_plus_tstep, sizeof(double)*9);
	} /* main while loop */
}

/*
 * convenience function, wraps fit_ellipsoid above, given the mils (lengths) along
 * each of the DIRECTIONS defined above
 * flop count = NUM_DIRS*3mults + C(fit_ellipsoid)
 */
void fit_ellipsoid_mils(const double *mils, double (*Q)[3][3])
{
	/* construct the points */
	double p[NUM_DIRECTIONS][3];
	for (int i=0; i<NUM_DIRECTIONS; i++) {
		for (int j=0; j<3; j++) {
			p[i][j] = mils[i] * DIRECTIONS_NORMALIZED[i][j];
		}
	}
#ifdef DEBUG
		ellipsoid_flop_count += NUM_DIRECTIONS*3;
#endif

	/* call the main routine */
	fit_ellipsoid(p, NUM_DIRECTIONS, Q);
}
