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

/*****************************************************************************
 *  Debug tools
 ****************************************************************************/

/*****************************************************************************
 *  Initial implementation
 ****************************************************************************/

/**
 * returns p^T Q p 
 * flop count = 12adds + 12mults
 *					  = 24
 */ 
static double _quadratic_form(const double p[3], const double (*Q)[3][3])
{
	/* helps a little bit */
	double qp0, qp1, qp2;

	qp0 = (*Q)[0][0]*p[0] + (*Q)[0][1]*p[1] + (*Q)[0][2]*p[2];
	qp1 = (*Q)[1][0]*p[0] + (*Q)[1][1]*p[1] + (*Q)[1][2]*p[2];
	qp2 = (*Q)[2][0]*p[0] + (*Q)[2][1]*p[1] + (*Q)[2][2]*p[2];

	double ret = qp0*p[0] + qp1*p[1] + qp2*p[2];

	return ret;
}

/**
 * returns \sum_i (p_i^T Q p_i - 1)^2
 * flop count = n*(C(_quadratic_Form)+1sub + 1add+1mult)
 *            = n*(24+3) = 27n
 */
static double _cost(const double (*p)[3], int n, const double (*Q)[3][3])
{
#if 0
	/* this improves w.r.t -O2 but worse w.r.t -O3.. */
	double ret0 = 0, ret1 = 0;
	int i;
	for (i=0; i<n-1; i+=2) { 
		double residual_i   = _quadratic_form(p[i], Q) - 1;
		double residual_ip1 = _quadratic_form(p[i+1], Q) - 1;

		ret0 += residual_i*residual_i;
		ret1 += residual_ip1*residual_ip1;
	}

	/* cleanup */
	ret0 += ret1;
	for (; i<n; i++) {
		double residual_i = _quadratic_form(p[i], Q) - 1;
		ret0 += residual_i*residual_i;
	}

	return ret0;
#else
	double ret = 0;

	for (int i=0; i<n; i++){
		double residual_i = _quadratic_form(p[i], Q) - 1;
		ret += residual_i*residual_i;
	}
	return ret;
#endif
}

/*
 * ||A||_F
 * flop count = 9*(1mult + 1add) + 1sqrt
 */
static double _normFroSq(const double A[3][3])
{
	double ret0, ret1, ret2;
	/* don't see a change in perf w.r.t -03 */
	ret0 = A[0][0]*A[0][0] + A[0][1]*A[0][1] + A[0][2]*A[0][2];
	ret1 = A[1][0]*A[1][0] + A[1][1]*A[1][1] + A[1][2]*A[1][2];
	ret2 = A[2][0]*A[2][0] + A[2][1]*A[2][1] + A[2][2]*A[2][2];
	return ret0+ret1+ret2;
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
void fit_ellipsoid_opt(const double (*p)[3], int n, double (*Q)[3][3])
{
	(*Q)[0][0] = 1; (*Q)[0][1] = 0; (*Q)[0][2] = 0;
	(*Q)[1][0] = 0; (*Q)[1][1] = 1; (*Q)[1][2] = 0;
	(*Q)[2][0] = 0; (*Q)[2][1] = 0; (*Q)[2][2] = 1;

	while(1) {
		/* determine gradient */


#if 0
		/* slightly worse w.r.t -O3 */
		
		double grad00=0, grad01=0, grad02=0, grad10=0, grad11=0, grad12=0, grad20=0, grad21=0, grad22=0;
		int i;
		for (i=0; i<n; i++) {
			double residual_i   = _quadratic_form(p[i], Q) - 1;

			grad00 += 2*residual_i*p[i][0]*p[i][0];
			grad01 += 2*residual_i*p[i][0]*p[i][1];
			grad02 += 2*residual_i*p[i][0]*p[i][2];
			grad10 += 2*residual_i*p[i][1]*p[i][0];
			grad11 += 2*residual_i*p[i][1]*p[i][1];
			grad12 += 2*residual_i*p[i][1]*p[i][2];
			grad20 += 2*residual_i*p[i][2]*p[i][0];
			grad21 += 2*residual_i*p[i][2]*p[i][1];
			grad22 += 2*residual_i*p[i][2]*p[i][2];
		}

		double ret = grad00*grad00 + grad01*grad01 + grad02*grad02 +
								 grad10*grad10 + grad11*grad11 + grad12*grad12 +
								 grad20*grad20 + grad21*grad21 + grad22*grad22;
		if (sqrt(ret) < EPSILON) break;

		/* take step direction to be negative gradient */
		
		// double step[3][3];
		// for (int i=0; i<3; i++) {
		// 	for (int j=0; j<3; j++) {
		// 		step[i][j] = -grad[i][j];
		// 	}
		// }

		/* backtracking line search */
		
		double t = 100;
		double Qk_plus_tstep[3][3] = {{0}};
		Qk_plus_tstep[0][0] = (*Q)[0][0] - t*grad00;
		Qk_plus_tstep[0][1] = (*Q)[0][1] - t*grad01;
		Qk_plus_tstep[0][2] = (*Q)[0][2] - t*grad02;
		Qk_plus_tstep[1][0] = (*Q)[1][0] - t*grad10;
		Qk_plus_tstep[1][1] = (*Q)[1][1] - t*grad11;
		Qk_plus_tstep[1][2] = (*Q)[1][2] - t*grad12;
		Qk_plus_tstep[2][0] = (*Q)[2][0] - t*grad20;
		Qk_plus_tstep[2][1] = (*Q)[2][1] - t*grad21;
		Qk_plus_tstep[2][2] = (*Q)[2][2] - t*grad22;
		
		double trace_gradstep = -ret;
		double _costQ         = _cost(p, n, Q);

		while (_cost(p, n, (double (*)[3][3]) Qk_plus_tstep) > _costQ + ALPHA*t*trace_gradstep) {
			t = t*BETA;

			Qk_plus_tstep[0][0] = (*Q)[0][0] - t*grad00;
			Qk_plus_tstep[0][1] = (*Q)[0][1] - t*grad01;
			Qk_plus_tstep[0][2] = (*Q)[0][2] - t*grad02;
			Qk_plus_tstep[1][0] = (*Q)[1][0] - t*grad10;
			Qk_plus_tstep[1][1] = (*Q)[1][1] - t*grad11;
			Qk_plus_tstep[1][2] = (*Q)[1][2] - t*grad12;
			Qk_plus_tstep[2][0] = (*Q)[2][0] - t*grad20;
			Qk_plus_tstep[2][1] = (*Q)[2][1] - t*grad21;
			Qk_plus_tstep[2][2] = (*Q)[2][2] - t*grad22;
		} /* back tracking while loop */
#else
		/* determine gradient */

		double grad[3][3] = {{0}};

		for (int i=0; i<n; i++) {
			double residual_i = _quadratic_form(p[i], Q) - 1;
			double jacobian_i[3][3] = {{0}};
		
			/* jacobian_i =  outer product */
			for (int j=0; j<3; j++) {
				for (int k=0; k<3; k++) {
					jacobian_i[j][k] = p[i][j]*p[i][k];
					grad[j][k] += 2*residual_i*jacobian_i[j][k];
				}
			}
		}

		if (_normFroSq(grad) < EPSILON*EPSILON) break;

		/* take step direction to be negative gradient */
		
		/* backtracking line search */
		
		double t = 100;
		double Qk_plus_tstep[3][3] = {{0}};
		for (int i=0; i<3; i++) {
			for (int j=0; j<3; j++) {
				Qk_plus_tstep[i][j] = (*Q)[i][j] - t*grad[i][j];
			}
		}
		double trace_gradstep = -_normFroSq(grad);
		double _costQ         = _cost(p, n, Q);

		while (_cost(p, n, (double (*)[3][3]) Qk_plus_tstep) > 
										_costQ + ALPHA*t*trace_gradstep) {
			t = t*BETA;
			
			for (int i=0; i<3; i++) {
				for (int j=0; j<3; j++) {
					Qk_plus_tstep[i][j] = (*Q)[i][j] - t*grad[i][j];
				}
			}
		} /* back tracking while loop */
#endif
		memcpy(Q, Qk_plus_tstep, sizeof(double)*9);
	} /* main while loop */
}