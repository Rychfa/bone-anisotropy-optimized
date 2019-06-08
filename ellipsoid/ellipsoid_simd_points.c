#include <string.h>
#include <math.h>
#include "ellipsoid.h"
#ifdef DEBUG
	#include <stdio.h>
#endif

#include <immintrin.h>
#include <stdlib.h>

/* some parameters */
static const double EPSILON = 1e-6;
static const double ALPHA   = 0.25;
static const double BETA    = 0.5;

/*****************************************************************************
 *  Preprocessor 
 ****************************************************************************/

#define VAR(var,i) var##_##i

#define RESIDUAL_ACC(k,i) \
	__m256d VAR(p0,k), VAR(p1,k), VAR(p2,k);\
	\
	VAR(p0,k) = _mm256_load_pd(&p0[i]);\
	VAR(p1,k) = _mm256_load_pd(&p1[i]);\
	VAR(p2,k) = _mm256_load_pd(&p2[i]);\
	\
	__m256d VAR(p0p0,k), VAR(p0p1,k), VAR(p0p2,k), VAR(p1p1,k), VAR(p1p2,k), VAR(p2p2,k);\
	VAR(p0p0,k) = _mm256_mul_pd(VAR(p0,k), VAR(p0,k));\
	VAR(p0p1,k) = _mm256_mul_pd(VAR(p0,k), VAR(p1,k));\
	VAR(p0p2,k) = _mm256_mul_pd(VAR(p0,k), VAR(p2,k));\
	VAR(p1p1,k) = _mm256_mul_pd(VAR(p1,k), VAR(p1,k));\
	VAR(p1p2,k) = _mm256_mul_pd(VAR(p1,k), VAR(p2,k));\
	VAR(p2p2,k) = _mm256_mul_pd(VAR(p2,k), VAR(p2,k));\
	\
	__m256d VAR(t0,k), VAR(t1,k);\
	VAR(t0,k) = _mm256_mul_pd(VAR(p0p0,k), Q00);\
	VAR(t1,k) = _mm256_mul_pd(VAR(p2p2,k), Q22);\
	VAR(t0,k) = _mm256_fmadd_pd(VAR(p0p1,k), Q01, VAR(t0,k));\
	VAR(t0,k) = _mm256_fmadd_pd(VAR(p0p2,k), Q02, VAR(t0,k));\
	VAR(t0,k) = _mm256_fmadd_pd(VAR(p0p1,k), Q01, VAR(t0,k));\
	VAR(t1,k) = _mm256_fmadd_pd(VAR(p1p2,k), Q12, VAR(t1,k));\
	VAR(t1,k) = _mm256_fmadd_pd(VAR(p0p2,k), Q02, VAR(t1,k));\
	VAR(t1,k) = _mm256_fmadd_pd(VAR(p1p2,k), Q12, VAR(t1,k));\
	VAR(t0,k) = _mm256_add_pd(VAR(t0,k), VAR(t1,k));\
	VAR(t1,k) = _mm256_fmadd_pd(VAR(p1p1,k), Q11, VAR(t0,k));\
	\
	__m256d VAR(residual,k) = _mm256_sub_pd(VAR(t1,k), vONES);\
	VAR(costQ,k) = _mm256_fmadd_pd(VAR(residual,k), VAR(residual,k), VAR(costQ,k));

/*****************************************************************************
 *  Private data
 ****************************************************************************/

static __m256d vONES, vBETA, vTWOS;
static __m128i vINDEX;

/*****************************************************************************
 *  Debug tools
 ****************************************************************************/

/*****************************************************************************
 *  Initial implementation
 ****************************************************************************/
/**
 * returns \sum_i (p_i^T Q p_i - 1)^2
 * flop count = n*(C(_quadratic_Form)+1sub + 1add+1mult)
 *            = n*(24+3) = 27n
 */
#if 1
static double _cost(const double *p0, const double *p1, const double *p2, int n,
		__m256d Q00, __m256d Q01, __m256d Q02, __m256d Q11, __m256d Q12, __m256d Q22)
{
	__m256d costQ_0 = _mm256_setzero_pd();

	int i;
	for (i=0; i<n; i+=4) {
		RESIDUAL_ACC(0, i);
	}

	return costQ_0[0] + costQ_0[1] + costQ_0[2] + costQ_0[3];
}

#else
/* unrolling cost doesn't help */
static double _cost(const double *p0, const double *p1, const double *p2, int n,
		__m256d Q00, __m256d Q01, __m256d Q02, __m256d Q11, __m256d Q12, __m256d Q22)
{
	__m256d costQ_0 = _mm256_setzero_pd();
	__m256d costQ_1 = _mm256_setzero_pd();

	int i;
	for (i=0; i<n; i+=4*2) {
		RESIDUAL_ACC(0, i);
		RESIDUAL_ACC(1, i+4);
	}
	costQ_0 = _mm256_add_pd(costQ_0, costQ_1);
	return costQ_0[0] + costQ_0[1] + costQ_0[2] + costQ_0[3];
}
#endif

void fit_ellipsoid_simd_points_init(void)
{
	vBETA  = _mm256_set1_pd(BETA);
	vONES  = _mm256_set1_pd(1);
	vTWOS  = _mm256_set1_pd(2);	
	vINDEX = _mm_setr_epi32(0, 3, 6, 9);
}

#if 1
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

void fit_ellipsoid_simd_points(const double (*p)[3], int n, double (*Q)[3][3])
{
	__m256d Q00, Q01, Q02, Q11, Q12, Q22;
	Q00 = _mm256_set1_pd(1);
	Q01 = _mm256_set1_pd(0);
	Q02 = _mm256_set1_pd(0);
	Q11 = _mm256_set1_pd(1);
	Q12 = _mm256_set1_pd(0);
	Q22 = _mm256_set1_pd(1);

	/* shuffle variables once and store */
	double p0[MAX_ARRAY_DIM] = {0};
	double p1[MAX_ARRAY_DIM] = {0};
	double p2[MAX_ARRAY_DIM] = {0};
	int i;
	for (i=0; i<n-4; i+=4) {
		__m256d p0_0, p1_0, p2_0;
		p0_0 = _mm256_i32gather_pd(&p[i][0], vINDEX, sizeof(double));
		p1_0 = _mm256_i32gather_pd(&p[i][1], vINDEX, sizeof(double));
		p2_0 = _mm256_i32gather_pd(&p[i][2], vINDEX, sizeof(double));
		
		_mm256_store_pd(&p0[i], p0_0);
		_mm256_store_pd(&p1[i], p1_0);
		_mm256_store_pd(&p2[i], p2_0);
	}
	for (;i<n; i++) {
		p0[i] = p[i][0];
		p1[i] = p[i][1];
		p2[i] = p[i][2];
	}

	while(1) {
		/* determine gradient */

		__m256d g00, g01, g02, g11, g12, g22;
		g00             = _mm256_setzero_pd();
		g01             = _mm256_setzero_pd();
		g02             = _mm256_setzero_pd();
		g11             = _mm256_setzero_pd();
		g12             = _mm256_setzero_pd();
		g22             = _mm256_setzero_pd();
		__m256d costQ_0 = _mm256_setzero_pd();
		
		for (i=0; i<n; i+=4) {
			RESIDUAL_ACC(0, i);
			__m256d t0;
			t0 = _mm256_mul_pd(residual_0, vTWOS); 

			g00 = _mm256_fmadd_pd(t0, p0p0_0, g00);
			g01 = _mm256_fmadd_pd(t0, p0p1_0, g01);
			g02 = _mm256_fmadd_pd(t0, p0p2_0, g02);
			g11 = _mm256_fmadd_pd(t0, p1p1_0, g11);
			g12 = _mm256_fmadd_pd(t0, p1p2_0, g12);
			g22 = _mm256_fmadd_pd(t0, p2p2_0, g22);
		}

		costQ_0[0] = costQ_0[0] + costQ_0[1] + costQ_0[2] + costQ_0[3];

		g00 = _mm256_set1_pd(g00[0]+g00[1]+g00[2]+g00[3]);
		g01 = _mm256_set1_pd(g01[0]+g01[1]+g01[2]+g01[3]);
		g02 = _mm256_set1_pd(g02[0]+g02[1]+g02[2]+g02[3]);
		g11 = _mm256_set1_pd(g11[0]+g11[1]+g11[2]+g11[3]);
		g12 = _mm256_set1_pd(g12[0]+g12[1]+g12[2]+g12[3]);
		g22 = _mm256_set1_pd(g22[0]+g22[1]+g22[2]+g22[3]);

		
		double normFro_sq = g00[0]*g00[0] + 2*g01[0]*g01[0] + 2*g02[0]*g02[0] + g11[0]*g11[0] + 2*g12[0]*g12[0] + g22[0]*g22[0];
		
		if (normFro_sq < EPSILON*EPSILON) break;

		/* take step direction to be negative gradient */
		
		/* backtracking line search */
		
		__m256d Qnext00, Qnext01, Qnext02, Qnext11, Qnext12, Qnext22;
		__m256d nt   = _mm256_set1_pd(-100);

		Qnext00 = _mm256_fmadd_pd(nt, g00, Q00);
		Qnext01 = _mm256_fmadd_pd(nt, g01, Q01);
		Qnext02 = _mm256_fmadd_pd(nt, g02, Q02);
		Qnext11 = _mm256_fmadd_pd(nt, g11, Q11);
		Qnext12 = _mm256_fmadd_pd(nt, g12, Q12);
		Qnext22 = _mm256_fmadd_pd(nt, g22, Q22);


		while (_cost(p0, p1, p2, n, Qnext00, Qnext01, Qnext02, Qnext11, Qnext12, Qnext22) > 
																					costQ_0[0] + ALPHA*nt[0]*normFro_sq) {
			nt = _mm256_mul_pd(nt, vBETA);

			Qnext00 = _mm256_fmadd_pd(nt, g00, Q00);
			Qnext01 = _mm256_fmadd_pd(nt, g01, Q01);
			Qnext02 = _mm256_fmadd_pd(nt, g02, Q02);
			Qnext11 = _mm256_fmadd_pd(nt, g11, Q11);
			Qnext12 = _mm256_fmadd_pd(nt, g12, Q12);
			Qnext22 = _mm256_fmadd_pd(nt, g22, Q22);
		} /* back tracking while loop */

		Q00 = Qnext00;
		Q01 = Qnext01;
		Q02 = Qnext02;
		Q11 = Qnext11;
		Q12 = Qnext12;
		Q22 = Qnext22;
	} /* main while loop */
	(*Q)[0][0] = Q00[0];
	(*Q)[0][1] = Q01[0];
	(*Q)[0][2] = Q02[0];
	(*Q)[1][0] = Q01[0];
	(*Q)[1][1] = Q11[0];
	(*Q)[1][2] = Q12[0];
	(*Q)[2][0] = Q02[0];
	(*Q)[2][1] = Q12[0];
	(*Q)[2][2] = Q22[0];
}
#else
/*****************************************************************************
 *  Unrolled: doesn't help
 ****************************************************************************/

void fit_ellipsoid_simd_points(const double (*p)[3], int n, double (*Q)[3][3])
{
	__m256d Q00, Q01, Q02, Q11, Q12, Q22;
	Q00 = _mm256_set1_pd(1);
	Q01 = _mm256_set1_pd(0);
	Q02 = _mm256_set1_pd(0);
	Q11 = _mm256_set1_pd(1);
	Q12 = _mm256_set1_pd(0);
	Q22 = _mm256_set1_pd(1);

	/* shuffle variables once and store */
	static double p0[MAX_ARRAY_DIM] = {0};
	static double p1[MAX_ARRAY_DIM] = {0};
	static double p2[MAX_ARRAY_DIM] = {0};
	
	__m128i vindex = _mm_setr_epi32(0, 3, 6, 9);
	for (int i=0; i<n; i+=4) {
		__m256d p0_0, p1_0, p2_0;
		p0_0 = _mm256_i32gather_pd(&p[i][0], vindex, sizeof(double));
		p1_0 = _mm256_i32gather_pd(&p[i][1], vindex, sizeof(double));
		p2_0 = _mm256_i32gather_pd(&p[i][2], vindex, sizeof(double));
		
		_mm256_store_pd(&p0[i], p0_0);
		_mm256_store_pd(&p1[i], p1_0);
		_mm256_store_pd(&p2[i], p2_0);
	}

	while(1) {
		/* determine gradient */

		__m256d g00, g01, g02, g11, g12, g22;
		g00             = _mm256_setzero_pd();
		g01             = _mm256_setzero_pd();
		g02             = _mm256_setzero_pd();
		g11             = _mm256_setzero_pd();
		g12             = _mm256_setzero_pd();
		g22             = _mm256_setzero_pd();
		__m256d costQ_0 = _mm256_setzero_pd();

		__m256d g00_1, g01_1, g02_1, g11_1, g12_1, g22_1;
		g00_1             = _mm256_setzero_pd();
		g01_1             = _mm256_setzero_pd();
		g02_1             = _mm256_setzero_pd();
		g11_1             = _mm256_setzero_pd();
		g12_1             = _mm256_setzero_pd();
		g22_1             = _mm256_setzero_pd();
		__m256d costQ_1 = _mm256_setzero_pd();
		
		int i;
		for (i=0; i<n; i+=4*2) {
			RESIDUAL_ACC(0, i);
			RESIDUAL_ACC(1, i+4);

			__m256d t0;
			t0 = _mm256_mul_pd(residual_0, vTWOS); 

			g00 = _mm256_fmadd_pd(t0, p0p0_0, g00);
			g01 = _mm256_fmadd_pd(t0, p0p1_0, g01);
			g02 = _mm256_fmadd_pd(t0, p0p2_0, g02);
			g11 = _mm256_fmadd_pd(t0, p1p1_0, g11);
			g12 = _mm256_fmadd_pd(t0, p1p2_0, g12);
			g22 = _mm256_fmadd_pd(t0, p2p2_0, g22);

			__m256d t1;
			t1 = _mm256_mul_pd(residual_1, vTWOS); 

			g00_1 = _mm256_fmadd_pd(t1, p0p0_1, g00_1);
			g01_1 = _mm256_fmadd_pd(t1, p0p1_1, g01_1);
			g02_1 = _mm256_fmadd_pd(t1, p0p2_1, g02_1);
			g11_1 = _mm256_fmadd_pd(t1, p1p1_1, g11_1);
			g12_1 = _mm256_fmadd_pd(t1, p1p2_1, g12_1);
			g22_1 = _mm256_fmadd_pd(t1, p2p2_1, g22_1);
		}

		costQ_0 = _mm256_add_pd(costQ_0, costQ_1);
		costQ_0[0] = costQ_0[0] + costQ_0[1] + costQ_0[2] + costQ_0[3];

		g00 = _mm256_add_pd(g00, g00_1);
		g01 = _mm256_add_pd(g01, g01_1);
		g02 = _mm256_add_pd(g02, g02_1);
		g11 = _mm256_add_pd(g11, g11_1);
		g12 = _mm256_add_pd(g12, g12_1);
		g22 = _mm256_add_pd(g22, g22_1);

		g00 = _mm256_set1_pd(g00[0]+g00[1]+g00[2]+g00[3]);
		g01 = _mm256_set1_pd(g01[0]+g01[1]+g01[2]+g01[3]);
		g02 = _mm256_set1_pd(g02[0]+g02[1]+g02[2]+g02[3]);
		g11 = _mm256_set1_pd(g11[0]+g11[1]+g11[2]+g11[3]);
		g12 = _mm256_set1_pd(g12[0]+g12[1]+g12[2]+g12[3]);
		g22 = _mm256_set1_pd(g22[0]+g22[1]+g22[2]+g22[3]);

		
		double normFro_sq = g00[0]*g00[0] + 2*g01[0]*g01[0] + 2*g02[0]*g02[0] + g11[0]*g11[0] + 2*g12[0]*g12[0] + g22[0]*g22[0];
		
		if (normFro_sq < EPSILON*EPSILON) break;

		/* take step direction to be negative gradient */
		
		/* backtracking line search */
		
		__m256d Qnext00, Qnext01, Qnext02, Qnext11, Qnext12, Qnext22;
		__m256d nt   = _mm256_set1_pd(-100);

		Qnext00 = _mm256_fmadd_pd(nt, g00, Q00);
		Qnext01 = _mm256_fmadd_pd(nt, g01, Q01);
		Qnext02 = _mm256_fmadd_pd(nt, g02, Q02);
		Qnext11 = _mm256_fmadd_pd(nt, g11, Q11);
		Qnext12 = _mm256_fmadd_pd(nt, g12, Q12);
		Qnext22 = _mm256_fmadd_pd(nt, g22, Q22);


		while (_cost(p0, p1, p2, n, Qnext00, Qnext01, Qnext02, Qnext11, Qnext12, Qnext22) > 
																					costQ_0[0] + ALPHA*nt[0]*normFro_sq) {
			nt = _mm256_mul_pd(nt, vBETA);

			Qnext00 = _mm256_fmadd_pd(nt, g00, Q00);
			Qnext01 = _mm256_fmadd_pd(nt, g01, Q01);
			Qnext02 = _mm256_fmadd_pd(nt, g02, Q02);
			Qnext11 = _mm256_fmadd_pd(nt, g11, Q11);
			Qnext12 = _mm256_fmadd_pd(nt, g12, Q12);
			Qnext22 = _mm256_fmadd_pd(nt, g22, Q22);
		} /* back tracking while loop */

		Q00 = Qnext00;
		Q01 = Qnext01;
		Q02 = Qnext02;
		Q11 = Qnext11;
		Q12 = Qnext12;
		Q22 = Qnext22;
	} /* main while loop */
	(*Q)[0][0] = Q00[0];
	(*Q)[0][1] = Q01[0];
	(*Q)[0][2] = Q02[0];
	(*Q)[1][0] = Q01[0];
	(*Q)[1][1] = Q11[0];
	(*Q)[1][2] = Q12[0];
	(*Q)[2][0] = Q02[0];
	(*Q)[2][1] = Q12[0];
	(*Q)[2][2] = Q22[0];
}

#endif


// void fit_ellipsoid_mils_simd(const double *mils, double (*Q)[3][3])
// {
//        /* construct the points */
//        double p[NUM_DIRECTIONS][3];
//        for (int i=0; i<NUM_DIRECTIONS; i++) {
//                for (int j=0; j<3; j++) {
//                        p[i][j] = mils[i] * DIRECTIONS_NORMALIZED[i][j];
//                }
//        }
//        /* call the main routine */
//        fit_ellipsoid_simd_points(p, NUM_DIRECTIONS, Q);
// }

void fit_ellipsoid_mils_simd(const double *mils, double (*Q)[3][3])
{
	__m256d Q00, Q01, Q02, Q11, Q12, Q22;
	Q00 = _mm256_set1_pd(1);
	Q01 = _mm256_set1_pd(0);
	Q02 = _mm256_set1_pd(0);
	Q11 = _mm256_set1_pd(1);
	Q12 = _mm256_set1_pd(0);
	Q22 = _mm256_set1_pd(1);

	/* shuffle variables once and store */
	double p0[MILS_ARRAY_DIM] = {0};
	double p1[MILS_ARRAY_DIM] = {0};
	double p2[MILS_ARRAY_DIM] = {0};
	int i;
	for (i=0; i<NUM_DIRECTIONS; i++) {
		p0[i] = mils[i]*DIRECTIONS_NORMALIZED[i][0];
		p1[i] = mils[i]*DIRECTIONS_NORMALIZED[i][1];
		p2[i] = mils[i]*DIRECTIONS_NORMALIZED[i][2];
	}

	while(1) {
		/* determine gradient */

		__m256d g00, g01, g02, g11, g12, g22;
		g00             = _mm256_setzero_pd();
		g01             = _mm256_setzero_pd();
		g02             = _mm256_setzero_pd();
		g11             = _mm256_setzero_pd();
		g12             = _mm256_setzero_pd();
		g22             = _mm256_setzero_pd();
		__m256d costQ_0 = _mm256_setzero_pd();
		
		for (i=0; i<NUM_DIRECTIONS; i+=4) {
			RESIDUAL_ACC(0, i);
			__m256d t0;
			t0 = _mm256_mul_pd(residual_0, vTWOS); 

			g00 = _mm256_fmadd_pd(t0, p0p0_0, g00);
			g01 = _mm256_fmadd_pd(t0, p0p1_0, g01);
			g02 = _mm256_fmadd_pd(t0, p0p2_0, g02);
			g11 = _mm256_fmadd_pd(t0, p1p1_0, g11);
			g12 = _mm256_fmadd_pd(t0, p1p2_0, g12);
			g22 = _mm256_fmadd_pd(t0, p2p2_0, g22);
		}

		costQ_0[0] = costQ_0[0] + costQ_0[1] + costQ_0[2] + costQ_0[3];

		g00 = _mm256_set1_pd(g00[0]+g00[1]+g00[2]+g00[3]);
		g01 = _mm256_set1_pd(g01[0]+g01[1]+g01[2]+g01[3]);
		g02 = _mm256_set1_pd(g02[0]+g02[1]+g02[2]+g02[3]);
		g11 = _mm256_set1_pd(g11[0]+g11[1]+g11[2]+g11[3]);
		g12 = _mm256_set1_pd(g12[0]+g12[1]+g12[2]+g12[3]);
		g22 = _mm256_set1_pd(g22[0]+g22[1]+g22[2]+g22[3]);

		
		double normFro_sq = g00[0]*g00[0] + 2*g01[0]*g01[0] + 2*g02[0]*g02[0] + g11[0]*g11[0] + 2*g12[0]*g12[0] + g22[0]*g22[0];
		
		if (normFro_sq < EPSILON*EPSILON) break;

		/* take step direction to be negative gradient */
		
		/* backtracking line search */
		
		__m256d Qnext00, Qnext01, Qnext02, Qnext11, Qnext12, Qnext22;
		__m256d nt   = _mm256_set1_pd(-100);

		Qnext00 = _mm256_fmadd_pd(nt, g00, Q00);
		Qnext01 = _mm256_fmadd_pd(nt, g01, Q01);
		Qnext02 = _mm256_fmadd_pd(nt, g02, Q02);
		Qnext11 = _mm256_fmadd_pd(nt, g11, Q11);
		Qnext12 = _mm256_fmadd_pd(nt, g12, Q12);
		Qnext22 = _mm256_fmadd_pd(nt, g22, Q22);


		while (_cost(p0, p1, p2, NUM_DIRECTIONS, Qnext00, Qnext01, Qnext02, Qnext11, Qnext12, Qnext22) > 
																					costQ_0[0] + ALPHA*nt[0]*normFro_sq) {
			nt = _mm256_mul_pd(nt, vBETA);

			Qnext00 = _mm256_fmadd_pd(nt, g00, Q00);
			Qnext01 = _mm256_fmadd_pd(nt, g01, Q01);
			Qnext02 = _mm256_fmadd_pd(nt, g02, Q02);
			Qnext11 = _mm256_fmadd_pd(nt, g11, Q11);
			Qnext12 = _mm256_fmadd_pd(nt, g12, Q12);
			Qnext22 = _mm256_fmadd_pd(nt, g22, Q22);
		} /* back tracking while loop */

		Q00 = Qnext00;
		Q01 = Qnext01;
		Q02 = Qnext02;
		Q11 = Qnext11;
		Q12 = Qnext12;
		Q22 = Qnext22;
	} /* main while loop */
	(*Q)[0][0] = Q00[0];
	(*Q)[0][1] = Q01[0];
	(*Q)[0][2] = Q02[0];
	(*Q)[1][0] = Q01[0];
	(*Q)[1][1] = Q11[0];
	(*Q)[1][2] = Q12[0];
	(*Q)[2][0] = Q02[0];
	(*Q)[2][1] = Q12[0];
	(*Q)[2][2] = Q22[0];
}