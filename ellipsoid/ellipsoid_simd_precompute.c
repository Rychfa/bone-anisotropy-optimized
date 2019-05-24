#include <string.h>
#include <math.h>
#include "ellipsoid.h"
#ifdef DEBUG
	#include <stdio.h>
#endif

#include <immintrin.h>
#include <stdlib.h>

/* some parameters */
static const double EPSILON = 1e-3;
static const double ALPHA   = 0.25;
static const double BETA    = 0.5;

/*****************************************************************************
 *  Preprocessor 
 ****************************************************************************/

#define VAR(var,i) var##_##i

/* assumes: pPtr, mask, costQ_i
 * declares: p_i, p0_i, p1_i, p2_i, t0_i, t1_i, t2_i, t3_i 
 */
#if 1
#define RESIDUAL_ACC(k,i) \
	__m256d VAR(ppT0,k), VAR(ppT1,k), VAR(ppT2,k);\
	\
	VAR(ppT0,k) = _mm256_load_pd(&ppT0[i+k][0]);\
	VAR(ppT1,k) = _mm256_load_pd(&ppT1[i+k][0]);\
	VAR(ppT2,k) = _mm256_load_pd(&ppT2[i+k][0]);\
	\
	__m256d VAR(QppT0,k);\
	VAR(QppT0,k) = _mm256_mul_pd(VAR(ppT0,k), Q0);\
	\
	__m256d VAR(t0,k), VAR(t1,k);\
	VAR(t0,k) = _mm256_fmadd_pd(VAR(ppT1,k), Q1, VAR(QppT0,k));\
	VAR(t1,k) = _mm256_fmadd_pd(VAR(ppT2,k), Q2, VAR(t0,k));\
	\
	double VAR(residual,k) = VAR(t1,k)[0] + VAR(t1,k)[1] + VAR(t1,k)[2] - 1;\
	VAR(costQ,k) += VAR(residual,k)*VAR(residual,k);

#else
	/* these two about the same */
#define RESIDUAL_ACC(k,i) \
	__m256d VAR(ppT0,k), VAR(ppT1,k), VAR(ppT2,k);\
	\
	VAR(ppT0,k) = _mm256_load_pd(&ppT0[i+k][0]);\
	VAR(ppT1,k) = _mm256_load_pd(&ppT1[i+k][0]);\
	VAR(ppT2,k) = _mm256_load_pd(&ppT2[i+k][0]);\
	\
	__m256d VAR(QppT0,k);\
	VAR(QppT0,k) = _mm256_mul_pd(VAR(ppT0,k), Q0);\
	\
	__m256d VAR(t0,k), VAR(t1,k);\
	VAR(t0,k) = _mm256_fmadd_pd(VAR(ppT1,k), Q1, VAR(QppT0,k));\
	VAR(t1,k) = _mm256_fmadd_pd(VAR(ppT2,k), Q2, VAR(t0,k));\
	\
	__m256d VAR(t8,k), VAR(t9,k), VAR(t10,k);\
	VAR(t8,k) = _mm256_permute4x64_pd(VAR(t1,k), 0b00000000);\
	VAR(t9,k) = _mm256_permute4x64_pd(VAR(t1,k), 0b01010101);\
	VAR(t10,k) = _mm256_permute4x64_pd(VAR(t1,k), 0b10101010);\
	VAR(t8,k) = _mm256_add_pd(VAR(t8,k), VAR(t9,k));\
	VAR(t8,k) = _mm256_add_pd(VAR(t8,k), VAR(t10,k));\
	\
	__m256d VAR(residual,k);\
	VAR(residual,k) = _mm256_sub_pd(VAR(t8,k), vONES);\
	VAR(costQ,k) += VAR(residual,k)[0]*VAR(residual,k)[0];
#endif

/*****************************************************************************
 *  Private data
 ****************************************************************************/

static __m256d vONES, vBETA, vTWOS;
static __m256i vMASK;

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
static double _cost(const double (*ppT0)[4], const double (*ppT1)[4],const double (*ppT2)[4], int n,
		__m256d Q0, __m256d Q1, __m256d Q2)
{
	double costQ_0 = 0;

	int i;
	for (i=0; i<n; i++) {
		RESIDUAL_ACC(0, i);
	}

	return costQ_0;
}

#else
static double _cost(const double (*ppT0)[4], const double (*ppT1)[4],const double (*ppT2)[4], int n,
		__m256d Q0, __m256d Q1, __m256d Q2)
{
	double costQ_0 = 0;
	double costQ_1 = 0;

	int i;
	for (i=0; i<n-2; i+=2) { /* unrolling twice: helps for large n > ~100 */
		RESIDUAL_ACC(0, i);
		RESIDUAL_ACC(1, i);
	}

	/* cleanup */
	for (; i<n; i++) {
		RESIDUAL_ACC(0, i);
	}

	return costQ_0 + costQ_1;
}
#endif

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

void fit_ellipsoid_simd_precompute_init(void)
{
	vMASK = _mm256_setr_epi64x(0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0);
	vBETA = _mm256_set1_pd(BETA);
	vONES = _mm256_set1_pd(1);
	vTWOS = _mm256_set1_pd(2);
}

void fit_ellipsoid_simd_precompute(const double (*p)[3], int n, double (*Q)[3][3])
{
	__m256d Q0, Q1, Q2;
	Q0 = _mm256_setr_pd(1, 0, 0, 0);
	Q1 = _mm256_setr_pd(0, 1, 0, 0);
	Q2 = _mm256_setr_pd(0, 0, 1, 0);

	/* pre compute reused variables and put on stack. note this modifies the 
	 * flop count!!
	 */
	static double ppT0[MAX_ARRAY_DIM][4] = {{0}};
	static double ppT1[MAX_ARRAY_DIM][4] = {{0}};
	static double ppT2[MAX_ARRAY_DIM][4] = {{0}};
	
	for (int i=0; i<n; i++) {
		__m256d p_i, p0, p1, p2;
		p_i = _mm256_maskload_pd(&p[i][0], vMASK);
		p0  = _mm256_permute4x64_pd(p_i, 0b00000000);
		p1  = _mm256_permute4x64_pd(p_i, 0b01010101);
		p2  = _mm256_permute4x64_pd(p_i, 0b10101010);
		// p0  = _mm256_broadcast_sd(&p[i][0]);
		// p1  = _mm256_broadcast_sd(&p[i][1]);
		// p2  = _mm256_broadcast_sd(&p[i][2]);

		__m256d ppT0_i, ppT1_i, ppT2_i;
		ppT0_i = _mm256_mul_pd(p_i, p0);
		ppT1_i = _mm256_mul_pd(p_i, p1);
		ppT2_i = _mm256_mul_pd(p_i, p2);

		_mm256_store_pd(&ppT0[i][0], ppT0_i);
		_mm256_store_pd(&ppT1[i][0], ppT1_i);
		_mm256_store_pd(&ppT2[i][0], ppT2_i);
	}

	while(1) {
		/* determine gradient */

		__m256d g0, g1, g2;
		g0 = _mm256_setzero_pd();
		g1 = _mm256_setzero_pd();
		g2 = _mm256_setzero_pd();

		double costQ_0 = 0;

		int i;
		for (i=0; i<n; i++) {
			RESIDUAL_ACC(0, i);

			__m256d t0;
			// t0 = _mm256_mul_pd(residual_0, vTWOS); 
			t0 = _mm256_set1_pd(residual_0*2);

			g0 = _mm256_fmadd_pd(t0, ppT0_0, g0);
			g1 = _mm256_fmadd_pd(t0, ppT1_0, g1);
			g2 = _mm256_fmadd_pd(t0, ppT2_0, g2);
		}

		__m256d g02, g12, g22;
		g02 = _mm256_mul_pd(g0, g0);
		g12 = _mm256_mul_pd(g1, g1);
		g22 = _mm256_mul_pd(g2, g2);

		g02 = _mm256_add_pd(g02, g12);
		g02 = _mm256_add_pd(g02, g22);

		double normFro_sq = g02[0] + g02[1] + g02[2];
		if (normFro_sq < EPSILON*EPSILON) break;

		/* take step direction to be negative gradient */
		
		/* backtracking line search */
		
		__m256d Qnext0, Qnext1, Qnext2;
		__m256d nt   = _mm256_set1_pd(-100);

		Qnext0 = _mm256_fmadd_pd(nt, g0, Q0);
		Qnext1 = _mm256_fmadd_pd(nt, g1, Q1);
		Qnext2 = _mm256_fmadd_pd(nt, g2, Q2);

		while (_cost(ppT0, ppT1, ppT2, n, Qnext0, Qnext1, Qnext2) > 
																					costQ_0 + ALPHA*nt[0]*normFro_sq) {
			nt = _mm256_mul_pd(nt, vBETA);

			Qnext0 = _mm256_fmadd_pd(nt, g0, Q0);	
			Qnext1 = _mm256_fmadd_pd(nt, g1, Q1);	
			Qnext2 = _mm256_fmadd_pd(nt, g2, Q2);

		} /* back tracking while loop */

		Q0 = Qnext0;
		Q1 = Qnext1;
		Q2 = Qnext2;

	} /* main while loop */
	_mm256_maskstore_pd(&(*Q)[0][0], vMASK, Q0);
	_mm256_maskstore_pd(&(*Q)[1][0], vMASK, Q1);
	_mm256_maskstore_pd(&(*Q)[2][0], vMASK, Q2);
}
#else
/*****************************************************************************
 *  Unrolled
 ****************************************************************************/

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
void fit_ellipsoid_simd(const double (*pPtr)[3], int n, double (*Q)[3][3])
{
	__m256d Q0, Q1, Q2;
	Q0 = _mm256_setr_pd(1, 0, 0, 0);
	Q1 = _mm256_setr_pd(0, 1, 0, 0);
	Q2 = _mm256_setr_pd(0, 0, 1, 0);
	
	vMASK = _mm256_setr_epi64x(0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0);
	vBETA = _mm256_set1_pd(BETA);
	vONES = _mm256_set1_pd(1);
	vTWOS = _mm256_set1_pd(2);

	while(1) {
		/* determine gradient */

		__m256d g0, g1, g2, g0_1, g1_1, g2_1;
		g0 = _mm256_setzero_pd();
		g1 = _mm256_setzero_pd();
		g2 = _mm256_setzero_pd();

		g0_1 = _mm256_setzero_pd();
		g1_1 = _mm256_setzero_pd();
		g2_1 = _mm256_setzero_pd();

		/* can comptue cost at same time. note: flops probably slightly reduce */
		double costQ_0 = 0, costQ_1 = 0;

		int i;
		for (i=0; i<n; i+=2) {
			RESIDUAL_ACC(0, i);
			RESIDUAL_ACC(1, i);

			__m256d t0, t1;
			t0 = _mm256_mul_pd(residual_0, vTWOS); 
			t1 = _mm256_mul_pd(residual_1, vTWOS); 
			// t0 = _mm256_set1_pd(residual_0*2);
			// t1 = _mm256_set1_pd(residual_1*2);

			__m256d t3, t4, t5, t6, t7, t8;
			t3 = _mm256_mul_pd(p0_0, p_0); 
			t4 = _mm256_mul_pd(p1_0, p_0);
			t5 = _mm256_mul_pd(p2_0, p_0);
			
			t6 = _mm256_mul_pd(p0_1, p_1); 
			t7 = _mm256_mul_pd(p1_1, p_1);
			t8 = _mm256_mul_pd(p2_1, p_1);
			
			g0   = _mm256_fmadd_pd(t0, t3, g0);
			g1   = _mm256_fmadd_pd(t0, t4, g1);
			g2   = _mm256_fmadd_pd(t0, t5, g2);
			g0_1 = _mm256_fmadd_pd(t1, t6, g0_1);
			g1_1 = _mm256_fmadd_pd(t1, t7, g1_1);
			g2_1 = _mm256_fmadd_pd(t1, t8, g2_1);
		}
		/* cleanup */
		costQ_0 += costQ_1;
		g0 = _mm256_add_pd(g0, g0_1);
		g1 = _mm256_add_pd(g1, g1_1);
		g2 = _mm256_add_pd(g2, g2_1);

		/* compute frobenius norm sq */
		__m256d g02, g12, g22;
		g02 = _mm256_mul_pd(g0, g0);
		g12 = _mm256_mul_pd(g1, g1);
		g22 = _mm256_mul_pd(g2, g2);

		g02 = _mm256_add_pd(g02, g12);
		g02 = _mm256_add_pd(g02, g22);

		double normFro_sq = g02[0] + g02[1] + g02[2];
		if (normFro_sq < EPSILON*EPSILON) break;

		/* take step direction to be negative gradient */
		
		/* backtracking line search */
		
		__m256d Qnext0, Qnext1, Qnext2;
		__m256d nt   = _mm256_set1_pd(-100);

		Qnext0 = _mm256_fmadd_pd(nt, g0, Q0);	
		Qnext1 = _mm256_fmadd_pd(nt, g1, Q1);	
		Qnext2 = _mm256_fmadd_pd(nt, g2, Q2);

		while (_cost(pPtr, n, Qnext0, Qnext1, Qnext2) > 
																					costQ_0 + ALPHA*nt[0]*normFro_sq) {
			nt = _mm256_mul_pd(nt, vBETA);

			Qnext0 = _mm256_fmadd_pd(nt, g0, Q0);	
			Qnext1 = _mm256_fmadd_pd(nt, g1, Q1);	
			Qnext2 = _mm256_fmadd_pd(nt, g2, Q2);

		} /* back tracking while loop */

		Q0 = Qnext0;
		Q1 = Qnext1;
		Q2 = Qnext2;

	} /* main while loop */
	_mm256_maskstore_pd(&(*Q)[0][0], vMASK, Q0);
	_mm256_maskstore_pd(&(*Q)[1][0], vMASK, Q1);
	_mm256_maskstore_pd(&(*Q)[2][0], vMASK, Q2);
}

#endif




#if 0

/*****************************************************************************
 *  Preprocessor 
 ****************************************************************************/

#define VAR(var,i) var##_##i

/* assumes: pPtr, mask, costQ_i
 * declares: p_i, p0_i, p1_i, p2_i, t0_i, t1_i, t2_i, t3_i 
 */
#if 1
#define RESIDUAL_ACC(k,i) \
	__m256d VAR(p,k), VAR(p0,k), VAR(p1,k), VAR(p2,k);\
	VAR(p ,k) = _mm256_maskload_pd(&pPtr[i+k][0], vMASK);\
	VAR(p0,k) = _mm256_broadcast_sd(&pPtr[i+k][0]);\
	VAR(p1,k) = _mm256_broadcast_sd(&pPtr[i+k][1]);\
	VAR(p2,k) = _mm256_broadcast_sd(&pPtr[i+k][2]);\
	\
	__m256d VAR(ppT0,k), VAR(ppT1,k), VAR(ppT2,k);\
	\
	VAR(ppT0,k) = _mm256_mul_pd(VAR(p0,k), VAR(p,k));\
	VAR(ppT1,k) = _mm256_mul_pd(VAR(p1,k), VAR(p,k));\
	VAR(ppT2,k) = _mm256_mul_pd(VAR(p2,k), VAR(p,k));\
	\
	__m256d VAR(QppT0,k);\
	VAR(QppT0,k) = _mm256_mul_pd(VAR(ppT0,k), Q0);\
	\
	__m256d VAR(t0,k), VAR(t1,k);\
	VAR(t0,k) = _mm256_fmadd_pd(VAR(ppT1,k), Q1, VAR(QppT0,k));\
	VAR(t1,k) = _mm256_fmadd_pd(VAR(ppT2,k), Q2, VAR(t0,k));\
	\
	double VAR(residual,k) = VAR(t1,k)[0] + VAR(t1,k)[1] + VAR(t1,k)[2] - 1;\
	VAR(costQ,k) += VAR(residual,k)*VAR(residual,k);

#else
	/* this version is slower */
#define RESIDUAL_ACC(k,i) \
	__m256d VAR(p,k), VAR(p0,k), VAR(p1,k), VAR(p2,k);\
	VAR(p ,k) = _mm256_load_pd(&pPtr[i+k][0]);\
	VAR(p, k) = _mm256_and_pd(VAR(p,k), _mm256_castsi256_pd(vMASK));\
	VAR(p0,k) = _mm256_permute4x64_pd(VAR(p,k), 0b00000000);\
	VAR(p1,k) = _mm256_permute4x64_pd(VAR(p,k), 0b01010101);\
	VAR(p2,k) = _mm256_permute4x64_pd(VAR(p,k), 0b10101010);\
	\
	__m256d VAR(ppT0,k), VAR(ppT1,k), VAR(ppT2,k);\
	\
	VAR(ppT0,k) = _mm256_mul_pd(VAR(p0,k), VAR(p,k));\
	VAR(ppT1,k) = _mm256_mul_pd(VAR(p1,k), VAR(p,k));\
	VAR(ppT2,k) = _mm256_mul_pd(VAR(p2,k), VAR(p,k));\
	\
	__m256d VAR(QppT0,k);\
	VAR(QppT0,k) = _mm256_mul_pd(VAR(ppT0,k), Q0);\
	\
	__m256d VAR(t0,k), VAR(t1,k);\
	VAR(t0,k) = _mm256_fmadd_pd(VAR(ppT1,k), Q1, VAR(QppT0,k));\
	VAR(t1,k) = _mm256_fmadd_pd(VAR(ppT2,k), Q2, VAR(t0,k));\
	\
	__m256d VAR(t8,k), VAR(t9,k), VAR(t10,k);\
	VAR(t8,k) = _mm256_permute4x64_pd(VAR(t1,k), 0b00000000);\
	VAR(t9,k) = _mm256_permute4x64_pd(VAR(t1,k), 0b01010101);\
	VAR(t10,k) = _mm256_permute4x64_pd(VAR(t1,k), 0b10101010);\
	VAR(t8,k) = _mm256_add_pd(VAR(t8,k), VAR(t9,k));\
	VAR(t8,k) = _mm256_add_pd(VAR(t8,k), VAR(t10,k));\
	\
	__m256d VAR(residual,k);\
	VAR(residual,k) = _mm256_sub_pd(VAR(t8,k), vONES);\
	VAR(costQ,k) += VAR(residual,k)[0]*VAR(residual,k)[0];
#endif

/*****************************************************************************
 *  Private data
 ****************************************************************************/

static __m256d vONES, vBETA, vTWOS;
static __m256i vMASK;

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
static double _cost(const double (*pPtr)[3], int n,
		__m256d Q0, __m256d Q1, __m256d Q2)
{
	double costQ_0 = 0;
	double costQ_1 = 0;

	int i;
	for (i=0; i<n-2; i+=2) { /* unrolling twice: helps a bit */
		RESIDUAL_ACC(0, i);
		RESIDUAL_ACC(1, i);
	}

	/* cleanup */
	for (; i<n; i++) {
		RESIDUAL_ACC(0, i);
	}

	return costQ_0 + costQ_1;
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
void fit_ellipsoid_simd(const double (*pPtr)[3], int n, double (*Q)[3][3])
{
	__m256d Q0, Q1, Q2;
	Q0 = _mm256_setr_pd(1, 0, 0, 0);
	Q1 = _mm256_setr_pd(0, 1, 0, 0);
	Q2 = _mm256_setr_pd(0, 0, 1, 0);
	
	vMASK = _mm256_setr_epi64x(0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0);
	vBETA = _mm256_set1_pd(BETA);
	vONES = _mm256_set1_pd(1);
	vTWOS = _mm256_set1_pd(2);

	while(1) {
		/* determine gradient */

		__m256d g0, g1, g2;
		g0 = _mm256_setzero_pd();
		g1 = _mm256_setzero_pd();
		g2 = _mm256_setzero_pd();

		/* can comptue cost at same time. note: flops probably slightly reduce */
		double costQ_0 = 0;

		int i;
		for (i=0; i<n; i++) {
			RESIDUAL_ACC(0, i);

			__m256d t0;
			// t0 = _mm256_mul_pd(residual_0, vTWOS); 
			t0 = _mm256_set1_pd(residual_0*2);

			g0 = _mm256_fmadd_pd(t0, ppT0_0, g0);
			g1 = _mm256_fmadd_pd(t0, ppT1_0, g1);
			g2 = _mm256_fmadd_pd(t0, ppT2_0, g2);
		}

		__m256d g02, g12, g22;
		g02 = _mm256_mul_pd(g0, g0);
		g12 = _mm256_mul_pd(g1, g1);
		g22 = _mm256_mul_pd(g2, g2);

		g02 = _mm256_add_pd(g02, g12);
		g02 = _mm256_add_pd(g02, g22);

		double normFro_sq = g02[0] + g02[1] + g02[2];
		if (normFro_sq < EPSILON*EPSILON) break;

		/* take step direction to be negative gradient */
		
		/* backtracking line search */
		
		__m256d Qnext0, Qnext1, Qnext2;
		__m256d nt   = _mm256_set1_pd(-100);

		Qnext0 = _mm256_fmadd_pd(nt, g0, Q0);
		Qnext1 = _mm256_fmadd_pd(nt, g1, Q1);
		Qnext2 = _mm256_fmadd_pd(nt, g2, Q2);

		while (_cost(pPtr, n, Qnext0, Qnext1, Qnext2) > 
																					costQ_0 + ALPHA*nt[0]*normFro_sq) {
			nt = _mm256_mul_pd(nt, vBETA);

			Qnext0 = _mm256_fmadd_pd(nt, g0, Q0);	
			Qnext1 = _mm256_fmadd_pd(nt, g1, Q1);	
			Qnext2 = _mm256_fmadd_pd(nt, g2, Q2);

		} /* back tracking while loop */

		Q0 = Qnext0;
		Q1 = Qnext1;
		Q2 = Qnext2;

	} /* main while loop */
	_mm256_maskstore_pd(&(*Q)[0][0], vMASK, Q0);
	_mm256_maskstore_pd(&(*Q)[1][0], vMASK, Q1);
	_mm256_maskstore_pd(&(*Q)[2][0], vMASK, Q2);
}
#else
/*****************************************************************************
 *  Unrolled
 ****************************************************************************/

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
void fit_ellipsoid_simd(const double (*pPtr)[3], int n, double (*Q)[3][3])
{
	__m256d Q0, Q1, Q2;
	Q0 = _mm256_setr_pd(1, 0, 0, 0);
	Q1 = _mm256_setr_pd(0, 1, 0, 0);
	Q2 = _mm256_setr_pd(0, 0, 1, 0);
	
	vMASK = _mm256_setr_epi64x(0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0);
	vBETA = _mm256_set1_pd(BETA);
	vONES = _mm256_set1_pd(1);
	vTWOS = _mm256_set1_pd(2);

	while(1) {
		/* determine gradient */

		__m256d g0, g1, g2, g0_1, g1_1, g2_1;
		g0 = _mm256_setzero_pd();
		g1 = _mm256_setzero_pd();
		g2 = _mm256_setzero_pd();

		g0_1 = _mm256_setzero_pd();
		g1_1 = _mm256_setzero_pd();
		g2_1 = _mm256_setzero_pd();

		/* can comptue cost at same time. note: flops probably slightly reduce */
		double costQ_0 = 0, costQ_1 = 0;

		int i;
		for (i=0; i<n; i+=2) {
			RESIDUAL_ACC(0, i);
			RESIDUAL_ACC(1, i);

			__m256d t0, t1;
			t0 = _mm256_mul_pd(residual_0, vTWOS); 
			t1 = _mm256_mul_pd(residual_1, vTWOS); 
			// t0 = _mm256_set1_pd(residual_0*2);
			// t1 = _mm256_set1_pd(residual_1*2);

			__m256d t3, t4, t5, t6, t7, t8;
			t3 = _mm256_mul_pd(p0_0, p_0); 
			t4 = _mm256_mul_pd(p1_0, p_0);
			t5 = _mm256_mul_pd(p2_0, p_0);
			
			t6 = _mm256_mul_pd(p0_1, p_1); 
			t7 = _mm256_mul_pd(p1_1, p_1);
			t8 = _mm256_mul_pd(p2_1, p_1);
			
			g0   = _mm256_fmadd_pd(t0, t3, g0);
			g1   = _mm256_fmadd_pd(t0, t4, g1);
			g2   = _mm256_fmadd_pd(t0, t5, g2);
			g0_1 = _mm256_fmadd_pd(t1, t6, g0_1);
			g1_1 = _mm256_fmadd_pd(t1, t7, g1_1);
			g2_1 = _mm256_fmadd_pd(t1, t8, g2_1);
		}
		/* cleanup */
		costQ_0 += costQ_1;
		g0 = _mm256_add_pd(g0, g0_1);
		g1 = _mm256_add_pd(g1, g1_1);
		g2 = _mm256_add_pd(g2, g2_1);

		/* compute frobenius norm sq */
		__m256d g02, g12, g22;
		g02 = _mm256_mul_pd(g0, g0);
		g12 = _mm256_mul_pd(g1, g1);
		g22 = _mm256_mul_pd(g2, g2);

		g02 = _mm256_add_pd(g02, g12);
		g02 = _mm256_add_pd(g02, g22);

		double normFro_sq = g02[0] + g02[1] + g02[2];
		if (normFro_sq < EPSILON*EPSILON) break;

		/* take step direction to be negative gradient */
		
		/* backtracking line search */
		
		__m256d Qnext0, Qnext1, Qnext2;
		__m256d nt   = _mm256_set1_pd(-100);

		Qnext0 = _mm256_fmadd_pd(nt, g0, Q0);	
		Qnext1 = _mm256_fmadd_pd(nt, g1, Q1);	
		Qnext2 = _mm256_fmadd_pd(nt, g2, Q2);

		while (_cost(pPtr, n, Qnext0, Qnext1, Qnext2) > 
																					costQ_0 + ALPHA*nt[0]*normFro_sq) {
			nt = _mm256_mul_pd(nt, vBETA);

			Qnext0 = _mm256_fmadd_pd(nt, g0, Q0);	
			Qnext1 = _mm256_fmadd_pd(nt, g1, Q1);	
			Qnext2 = _mm256_fmadd_pd(nt, g2, Q2);

		} /* back tracking while loop */

		Q0 = Qnext0;
		Q1 = Qnext1;
		Q2 = Qnext2;

	} /* main while loop */
	_mm256_maskstore_pd(&(*Q)[0][0], vMASK, Q0);
	_mm256_maskstore_pd(&(*Q)[1][0], vMASK, Q1);
	_mm256_maskstore_pd(&(*Q)[2][0], vMASK, Q2);
}

#endif

#endif