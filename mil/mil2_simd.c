/**    ____                   __  ___
*     / __ )____  ____  ___  /  |/  /___ _____
*    / __  / __ \/ __ \/ _ \/ /|_/ / __ `/ __ \
*   / /_/ / /_/ / / / /  __/ /  / / /_/ / /_/ /
*  /_____/\____/_/ /_/\___/_/  /_/\__,_/ .___/
*                                     /_/
*  Bone Anisotropy Mapping.
*
*  Copyright (C) 2019  Jarunan Panyasantisuk (jarunan@ethz.ch)
*                      Rajan Gill            (rgill@ethz.ch)
*                      Ryan Cherifa          (rcherifa@student.ethz.ch)
*                      Joao Rivera           (hector.rivera@inf.ethz.ch)
*
*  This program is free software: you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation, either version 3 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/
#include "mil2_simd.h"
#include "ellipsoid.h"
#include <immintrin.h>
//#include <stdio.h>

//#define DEBUG
int gBone1, gBone2, gInter1, gInter2;
unsigned int count_simd = 0;
unsigned int already_tested_simd[13] = {0};

#if 0
void simd_mil_test_all_old(const double *hr_sphere_region, int n, double *directions_vectors_mil) {

    double bone_length[NUM_DIRECTIONS];
    int intercepts[NUM_DIRECTIONS];

    for (int i = 0; i < NUM_DIRECTIONS; ++i) {
        bone_length[i] = 0.0;
        intercepts[i]  = 0;
    }

    for (int kk = 0; kk < n; kk+=BLOCK_SIZE) {
        for (int jj = 0; jj < n; jj+=BLOCK_SIZE) {
            for (int ii = 0; ii < n; ii+=BLOCK_SIZE) {
                int intercept_blk;

                bone_length[0] += simd_mil_1D(hr_sphere_region, &intercept_blk, n, kk, jj, ii, 1);
                intercepts[0]  += intercept_blk;
                bone_length[1] += simd_mil_1D(hr_sphere_region, &intercept_blk, n, kk, ii, jj, 2);
                intercepts[1]  += intercept_blk;
                bone_length[2] += simd_mil_1D(hr_sphere_region, &intercept_blk, n, jj, ii, kk, 3);
                intercepts[2]  += intercept_blk;


                bone_length[3] += simd_mil_2D_pos(hr_sphere_region, &intercept_blk, n, kk, jj, ii, 4);
                intercepts[3]  += intercept_blk;
                bone_length[4] += simd_mil_2D_pos(hr_sphere_region, &intercept_blk, n, jj, kk, ii, 5);
                intercepts[4]  += intercept_blk;
                bone_length[5] += simd_mil_2D_pos(hr_sphere_region, &intercept_blk, n, ii, jj, kk, 6);
                intercepts[5]  += intercept_blk;

                bone_length[6] += simd_mil_2D_neg(hr_sphere_region, &intercept_blk, n, kk, jj, ii, 7);
                intercepts[6]  += intercept_blk;
                bone_length[7] += simd_mil_2D_neg(hr_sphere_region, &intercept_blk, n, jj, kk, ii, 8);
                intercepts[7]  += intercept_blk;
                bone_length[8] += simd_mil_2D_neg(hr_sphere_region, &intercept_blk, n, ii, jj, kk, 9);
                intercepts[8]  += intercept_blk;

            }
        }
    }

    for (int i = 0; i < NUM_DIRECTIONS; ++i) {
        if (intercepts[i] == 0) {
            intercepts[i] = 1;
        }
        directions_vectors_mil[i] = bone_length[i] / intercepts[i];
    }

#ifdef DEBUG
    printf("SIMD 0 - BONE LENGTH = %.8f, INTERCEPTS = %d\n", bone_length[0], intercepts[0]);
    printf("SIMD 1 - BONE LENGTH = %.8f, INTERCEPTS = %d\n", bone_length[1], intercepts[1]);
    printf("SIMD 2 - BONE LENGTH = %.8f, INTERCEPTS = %d\n", bone_length[2], intercepts[2]);
    printf("SIMD 3 - BONE LENGTH = %.8f, INTERCEPTS = %d\n", bone_length[3], intercepts[3]);
    printf("SIMD 4 - BONE LENGTH = %.8f, INTERCEPTS = %d\n", bone_length[4], intercepts[4]);
    printf("SIMD 5 - BONE LENGTH = %.8f, INTERCEPTS = %d\n", bone_length[5], intercepts[5]);
    printf("SIMD 6 - BONE LENGTH = %.8f, INTERCEPTS = %d\n", bone_length[6], intercepts[6]);
    printf("SIMD 7 - BONE LENGTH = %.8f, INTERCEPTS = %d\n", bone_length[7], intercepts[7]);
    printf("SIMD 8 - BONE LENGTH = %.8f, INTERCEPTS = %d\n", bone_length[8], intercepts[8]);
#endif
    gBone2 = bone_length[0];
    gInter2 = intercepts[0];

}
#endif

inline double horizontal_add (__m256d a) {
    __m256d t1 = _mm256_hadd_pd(a,a);
    __m256d t2 = _mm256_permute4x64_pd(t1, 0b00011011);
    __m256d t3 = _mm256_add_pd(t2,t1);
    __m128d t4 = _mm256_extractf128_pd(t3,0);
    return _mm_cvtsd_f64(t4);
}

inline int horizontal_addi(__m256i a) {
    __int64_t x = _mm256_extract_epi64(a, 0);
    __int64_t x1 = _mm256_extract_epi64(a, 1);
    __int64_t x2 = _mm256_extract_epi64(a, 2);
    __int64_t x3 = _mm256_extract_epi64(a, 3);
    __int64_t x4 = x + x1;
    __int64_t x5 = x2 + x3;
    __int64_t x6 = x5 + x4;
    return x6;
}
#if 0
inline double simd_mil_1D(const double *hr_sphere_region, int* intercepts, int n, const int kk, const int jj, const int ii,  const int vecID) {
    double bone_length;

    /* Init accumulators */
    __m256d bone_count = _mm256_setzero_pd();
    __m256d bone_count2 = _mm256_setzero_pd();

    __m256i edge_count = _mm256_set1_epi64x(0);
    __m256i edge_count2 = _mm256_set1_epi64x(0);

    int nn = n * n;
    int i_prev = (ii > 0) ? ii - 1 : 0;

    for (int k = kk + 1; k < kk + BLOCK_SIZE; k += STRIDE) {
        int knn = k * nn;
        for (int j = jj + 1; j < jj + BLOCK_SIZE; j += STRIDE*NUM_ACC*2) {
            __m256i prev_mask, prev_mask2;

            SIMD_LOAD_PREV_1D

            for (int i = ii; i < ii + BLOCK_SIZE; ++i) {
                int inn = i * nn;
                __m256d region1;
                __m256d region2;

                /* Load working set */
                SIMD_LOAD_DATA_SET_1D

                /* Perform computation */
                SIMD_COMPUTATION

                /* Update state of prev_mask */
                prev_mask = curr_mask;
                prev_mask2 = curr_mask2;
#ifdef  DEBUG
                count_simd++;
#endif
            }
        }
    }
    bone_count = _mm256_add_pd(bone_count, bone_count2);
    edge_count = _mm256_add_epi64(edge_count, edge_count2);

    bone_length  = horizontal_add(bone_count);
    *intercepts  = horizontal_addi(edge_count);

#ifdef  DEBUG
    if (!already_tested_simd[vecID]) {
        printf("%d\n", count_simd);
        already_tested_simd[vecID] = 1;
    }
    count_simd = 0;
#endif

    return bone_length;
}

double simd_mil_2D_pos(const double *hr_sphere_region, int* intercepts, int n, const int kk, const int jj, const int ii,  const int vecID) {
    double bone_length;

    /* Init accumulators */
    __m256d bone_count = _mm256_setzero_pd();
    __m256d bone_count2 = _mm256_setzero_pd();

    __m256i edge_count = _mm256_set1_epi64x(0);
    __m256i edge_count2 = _mm256_set1_epi64x(0);


    for (int k = kk + 1; k < kk + BLOCK_SIZE; k += STRIDE * NUM_ACC) {
        for (int ij = 0; ij < BLOCK_SIZE; ij += STRIDE) {
            unsigned int i1_prev, i2_prev, j1_prev, j2_prev;
            __m256i prev_mask, prev_mask2;

            int i1 = ii + ij;
            int j1 = jj;
            int i2 = ii;
            int j2 = jj + ij;

            /* Initialise previous mask */
            SIMD_LOAD_PREV_2D_POS

            while (i1 + 1 < ii + BLOCK_SIZE && j2 + 1 < jj + BLOCK_SIZE) {
                __m256d region1;
                __m256d region2;

                /* Load working set */
                LOAD_DATA_SET_2D_POS

                /* Perform computation */
                SIMD_COMPUTATION

                /* Update state of prev_mask */
                prev_mask = curr_mask;
                prev_mask2 = curr_mask2;
#ifdef  DEBUG
                count_simd++;
#endif

                ++i1;
                ++j1;
                ++i2;
                ++j2;
            }
        }
    } /* End iteration over dimension k */

    bone_count = _mm256_add_pd(bone_count, bone_count2);
    edge_count = _mm256_add_epi64(edge_count, edge_count2);

    bone_length  = horizontal_add(bone_count);
    *intercepts  = horizontal_addi(edge_count);

#ifdef  DEBUG
    if (!already_tested_simd[vecID]) {
        printf("%d\n", count_simd);
        already_tested_simd[vecID] = 1;
    }
    count_simd = 0;
#endif

    return bone_length;

}

double simd_mil_2D_neg(const double *hr_sphere_region, int* intercepts, int n, const int kk, const int jj, const int ii,  const int vecID) {
    double bone_length;

    /* Init accumulators */
    __m256d bone_count = _mm256_setzero_pd();
    __m256d bone_count2 = _mm256_setzero_pd();

    __m256i edge_count = _mm256_set1_epi64x(0);
    __m256i edge_count2 = _mm256_set1_epi64x(0);


    for (int k = kk + 1; k < kk + BLOCK_SIZE; k += STRIDE * NUM_ACC) {
        for (int ij = 0; ij < BLOCK_SIZE; ij += STRIDE) {
            unsigned int i1_prev, i2_prev, j1_prev, j2_prev;
            __m256i prev_mask, prev_mask2;

            int i1 = ii + (BLOCK_SIZE-1) - ij;  /* Start at the end of the row in the block */
            int j1 = jj;
            int i2 = ii + (BLOCK_SIZE-1);
            int j2 = jj + ij;

            /* Initialise previous mask */
            LOAD_PREV_2D_NEG

            while (j2 + 1 < jj + BLOCK_SIZE) {
                __m256d region1;
                __m256d region2;

                /* Load working set */
                LOAD_DATA_SET_2D_NEG

                /* Perform computation */
                SIMD_COMPUTATION

                /* Update state of prev_mask */
                prev_mask = curr_mask;
                prev_mask2 = curr_mask2;

                --i1;
                ++j1;
                --i2;
                ++j2;

#ifdef  DEBUG
                count_simd++;
#endif
            }
        }
    } /* End iteration over dimension k */

    bone_count = _mm256_add_pd(bone_count, bone_count2);
    edge_count = _mm256_add_epi64(edge_count, edge_count2);

    bone_length  = horizontal_add(bone_count);
    *intercepts  = horizontal_addi(edge_count);

#ifdef  DEBUG
    if (!already_tested_simd[vecID]) {
        printf("%d\n", count_simd);
        already_tested_simd[vecID] = 1;
    }
    count_simd = 0;
#endif

    return bone_length;

}

#endif
///
/// Test all vectors with SIMD.
///
void mil2_simd(const double *hr_sphere_region, int n, double *directions_vectors_mil) {

    double bone_length[13] = {0.0};
    int intercepts[13] = {0};
    for (int kk_b = 0; kk_b < n; kk_b+=BLOCK_SIZE) {
        for (int jj_b = 0; jj_b < n; jj_b+=BLOCK_SIZE) {
            for (int ii_b = 0; ii_b < n; ii_b+=BLOCK_SIZE) {

//                for (int v = 0; v < 2; ++v) {
                    BLOCK_KERNEL_1D_SIMD(1, kk_b, jj_b, ii_b)
                    BLOCK_KERNEL_1D_SIMD(2, kk_b, ii_b, jj_b)
                    BLOCK_KERNEL_1D_SIMD(3, jj_b, ii_b, kk_b)
                    BLOCK_KERNEL_2D_POS_SIMD(4, kk_b, jj_b, ii_b)
                    BLOCK_KERNEL_2D_POS_SIMD(5, jj_b, kk_b, ii_b)
                    BLOCK_KERNEL_2D_POS_SIMD(6, ii_b, jj_b, kk_b)
//                }
            }
        }
    }

    for (int i = 0; i < 13; ++i) {
        if (intercepts[i] == 0) {
            intercepts[i] = 1;
        }
        directions_vectors_mil[i] = bone_length[i] / intercepts[i];
    }
}