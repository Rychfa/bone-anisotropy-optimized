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
#include "mil2.h"
#include "ellipsoid.h"
#include <immintrin.h>
#include <stdio.h>

//#define DEBUG
int gBone1, gBone2, gInter1, gInter2;
unsigned int count = 0;
unsigned int already_tested[13] = {0};

static int facesVectors[3][3] =
        {
             /* { i,  j,  k} */
                { 0,  1,  1},
                { 1,  0,  1},
                { 1,  1,  0},
        };

int get_iterator_vectors (const int directionVector[3], int iteratorVectors[][3]) {

    int n = 0;
    for (int i = 0; i < 3; ++i) {
        if (directionVector[i] != 0) {
            iteratorVectors[n][0] = facesVectors[i][0];
            iteratorVectors[n][1] = facesVectors[i][1];
            iteratorVectors[n][2] = facesVectors[i][2];
            ++n;
        }
    }

    return n;
}

void mil2(const double *hr_sphere_region, int n, double *directions_vectors_mil) {

    const int n_vectors = NUM_DIRECTIONS;

    double directions_vectors_bone_length[n_vectors], directions_vectors_intercepts[n_vectors];
    int iteratorVectors[3][3];

    for (int j = 0; j < n_vectors; ++j) {
        directions_vectors_bone_length[j] = 0.0;
        directions_vectors_intercepts[j]  = 1.0;
    }


    /* for every direction vector */
    for (int v = 0; v < n_vectors; ++v) {

        int numIterVecs = 0;
        /* get array of iterator vectors */
        numIterVecs = get_iterator_vectors(DIRECTIONS[v], iteratorVectors);

        /* Iterate over all initial faces */
        for (int f = 0; f < numIterVecs;  ++f) {
            /* Iterate through initial points in the face */
            for (int kk = 0 ; kk <= (n-1)*iteratorVectors[f][2]; kk += STRIDE) {
                for (int jj = 0 ; jj <= (n-1)*iteratorVectors[f][1]; jj += STRIDE) {
                    for (int ii = 0 ; ii <= (n-1)*iteratorVectors[f][0]; ii += STRIDE) {
                        int k = (DIRECTIONS[v][2] == -1 && iteratorVectors[f][2] == 0)? n - 1 : kk;
                        int j = (DIRECTIONS[v][1] == -1 && iteratorVectors[f][1] == 0)? n - 1 : jj;
                        int i = (DIRECTIONS[v][0] == -1 && iteratorVectors[f][0] == 0)? n - 1 : ii;
                        unsigned int current_point;
                        unsigned int prev_point = hr_sphere_region[ k*n*n + j*n + i];
                        while ((k < n && k >= 0) && (j < n && j >= 0) && (i < n && i >= 0)) {

                            current_point = hr_sphere_region[ k*n*n + j*n + i];

                            directions_vectors_bone_length[v] += current_point;
                            directions_vectors_intercepts[v]  += (double) (current_point ^ prev_point);
                            prev_point = current_point;

                            k += DIRECTIONS[v][2];
                            j += DIRECTIONS[v][1];
                            i += DIRECTIONS[v][0];
                        }
                    }
                }
            }

        }

    }

    for (int i = 0; i < n_vectors; ++i) {
        directions_vectors_mil[i] = directions_vectors_bone_length[i] / directions_vectors_intercepts[i];
    }

}

int get_start_position (int ii, int n, int dir, int iter) {
    int i = ii + iter;
    if (dir == -1) {
        if (iter == 0) {
            i = n - 1;
        } else {
            i = ii;
        }
    }
    return i;
}

void mil2_baseline(const double *hr_sphere_region, int n, double *directions_vectors_mil) {

    const int n_vectors = NUM_DIRECTIONS;

    double directions_vectors_bone_length[n_vectors];
    int directions_vectors_intercepts[n_vectors];
    int iteratorVectors[3][3];

    for (int j = 0; j < n_vectors; ++j) {
        directions_vectors_bone_length[j] = 0.0;
        directions_vectors_intercepts[j]  = 0;
    }

    /* for every direction vector */
    for (int v = 0; v < n_vectors; ++v) {

        int numIterVecs = 0;
        /* get array of iterator vectors */
        numIterVecs = get_iterator_vectors(DIRECTIONS[v], iteratorVectors);

        /* Iterate over all initial faces */
        for (int f = 0; f < numIterVecs;  ++f) {
            /* Iterate through initial points in the face */
            for (int kk = 0 ; kk <= (n-1)*iteratorVectors[f][2]; kk += STRIDE) {
                for (int jj = 0 ; jj <= (n-1)*iteratorVectors[f][1]; jj += STRIDE) {
                    for (int ii = 0 ; ii <= (n-1)*iteratorVectors[f][0]; ii += STRIDE) {
                        int k = get_start_position(kk, n, DIRECTIONS[v][2], iteratorVectors[f][2]);
                        int j = get_start_position(jj, n, DIRECTIONS[v][1], iteratorVectors[f][1]);
                        int i = get_start_position(ii, n, DIRECTIONS[v][0], iteratorVectors[f][0]);

                        unsigned int current_mask;
                        unsigned int prev_mask = hr_sphere_region[ k*n*n + j*n + i] > 0.5;
                        while ((k < n && k >= 0) && (j < n && j >= 0) && (i < n && i >= 0)) {

                            directions_vectors_bone_length[v] += hr_sphere_region[ k*n*n + j*n + i];
                            current_mask = hr_sphere_region[ k*n*n + j*n + i] > 0.5;
                            directions_vectors_intercepts[v]  += current_mask ^ prev_mask;
                            prev_mask = current_mask;

                            i += DIRECTIONS[v][0];
                            j += DIRECTIONS[v][1];
                            k += DIRECTIONS[v][2];
                        }
                    }
                }
            }

        }

    }

    for (int i = 0; i < n_vectors; ++i) {
        if (directions_vectors_intercepts[i] == 0) {
            directions_vectors_intercepts[i] = 1;
        }

        directions_vectors_mil[i] = directions_vectors_bone_length[i] / directions_vectors_intercepts[i];
    }

#ifdef DEBUG
    printf("TRUE 0 - BONE LENGTH = %.8f, INTERCEPTS = %d\n", directions_vectors_bone_length[0], directions_vectors_intercepts[0]);
    printf("TRUE 1 - BONE LENGTH = %.8f, INTERCEPTS = %d\n", directions_vectors_bone_length[1], directions_vectors_intercepts[1]);
    printf("TRUE 2 - BONE LENGTH = %.8f, INTERCEPTS = %d\n", directions_vectors_bone_length[2], directions_vectors_intercepts[2]);
#endif
    gBone1 = directions_vectors_bone_length[3];
    gInter1 = directions_vectors_intercepts[3];

}

#if false
///
/// Blocking version for vectors (1,0,0), (0,1,0), (0,0,1).
/// Using accumulators and unrolling.
///
inline double mil_1D(const double *hr_sphere_region, int* intercepts, int n, const int kk, const int jj, const int ii,  const int vecID) {
    double bone_length;

    /* Init accumulators */
    double acc1 = 0.0;
    double acc2 = 0.0;
    double acc3 = 0.0;
    double acc4 = 0.0;

    unsigned int edge_count1 = 0;
    unsigned int edge_count2 = 0;
    unsigned int edge_count3 = 0;
    unsigned int edge_count4 = 0;

    for (int k = kk + 1; k < kk + BLOCK_SIZE; k += STRIDE) {
        for (int j = jj + 1; j < jj + BLOCK_SIZE; j += STRIDE*NUM_ACC) {
            unsigned int i_prev, prev_mask1, prev_mask2, prev_mask3, prev_mask4;
            LOAD_PREV_1D

            for (int i = ii; i < ii + BLOCK_SIZE; ++i) {
                double r1, r2, r3, r4;

                /* Load working set */
                LOAD_DATA_SET_1D

                /* Perform computation */
                COMPUTATION

                /* Update state of prev_mask */
                prev_mask1 = curr_mask1;
                prev_mask2 = curr_mask2;
                prev_mask3 = curr_mask3;
                prev_mask4 = curr_mask4;
            }
        }
    }

    acc1 += acc2;
    acc3 += acc4;
    bone_length  = acc1 + acc3;

    edge_count1 += edge_count2;
    edge_count3 += edge_count4;
    *intercepts  = edge_count1 + edge_count3;

#ifdef  DEBUG
    if (!already_tested[vecID]) {
        printf("%d\n", count);
        already_tested[vecID] = 1;
    }
    count = 0;
#endif


    return bone_length;

}

///
/// Blocking version for vectors (1,1,0), (0,1,1), (1,0,1).
/// Using unrolling and uneven accumulators.
///
double mil_2D_pos_uneven(const double *hr_sphere_region, int* intercepts, int n, const int kk, const int jj, const int ii,  const int vecID) {
    double bone_length;

    /* Init accumulators */
    double acc1 = 0.0;
    double acc2 = 0.0;
    double acc3 = 0.0;
    double acc4 = 0.0;

    unsigned int edge_count1 = 0;
    unsigned int edge_count2 = 0;
    unsigned int edge_count3 = 0;
    unsigned int edge_count4 = 0;

    for (int k = kk + 1; k < kk + BLOCK_SIZE; k += STRIDE) {
        for (int ij = 0; ij < BLOCK_SIZE; ij += STRIDE*NUM_ACC/2) {
            unsigned int i1_prev, i2_prev, j1_prev, j2_prev;
            unsigned int prev_mask1, prev_mask2, prev_mask3, prev_mask4;

            int i1 = ii + ij;
            int j1 = jj;
            int i2 = ii;
            int j2 = jj + ij;

            /* Initialise previous mask */
            SIMD_LOAD_PREV_2D_POS

            while (i1 + 1 < ii + BLOCK_SIZE && j2 + 1 < jj + BLOCK_SIZE) {
                double r1, r2, r3, r4;

                /* Load working set */
                LOAD_DATA_SET_2D_POS

                /* Perform computation */
                COMPUTATION

                /* Update state of prev_mask */
                prev_mask1 = curr_mask1;
                prev_mask2 = curr_mask2;
                prev_mask3 = curr_mask3;
                prev_mask4 = curr_mask4;
                count++;


                ++i1;
                ++j1;
                ++i2;
                ++j2;
            }
        }
    } /* End iteration over dimension k */

    /* Sum up accumulators */
    acc1 += acc2;
    acc3 += acc4;
    bone_length = acc1 + acc3;

    edge_count1 += edge_count2;
    edge_count3 += edge_count4;
    *intercepts  = edge_count1 + edge_count3;

    return bone_length;
}

///
/// Blocking version for vectors (1,1,0), (0,1,1), (1,0,1).
/// Using accumulators and unrolling over all dimensions.
///
double mil_2D_pos(const double *hr_sphere_region, int* intercepts, int n, const int kk, const int jj, const int ii,  const int vecID) {
    double bone_length;

    /* Init accumulators */
    double acc1 = 0.0;
    double acc2 = 0.0;
    double acc3 = 0.0;
    double acc4 = 0.0;

    unsigned int edge_count1 = 0;
    unsigned int edge_count2 = 0;
    unsigned int edge_count3 = 0;
    unsigned int edge_count4 = 0;

    for (int k = kk + 1; k < kk + BLOCK_SIZE; k += STRIDE*NUM_ACC/2) {
        for (int ij = 0; ij < BLOCK_SIZE; ij += STRIDE) {
            unsigned int i1_prev, i2_prev, j1_prev, j2_prev;
            unsigned int prev_mask1, prev_mask2, prev_mask3, prev_mask4;

            int i1 = ii + ij;
            int j1 = jj;
            int i2 = ii;
            int j2 = jj + ij;

            /* Initialise previous mask */
            SIMD_LOAD_PREV_2D_POS

            while (i1 + 1 < ii + BLOCK_SIZE && j2 + 1 < jj + BLOCK_SIZE) {
                double r1, r2, r3, r4;

                /* Load working set */
                LOAD_DATA_SET_2D_POS

                /* Perform computation */
                COMPUTATION

                /* Update state of prev_mask */
                prev_mask1 = curr_mask1;
                prev_mask2 = curr_mask2;
                prev_mask3 = curr_mask3;
                prev_mask4 = curr_mask4;
                count++;


                ++i1;
                ++j1;
                ++i2;
                ++j2;
            }
        }
    } /* End iteration over dimension k */

    /* Sum up accumulators */
    acc1 += acc2;
    acc3 += acc4;
    bone_length = acc1 + acc3;

    edge_count1 += edge_count2;
    edge_count3 += edge_count4;
    *intercepts  = edge_count1 + edge_count3;

    return bone_length;
}

///
/// Blocking version for vectors (-1,1,0), (0,-1,1), (-1,0,1).
/// Using accumulators and unrolling.
///
double mil_2D_neg(const double *hr_sphere_region, int* intercepts, int n, const int kk, const int jj, const int ii,  const int vecID) {
    double bone_length;

    /* Init accumulators */
    double acc1 = 0.0, acc5 = 0.0;
    double acc2 = 0.0, acc6 = 0.0;
    double acc3 = 0.0, acc7 = 0.0;
    double acc4 = 0.0, acc8 = 0.0;

    unsigned int edge_count1 = 0, edge_count5 = 0;
    unsigned int edge_count2 = 0, edge_count6 = 0;
    unsigned int edge_count3 = 0, edge_count7 = 0;
    unsigned int edge_count4 = 0, edge_count8 = 0;

    for (int k = kk + 1; k < kk + BLOCK_SIZE; k += STRIDE) {
        for (int ij = 0; ij < BLOCK_SIZE; ij += STRIDE*NUM_ACC/2) {
            unsigned int i1_prev, i2_prev, j1_prev, j2_prev;
            unsigned int prev_mask1, prev_mask2, prev_mask3, prev_mask4;

            int i1 = ii + (BLOCK_SIZE-1) - ij;  /* Start at the end of the row in the block */
            int j1 = jj;
            int i2 = ii + (BLOCK_SIZE-1);
            int j2 = jj + ij;

            /* Initialise previous mask */
            LOAD_PREV_2D_NEG

            while (j2 + 1 < jj + BLOCK_SIZE) {
                double r1, r2, r3, r4;

                /* Load working set */
                LOAD_DATA_SET_2D_NEG

                /* Perform computation */
                COMPUTATION

                /* Update state of prev_mask */
                prev_mask1 = curr_mask1;
                prev_mask2 = curr_mask2;
                prev_mask3 = curr_mask3;
                prev_mask4 = curr_mask4;
                count++;


                --i1;
                ++j1;
                --i2;
                ++j2;
            }
        }
    } /* End iteration over dimension k */

    acc1 += acc2;   edge_count1 += edge_count2;
    acc3 += acc4;   edge_count3 += edge_count4;
    acc5 += acc6;   edge_count5 += edge_count6;
    acc7 += acc8;   edge_count7 += edge_count8;
    acc1 += acc3;   edge_count1 += edge_count3;
    acc5 += acc7;   edge_count5 += edge_count7;
    bone_length = acc1 + acc5;
    *intercepts  = edge_count1 + edge_count5;

#ifdef  DEBUG
    if (!already_tested[vecID]) {
        printf("%d\n", count);
        already_tested[vecID] = 1;
    }
    count = 0;
#endif

    return bone_length;
}
#endif

///
/// Test all vectors.
///
void mil_test_all(const double *hr_sphere_region, int n, double *directions_vectors_mil) {

    double bone_length[13] = {0.0};
    int intercepts[13] = {0};
    for (int kk_b = 0; kk_b < n; kk_b+=BLOCK_SIZE) {
        for (int jj_b = 0; jj_b < n; jj_b+=BLOCK_SIZE) {
            for (int ii_b = 0; ii_b < n; ii_b+=BLOCK_SIZE) {

                for (int v = 0; v < 2; ++v) {
                    BLOCK_KERNEL_1D(1, kk_b, jj_b, ii_b)
                    BLOCK_KERNEL_1D(2, kk_b, ii_b, jj_b)
                    BLOCK_KERNEL_1D(3, jj_b, ii_b, kk_b)
                    BLOCK_KERNEL_2D(4, kk_b, jj_b, ii_b)
                    BLOCK_KERNEL_2D(5, jj_b, kk_b, ii_b)
                    BLOCK_KERNEL_2D(6, ii_b, jj_b, kk_b)
                }
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