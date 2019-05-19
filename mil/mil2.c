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
#include <stdio.h>

#define STRIDE 2
#define DEBUG
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

void mil2(const int *hr_sphere_region, int n, double *directions_vectors_mil) {

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

void mil2_baseline(const float *hr_sphere_region, int n, float *directions_vectors_mil) {

    const int n_vectors = NUM_DIRECTIONS;

    float directions_vectors_bone_length[n_vectors];
    int directions_vectors_intercepts[n_vectors];
    int iteratorVectors[3][3];

    for (int j = 0; j < n_vectors; ++j) {
        directions_vectors_bone_length[j] = 0.0;
        directions_vectors_intercepts[j]  = 1;
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
                        float current_point;
                        float prev_point = hr_sphere_region[ k*n*n + j*n + i];
                        while ((k < n && k >= 0) && (j < n && j >= 0) && (i < n && i >= 0)) {

                            current_point = hr_sphere_region[ k*n*n + j*n + i];

                            directions_vectors_bone_length[v] += current_point;
                            directions_vectors_intercepts[v]  += (current_point > 0.5) ^ (prev_point > 0.5);
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

///
/// Scalar replacement to remove aliasing
///
void mil2_o1(const double *hr_sphere_region, int n, double *directions_vectors_mil) {

    const int n_vectors = NUM_DIRECTIONS;

    int iteratorVectors[3][3];

    /* for every direction vector */
    for (int v = 0; v < n_vectors; ++v) {

        int numIterVecs = 0;
        /* get array of iterator vectors */
        numIterVecs = get_iterator_vectors(DIRECTIONS[v], iteratorVectors);
        int d0 = DIRECTIONS[v][0];
        int d1 = DIRECTIONS[v][1];
        int d2 = DIRECTIONS[v][2];
        double sum = 0.0;
        int intersects = 1;

        /* Iterate over all initial faces */
        for (int f = 0; f < numIterVecs;  ++f) {
            int iter0 = iteratorVectors[f][0];
            int iter1 = iteratorVectors[f][1];
            int iter2 = iteratorVectors[f][2];
            /* Iterate through initial points in the face */
            for (int kk = 0 ; kk <= (n-1)*iter2; kk += STRIDE) {
                for (int jj = 0 ; jj <= (n-1)*iter1; jj += STRIDE) {
                    for (int ii = 0 ; ii <= (n-1)*iter0; ii += STRIDE) {
                        int k = (d2 == -1 && iter2 == 0)? n - 1 : kk;
                        int j = (d1 == -1 && iter1 == 0)? n - 1 : jj;
                        int i = (d0 == -1 && iter0 == 0)? n - 1 : ii;
                        double current_point;
                        unsigned int current_masked;
                        unsigned int prev_masked = hr_sphere_region[ k*n*n + j*n + i] > 0.5;
                        while ((i < n && i >= 0) && (j < n && j >= 0) && (k < n && k >= 0) ) {

                            current_point = hr_sphere_region[ k*n*n + j*n + i];
                            current_masked = current_point > 0.5;

                            sum += current_point;
                            intersects += (current_masked) ^ (prev_masked);
                            prev_masked = current_masked;

                            k += d2;
                            j += d1;
                            i += d0;
                        }
                    }
                }
            }

        }
        directions_vectors_mil[v] = sum / intersects;
    }

}

#define BLOCK_SIZE 16
#define NUM_ACC 4
void dummy0(const float *hr_sphere_region, int n, float *directions_vectors_mil) {

    const int n_vectors = NUM_DIRECTIONS;
    float directions_vectors_bone_length[n_vectors];
    unsigned int directions_vectors_intercepts[n_vectors];
    int count = 0;

#if 1
    for (int kk = 0; kk < n; kk+=BLOCK_SIZE) {
        for (int jj = 0; jj < n; jj+=BLOCK_SIZE) {
            for (int ii = 0; ii < n; ii+=BLOCK_SIZE) {
                /* for every direction vector */
                for (int v = 0; v < n_vectors; ++v) {

                    /* Init accumulators */
                    float acc1 = 0.0, acc5 = 0.0;
                    float acc2 = 0.0, acc6 = 0.0;
                    float acc3 = 0.0, acc7 = 0.0;
                    float acc4 = 0.0, acc8 = 0.0;

                    unsigned int edge_count1 = 1, edge_count5 = 1;
                    unsigned int edge_count2 = 1, edge_count6 = 1;
                    unsigned int edge_count3 = 1, edge_count7 = 1;
                    unsigned int edge_count4 = 1, edge_count8 = 1;

                    for (int k = kk; k < kk + BLOCK_SIZE; k += STRIDE) {
                        for (int j = jj; j < jj + BLOCK_SIZE; j += STRIDE*NUM_ACC) {

                            unsigned int i_prev = (ii > 0) ? ii - 1 : 0;

                            unsigned int prev_mask1 = 0; //hr_sphere_region[ k*n*n + (j+0*STRIDE)*n + i_prev ]  > 0.5;
                            unsigned int prev_mask2 = 1; //hr_sphere_region[ k*n*n + (j+1*STRIDE)*n + i_prev ]  > 0.5;
                            unsigned int prev_mask3 = 0; //hr_sphere_region[ k*n*n + (j+2*STRIDE)*n + i_prev ]  > 0.5;
                            unsigned int prev_mask4 = 1; //hr_sphere_region[ k*n*n + (j+3*STRIDE)*n + i_prev ]  > 0.5;

                            for (int i = ii; i < ii + BLOCK_SIZE; ++i) {

                                unsigned int curr_mask1; //, curr_mask5;
                                unsigned int curr_mask2; //, curr_mask6;
                                unsigned int curr_mask3; //, curr_mask7;
                                unsigned int curr_mask4; //, curr_mask8;

                                /* Load working set */
                                float r1 = hr_sphere_region[ k*n*n + (j+0*STRIDE)*n + i];  //float r5 = *ptr_hr++;
                                float r2 = hr_sphere_region[ k*n*n + (j+1*STRIDE)*n + i];  //float r6 = *ptr_hr++;
                                float r3 = hr_sphere_region[ k*n*n + (j+2*STRIDE)*n + i];  //float r7 = *ptr_hr++;
                                float r4 = hr_sphere_region[ k*n*n + (j+3*STRIDE)*n + i];  //float r8 = *ptr_hr;

                                /* Perform computation */
                                acc1 += r1;  //acc5 += r5;
                                acc2 += r2;  //acc6 += r6;
                                acc3 += r3;  //acc7 += r7;
                                acc4 += r4;  //acc8 += r8;

                                /* Calculate masks */
                                curr_mask1 = r1 > 0.5;  //curr_mask5 = r5 > 0.5;
                                curr_mask2 = r2 > 0.5;  //curr_mask6 = r6 > 0.5;
                                curr_mask3 = r3 > 0.5;  //curr_mask7 = r7 > 0.5;
                                curr_mask4 = r4 > 0.5;  //curr_mask8 = r8 > 0.5;

                                /* Detect edge and add to counter */
                                edge_count1 += curr_mask1 ^ prev_mask1;  //edge_count5 += curr_mask5 ^ curr_mask4;
                                edge_count2 += curr_mask2 ^ prev_mask2;  //edge_count6 += curr_mask6 ^ curr_mask5;
                                edge_count3 += curr_mask3 ^ prev_mask3;  //edge_count7 += curr_mask7 ^ curr_mask6;
                                edge_count4 += curr_mask4 ^ prev_mask4;  //edge_count8 += curr_mask8 ^ curr_mask7;

                                /* Update state of prev_mask */
                                prev_mask1 = curr_mask1;
                                prev_mask2 = curr_mask2;
                                prev_mask3 = curr_mask3;
                                prev_mask4 = curr_mask4;
                            }
                        }
                    }

                    acc1 += acc2;   edge_count1 += edge_count2;
                    acc3 += acc4;   edge_count3 += edge_count4;
                    acc5 += acc6;   edge_count5 += edge_count6;
                    acc7 += acc8;   edge_count7 += edge_count8;
                    acc1 += acc3;   edge_count1 += edge_count3;
                    acc5 += acc7;   edge_count5 += edge_count7;
                    directions_vectors_bone_length[v] = acc1 + acc5;
                    directions_vectors_intercepts[v]  = edge_count1 + edge_count5;
                }
            }
        }
    }


    for (int v = 0; v < n_vectors; ++v) {
        directions_vectors_mil[v] = directions_vectors_bone_length[v] / directions_vectors_intercepts[v];
    }


#else

    /* for every direction vector */
    for (int v = 0; v < n_vectors; ++v) {
        double acc1 = 0.0;  double acc5 = 0.0;
        double acc2 = 0.0;  double acc6 = 0.0;
        double acc3 = 0.0;  double acc7 = 0.0;
        double acc4 = 0.0;  double acc8 = 0.0;

        for (int kk = 0; kk < n; kk += STRIDE) {
            for (int jj = 0; jj < n; jj += STRIDE) {
                for (int i = 0; i < n; i += NUM_ACC) {
                    acc1 += hr_sphere_region[kk * n * n + jj * n + i];
                    acc2 += hr_sphere_region[kk * n * n + jj * n + i + 1];
                    acc3 += hr_sphere_region[kk * n * n + jj * n + i + 2];
                    acc4 += hr_sphere_region[kk * n * n + jj * n + i + 3];
                    acc5 += hr_sphere_region[kk * n * n + jj * n + i + 4];
                    acc6 += hr_sphere_region[kk * n * n + jj * n + i + 5];
                    acc7 += hr_sphere_region[kk * n * n + jj * n + i + 6];
                    acc8 += hr_sphere_region[kk * n * n + jj * n + i + 7];
//                    count++;
                }
            }
        }
//
//        for (int i = 0; i < n*n*n/(STRIDE*STRIDE); i += NUM_ACC) {
//            acc1 += hr_sphere_region[i];
//            acc2 += hr_sphere_region[i + 1];
//            acc3 += hr_sphere_region[i + 2];
//            acc4 += hr_sphere_region[i + 3];
//            acc5 += hr_sphere_region[i + 4];
//            acc6 += hr_sphere_region[i + 5];
//            acc7 += hr_sphere_region[i + 6];
//            acc8 += hr_sphere_region[i + 7];
////            count++;
//        }

        acc1 += acc2;
        acc3 += acc4;
        acc5 += acc6;
        acc7 += acc8;
        acc1 += acc3;
        acc5 += acc7;
        directions_vectors_mil[v] = acc1 + acc5;
    }

//    printf("%d\n", count);
#endif

}

void working (int vecID, int* hr_sphere_region) {
    switch (vecID) {
//        case 1: /* Vector (1,0,0) */
//            /* Unroll over dimension y */
//            r1 = hr_sphere_region[ k*n*n + (j+0*STRIDE)*n + i];
//            r2 = hr_sphere_region[ k*n*n + (j+1*STRIDE)*n + i];
//            r3 = hr_sphere_region[ k*n*n + (j+2*STRIDE)*n + i];
//            r4 = hr_sphere_region[ k*n*n + (j+3*STRIDE)*n + i];
//            break;
//        case 2: /* Vector (0,1,0) */
//            /* Unroll over dimension x */
//            r1 = hr_sphere_region[ k*n*n + i*n + (j+0*STRIDE)];
//            r2 = hr_sphere_region[ k*n*n + i*n + (j+1*STRIDE)];
//            r3 = hr_sphere_region[ k*n*n + i*n + (j+2*STRIDE)];
//            r4 = hr_sphere_region[ k*n*n + i*n + (j+3*STRIDE)];
//            break;
//        case 3: /* Vector (0,0,1) */
//            /* Unroll over dimension x */
//            r1 = hr_sphere_region[ i*n*n + k*n + (j+0*STRIDE)];
//            r2 = hr_sphere_region[ i*n*n + k*n + (j+1*STRIDE)];
//            r3 = hr_sphere_region[ i*n*n + k*n + (j+2*STRIDE)];
//            r4 = hr_sphere_region[ i*n*n + k*n + (j+3*STRIDE)];
//            break;
//        case 4: /* Vector (1,1,0) */
//            /* Unroll over dimensions x,y */
//            r1 = hr_sphere_region[ k*n*n + j1*n + (i1+1+0*STRIDE) ];
//            r2 = hr_sphere_region[ k*n*n + (j2+1+0*STRIDE)*n + i2 ];
//            if (i1+1+1*STRIDE < ii + BLOCK_SIZE) {
//                r3 = hr_sphere_region[ k*n*n + j1*n + (i1+1+1*STRIDE) ];
//                r4 = hr_sphere_region[ k*n*n + (j2+1+1*STRIDE)*n + i2 ];
//            } else {
//                r3 = 0.0f;
//                r4 = 0.0f;
//            }
//            break;

//        r1 = hr_sphere_region[ k*n*n + j1*n + (i1-1-0*STRIDE) ];
//        r2 = hr_sphere_region[ k*n*n + (j2+1+0*STRIDE)*n + i2 ];
//        if (j2+1+1*STRIDE < jj + BLOCK_SIZE) {
//            r3 = hr_sphere_region[ k*n*n + j1*n + (i1-1-1*STRIDE) ];
//            r4 = hr_sphere_region[ k*n*n + (j2+1+1*STRIDE)*n + i2 ];
//        } else {
//            r3 = 0.0f;
//            r4 = 0.0f;
//        }

        default: ;
    }
}

#define LOAD_DATA_SET_1D                                            \
        switch (vecID) {                                            \
            case 1: /* Vector (1,0,0) */                            \
                /* Unroll over dimension y */                       \
                r1 = hr_sphere_region[ k*n*n + (j+0*STRIDE)*n + i]; \
                r2 = hr_sphere_region[ k*n*n + (j+1*STRIDE)*n + i]; \
                r3 = hr_sphere_region[ k*n*n + (j+2*STRIDE)*n + i]; \
                r4 = hr_sphere_region[ k*n*n + (j+3*STRIDE)*n + i]; \
                break;                                              \
            case 2: /* Vector (0,1,0) */                            \
                /* Unroll over dimension x */                       \
                r1 = hr_sphere_region[ k*n*n + i*n + (j+0*STRIDE)]; \
                r2 = hr_sphere_region[ k*n*n + i*n + (j+1*STRIDE)]; \
                r3 = hr_sphere_region[ k*n*n + i*n + (j+2*STRIDE)]; \
                r4 = hr_sphere_region[ k*n*n + i*n + (j+3*STRIDE)]; \
                break;                                              \
            case 3: /* Vector (0,0,1) */                            \
                /* Unroll over dimension x */                       \
                r1 = hr_sphere_region[ i*n*n + k*n + (j+0*STRIDE)]; \
                r2 = hr_sphere_region[ i*n*n + k*n + (j+1*STRIDE)]; \
                r3 = hr_sphere_region[ i*n*n + k*n + (j+2*STRIDE)]; \
                r4 = hr_sphere_region[ i*n*n + k*n + (j+3*STRIDE)]; \
                break;                                              \
            default: ;                                              \
        }

#define LOAD_DATA_SET_2D                                                 \
    switch (vecID) {                                                     \
        case 4: /* Vector (1,1,0) */                                     \
            /* Unroll over dimensions x,y */                             \
            r1 = hr_sphere_region[ k*n*n + j1*n + (i1+1+0*STRIDE) ];     \
            r2 = hr_sphere_region[ k*n*n + (j2+1+0*STRIDE)*n + i2 ];     \
            if (i1+1+1*STRIDE < ii + BLOCK_SIZE) {                       \
                r3 = hr_sphere_region[ k*n*n + j1*n + (i1+1+1*STRIDE) ]; \
                r4 = hr_sphere_region[ k*n*n + (j2+1+1*STRIDE)*n + i2 ]; \
            } else {                                                     \
                r3 = 0.0f;                                               \
                r4 = 0.0f;                                               \
            }                                                            \
            break;                                                       \
        case 5: /* Vector (1,0,1) */                                     \
            /* Unroll over dimensions x,z */                             \
            r1 = hr_sphere_region[ j1*n*n + k*n + (i1+1+0*STRIDE) ];     \
            r2 = hr_sphere_region[ (j2+1+0*STRIDE)*n*n + k*n + i2 ];     \
            if (i1+1+1*STRIDE < ii + BLOCK_SIZE) {                       \
                r3 = hr_sphere_region[ j1*n*n + k*n + (i1+1+1*STRIDE) ]; \
                r4 = hr_sphere_region[ (j2+1+1*STRIDE)*n*n + k*n + i2 ]; \
            } else {                                                     \
                r3 = 0.0f;                                               \
                r4 = 0.0f;                                               \
            }                                                            \
            break;                                                       \
        case 6: /* Vector (0,1,1) */                                     \
            /* Unroll over dimensions y,z */                             \
            r1 = hr_sphere_region[ (i1+1+0*STRIDE)*n*n + j1*n + k ];     \
            r2 = hr_sphere_region[ i2*n*n + (j2+1+0*STRIDE)*n + k ];     \
            if (i1+1+1*STRIDE < ii + BLOCK_SIZE) {                       \
                r3 = hr_sphere_region[ (i1+1+1*STRIDE)*n*n + j1*n + k ]; \
                r4 = hr_sphere_region[ i2*n*n + (j2+1+1*STRIDE)*n + k ]; \
            } else {                                                     \
                r3 = 0.0f;                                               \
                r4 = 0.0f;                                               \
            }                                                            \
            break;                                                       \
        case 7: /* Vector (-1,1,0) */                                    \
            /* Unroll over dimensions x,y */                             \
            r1 = hr_sphere_region[ k*n*n + j1*n + (i1-1-0*STRIDE) ];     \
            r2 = hr_sphere_region[ k*n*n + (j2+1+0*STRIDE)*n + i2 ];     \
            if (j2+1+1*STRIDE < jj + BLOCK_SIZE) {                       \
                r3 = hr_sphere_region[ k*n*n + j1*n + (i1-1-1*STRIDE) ]; \
                r4 = hr_sphere_region[ k*n*n + (j2+1+1*STRIDE)*n + i2 ]; \
            } else {                                                     \
                r3 = 0.0f;                                               \
                r4 = 0.0f;                                               \
            }                                                            \
            break;                                                       \
        case 8: /* Vector (-1,0,1) */                                    \
            /* Unroll over dimensions x,y */                             \
            r1 = hr_sphere_region[ j1*n*n + k*n + (i1-1-0*STRIDE) ];     \
            r2 = hr_sphere_region[ (j2+1+0*STRIDE)*n*n + k*n + i2 ];     \
            if (j2+1+1*STRIDE < jj + BLOCK_SIZE) {                       \
                r3 = hr_sphere_region[ j1*n*n + k*n + (i1-1-1*STRIDE) ]; \
                r4 = hr_sphere_region[ (j2+1+1*STRIDE)*n*n + k*n + i2 ]; \
            } else {                                                     \
                r3 = 0.0f;                                               \
                r4 = 0.0f;                                               \
            }                                                            \
            break;                                                       \
        case 9: /* Vector (0,1,-1) */                                    \
            /* Unroll over dimensions x,y */                             \
            r1 = hr_sphere_region[ (i1-1-0*STRIDE)*n*n + j1*n + k ];     \
            r2 = hr_sphere_region[ i2*n*n + (j2+1+0*STRIDE)*n + k ];     \
            if (j2+1+1*STRIDE < jj + BLOCK_SIZE) {                       \
                r3 = hr_sphere_region[ (i1-1-1*STRIDE)*n*n + j1*n + k ]; \
                r4 = hr_sphere_region[ i2*n*n + (j2+1+1*STRIDE)*n + k ]; \
            } else {                                                     \
                r3 = 0.0f;                                               \
                r4 = 0.0f;                                               \
            }                                                            \
            break;                                                       \
        default: ;                                                       \
    }


#define COMPUTATION \
        unsigned int curr_mask1;                \
        unsigned int curr_mask2;                \
        unsigned int curr_mask3;                \
        unsigned int curr_mask4;                \
        acc1 += r1;                             \
        acc2 += r2;                             \
        acc3 += r3;                             \
        acc4 += r4;                             \
                                                \
        /* Calculate masks */                   \
        curr_mask1 = r1 > 0.5f;                 \
        curr_mask2 = r2 > 0.5f;                 \
        curr_mask3 = r3 > 0.5f;                 \
        curr_mask4 = r4 > 0.5f;                 \
                                                \
        /* Detect edge and add to counter */    \
        edge_count1 += curr_mask1 ^ prev_mask1; \
        edge_count2 += curr_mask2 ^ prev_mask2; \
        edge_count3 += curr_mask3 ^ prev_mask3; \
        edge_count4 += curr_mask4 ^ prev_mask4;

///
/// Blocking version for vectors (1,0,0), (0,1,0), (0,0,1).
/// Using accumulators and unrolling.
///
inline float mil_1D(const float *hr_sphere_region, const int n, const int kk, const int jj, const int ii, const int vecID) {
    float bone_length;
    unsigned int intercepts;

    /* Init accumulators */
    float acc1 = 0.0, acc5 = 0.0;
    float acc2 = 0.0, acc6 = 0.0;
    float acc3 = 0.0, acc7 = 0.0;
    float acc4 = 0.0, acc8 = 0.0;

    unsigned int edge_count1 = 1, edge_count5 = 1;
    unsigned int edge_count2 = 1, edge_count6 = 1;
    unsigned int edge_count3 = 1, edge_count7 = 1;
    unsigned int edge_count4 = 1, edge_count8 = 1;

    for (int k = kk; k < kk + BLOCK_SIZE; k += STRIDE) {
        for (int j = jj; j < jj + BLOCK_SIZE; j += STRIDE*NUM_ACC) {

            //unsigned int i_prev = (ii > 0) ? ii - 1 : 0;

            unsigned int prev_mask1 = 0; //hr_sphere_region[ k*n*n + (j+0*STRIDE)*n + i_prev ]  > 0.5;
            unsigned int prev_mask2 = 0; //hr_sphere_region[ k*n*n + (j+1*STRIDE)*n + i_prev ]  > 0.5;
            unsigned int prev_mask3 = 0; //hr_sphere_region[ k*n*n + (j+2*STRIDE)*n + i_prev ]  > 0.5;
            unsigned int prev_mask4 = 0; //hr_sphere_region[ k*n*n + (j+3*STRIDE)*n + i_prev ]  > 0.5;

            for (int i = ii; i < ii + BLOCK_SIZE; ++i) {
                float r1, r2, r3, r4;

                /* Load working set */
                LOAD_DATA_SET_1D

                /* Perform computation */
                COMPUTATION

                /* Update state of prev_mask */
                prev_mask1 = curr_mask1;
                prev_mask2 = curr_mask2;
                prev_mask3 = curr_mask3;
                prev_mask4 = curr_mask4;
#ifdef  DEBUG
                count++;
#endif
            }
        }
    }

    acc1 += acc2;   edge_count1 += edge_count2;
    acc3 += acc4;   edge_count3 += edge_count4;
    acc5 += acc6;   edge_count5 += edge_count6;
    acc7 += acc8;   edge_count7 += edge_count8;
    acc1 += acc3;   edge_count1 += edge_count3;
    acc5 += acc7;   edge_count5 += edge_count7;
    bone_length = acc1 + acc5;
    intercepts  = edge_count1 + edge_count5;

#ifdef  DEBUG
    if (!already_tested[vecID]) {
        printf("%d\n", count);
        already_tested[vecID] = 1;
    }
    count = 0;
#endif

    return bone_length / intercepts;

}

///
/// Blocking version for vectors (1,1,0), (0,1,1), (1,0,1).
/// Using accumulators and unrolling.
///
inline float mil_2D_pos(const float *hr_sphere_region, const int n, const int kk, const int jj, const int ii, const int vecID) {
    float bone_length;
    unsigned int intercepts;

    /* Init accumulators */
    float acc1 = 0.0, acc5 = 0.0;
    float acc2 = 0.0, acc6 = 0.0;
    float acc3 = 0.0, acc7 = 0.0;
    float acc4 = 0.0, acc8 = 0.0;

    unsigned int edge_count1 = 1, edge_count5 = 1;
    unsigned int edge_count2 = 1, edge_count6 = 1;
    unsigned int edge_count3 = 1, edge_count7 = 1;
    unsigned int edge_count4 = 1, edge_count8 = 1;

    for (int k = kk; k < kk + BLOCK_SIZE; k += STRIDE) {
        for (int ij = 0; ij < BLOCK_SIZE; ij += STRIDE*NUM_ACC/2) {

            unsigned int prev_mask1 = 0; // change this to previous of the initial val
            unsigned int prev_mask2 = 0;
            unsigned int prev_mask3 = 0;
            unsigned int prev_mask4 = 0;
            int i1 = ii + ij;
            int j1 = jj;
            int i2 = ii;
            int j2 = jj + ij;
            while (i1 + 1 < ii + BLOCK_SIZE /*&& j2+1 < jj + BLOCK_SIZE - 1*/) {
                float r1, r2, r3, r4;

                /* Load working set */
                LOAD_DATA_SET_2D

                /* Perform computation */
                COMPUTATION

                /* Update state of prev_mask */
                prev_mask1 = curr_mask1;
                prev_mask2 = curr_mask2;
                if (i1+1+1*STRIDE < ii + BLOCK_SIZE) {
                    prev_mask3 = curr_mask3;
                    prev_mask4 = curr_mask4;
                }
                count++;


                ++i1;
                ++j1;
                ++i2;
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
    intercepts  = edge_count1 + edge_count5;

#ifdef  DEBUG
    if (!already_tested[vecID]) {
        printf("%d\n", count);
        already_tested[vecID] = 1;
    }
    count = 0;
#endif

    return bone_length / intercepts;

}

///
/// Blocking version for vectors (-1,1,0), (0,-1,1), (-1,0,1).
/// Using accumulators and unrolling.
///
inline float mil_2D_neg(const float *hr_sphere_region, const int n, const int kk, const int jj, const int ii, const int vecID) {
    float bone_length;
    unsigned int intercepts;

    /* Init accumulators */
    float acc1 = 0.0, acc5 = 0.0;
    float acc2 = 0.0, acc6 = 0.0;
    float acc3 = 0.0, acc7 = 0.0;
    float acc4 = 0.0, acc8 = 0.0;

    unsigned int edge_count1 = 1, edge_count5 = 1;
    unsigned int edge_count2 = 1, edge_count6 = 1;
    unsigned int edge_count3 = 1, edge_count7 = 1;
    unsigned int edge_count4 = 1, edge_count8 = 1;

    for (int k = kk; k < kk + BLOCK_SIZE; k += STRIDE) {
        for (int ij = 0; ij < BLOCK_SIZE; ij += STRIDE*NUM_ACC/2) {

            unsigned int prev_mask1 = 0; // change this to previous of the initial val
            unsigned int prev_mask2 = 0;
            unsigned int prev_mask3 = 0;
            unsigned int prev_mask4 = 0;
            int i1 = ii + (BLOCK_SIZE-1) - ij;  /* Start at the end of the row in the block */
            int j1 = jj;
            int i2 = ii + (BLOCK_SIZE-1);
            int j2 = jj + ij;
            while (j2 + 1 < jj + BLOCK_SIZE) {
                float r1, r2, r3, r4;

                /* Load working set */
                LOAD_DATA_SET_2D

                /* Perform computation */
                COMPUTATION

                /* Update state of prev_mask */
                prev_mask1 = curr_mask1;
                prev_mask2 = curr_mask2;
                if (j2+1+1*STRIDE < jj + BLOCK_SIZE) {
                    prev_mask3 = curr_mask3;
                    prev_mask4 = curr_mask4;
                }
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
    intercepts  = edge_count1 + edge_count5;

#ifdef  DEBUG
    if (!already_tested[vecID]) {
        printf("%d\n", count);
        already_tested[vecID] = 1;
    }
    count = 0;
#endif

    return bone_length / intercepts;

}

///
/// Test for first vector only (1,0,0).
///
void mil_test_v1(const float *hr_sphere_region, int n, float *directions_vectors_mil) {

    const int n_vectors = NUM_DIRECTIONS;
    for (int kk = 0; kk < n; kk+=BLOCK_SIZE) {
        for (int jj = 0; jj < n; jj+=BLOCK_SIZE) {
            for (int ii = 0; ii < n; ii+=BLOCK_SIZE) {
                /* for every direction vector */
                for (int v = 0; v < n_vectors; ++v) {

                    directions_vectors_mil[v] = mil_1D(hr_sphere_region, n, kk, jj, ii, 1);

                }
            }
        }
    }
}

///
/// Test second vector (0,1,0).
///
void mil_test_v2(const float *hr_sphere_region, int n, float *directions_vectors_mil) {

    const int n_vectors = NUM_DIRECTIONS;
    for (int kk = 0; kk < n; kk+=BLOCK_SIZE) {
        for (int jj = 0; jj < n; jj+=BLOCK_SIZE) {
            for (int ii = 0; ii < n; ii+=BLOCK_SIZE) {
                /* for every direction vector */
                for (int v = 0; v < n_vectors; ++v) {

                    directions_vectors_mil[v] = mil_1D(hr_sphere_region, n, kk, jj, ii, 2);

                }
            }
        }
    }

}

///
/// Test third vector (0,0,1).
///
void mil_test_v3(const float *hr_sphere_region, int n, float *directions_vectors_mil) {

    const int n_vectors = NUM_DIRECTIONS;
    for (int kk = 0; kk < n; kk+=BLOCK_SIZE) {
        for (int jj = 0; jj < n; jj+=BLOCK_SIZE) {
            for (int ii = 0; ii < n; ii+=BLOCK_SIZE) {
                /* for every direction vector */
                for (int v = 0; v < n_vectors; ++v) {

                    directions_vectors_mil[v] = mil_1D(hr_sphere_region, n, kk, jj, ii, 3);

                }
            }
        }
    }

}

///
/// Test third vector (1,1,0).
///
void mil_test_v4(const float *hr_sphere_region, int n, float *directions_vectors_mil) {

    const int n_vectors = NUM_DIRECTIONS;
    for (int kk = 0; kk < n; kk+=BLOCK_SIZE) {
        for (int jj = 0; jj < n; jj+=BLOCK_SIZE) {
            for (int ii = 0; ii < n; ii+=BLOCK_SIZE) {
                /* for every direction vector */
                for (int v = 0; v < n_vectors; ++v) {

                    directions_vectors_mil[v] = mil_2D_pos(hr_sphere_region, n, kk, jj, ii, 4);

                }
            }
        }
    }

}

///
/// Test third vector (1,0,1).
///
void mil_test_v5(const float *hr_sphere_region, int n, float *directions_vectors_mil) {

    const int n_vectors = NUM_DIRECTIONS;
    for (int kk = 0; kk < n; kk+=BLOCK_SIZE) {
        for (int jj = 0; jj < n; jj+=BLOCK_SIZE) {
            for (int ii = 0; ii < n; ii+=BLOCK_SIZE) {
                /* for every direction vector */
                for (int v = 0; v < n_vectors; ++v) {

                    directions_vectors_mil[v] = mil_2D_pos(hr_sphere_region, n, kk, jj, ii, 5);

                }
            }
        }
    }

}

///
/// Test third vector (0,1,1).
///
void mil_test_v6(const float *hr_sphere_region, int n, float *directions_vectors_mil) {

    const int n_vectors = NUM_DIRECTIONS;
    for (int kk = 0; kk < n; kk+=BLOCK_SIZE) {
        for (int jj = 0; jj < n; jj+=BLOCK_SIZE) {
            for (int ii = 0; ii < n; ii+=BLOCK_SIZE) {
                /* for every direction vector */
                for (int v = 0; v < n_vectors; ++v) {

                    directions_vectors_mil[v] = mil_2D_pos(hr_sphere_region, n, kk, jj, ii, 6);

                }
            }
        }
    }

}

///
/// Test third vector (-1,1,0).
///
void mil_test_v7(const float *hr_sphere_region, int n, float *directions_vectors_mil) {

    const int n_vectors = NUM_DIRECTIONS;
    for (int kk = 0; kk < n; kk+=BLOCK_SIZE) {
        for (int jj = 0; jj < n; jj+=BLOCK_SIZE) {
            for (int ii = 0; ii < n; ii+=BLOCK_SIZE) {
                /* for every direction vector */
                for (int v = 0; v < n_vectors; ++v) {

                    directions_vectors_mil[v] = mil_2D_neg(hr_sphere_region, n, kk, jj, ii, 7);

                }
            }
        }
    }

}

///
/// Test third vector (-1,0,1).
///
void mil_test_v8(const float *hr_sphere_region, int n, float *directions_vectors_mil) {

    const int n_vectors = NUM_DIRECTIONS;
    for (int kk = 0; kk < n; kk+=BLOCK_SIZE) {
        for (int jj = 0; jj < n; jj+=BLOCK_SIZE) {
            for (int ii = 0; ii < n; ii+=BLOCK_SIZE) {
                /* for every direction vector */
                for (int v = 0; v < n_vectors; ++v) {

                    directions_vectors_mil[v] = mil_2D_neg(hr_sphere_region, n, kk, jj, ii, 8);

                }
            }
        }
    }

}

///
/// Test third vector (0,1,-1).
///
void mil_test_v9(const float *hr_sphere_region, int n, float *directions_vectors_mil) {

    const int n_vectors = NUM_DIRECTIONS;
    for (int kk = 0; kk < n; kk+=BLOCK_SIZE) {
        for (int jj = 0; jj < n; jj+=BLOCK_SIZE) {
            for (int ii = 0; ii < n; ii+=BLOCK_SIZE) {
                /* for every direction vector */
                for (int v = 0; v < n_vectors; ++v) {

                    directions_vectors_mil[v] = mil_2D_neg(hr_sphere_region, n, kk, jj, ii, 9);

                }
            }
        }
    }

}

///
/// Blocking version for vector (0,1,0).
/// Using accumulators and unrolling over dimension x.
///
void mil_block_v2(const float *hr_sphere_region, int n, float *directions_vectors_mil) {

    const int n_vectors = NUM_DIRECTIONS;
    float directions_vectors_bone_length[n_vectors];
    unsigned int directions_vectors_intercepts[n_vectors];
    unsigned int count = 0;
    static unsigned int flag = 0;

    for (int kk = 0; kk < n; kk+=BLOCK_SIZE) {
        for (int jj = 0; jj < n; jj+=BLOCK_SIZE) {
            for (int ii = 0; ii < n; ii+=BLOCK_SIZE) {
                /* for every direction vector */
                for (int v = 0; v < n_vectors; ++v) {

                    /* Init accumulators */
                    float acc1 = 0.0, acc5 = 0.0;
                    float acc2 = 0.0, acc6 = 0.0;
                    float acc3 = 0.0, acc7 = 0.0;
                    float acc4 = 0.0, acc8 = 0.0;

                    unsigned int edge_count1 = 1, edge_count5 = 1;
                    unsigned int edge_count2 = 1, edge_count6 = 1;
                    unsigned int edge_count3 = 1, edge_count7 = 1;
                    unsigned int edge_count4 = 1, edge_count8 = 1;

                    for (int k = kk; k < kk + BLOCK_SIZE; k += STRIDE) {
                        for (int i = ii; i < ii + BLOCK_SIZE; i += STRIDE*NUM_ACC) {

                            //unsigned int i_prev = (ii > 0) ? ii - 1 : 0;

                            unsigned int prev_mask1 = 0; //hr_sphere_region[ k*n*n + (j+0*STRIDE)*n + i_prev ]  > 0.5;
                            unsigned int prev_mask2 = 0; //hr_sphere_region[ k*n*n + (j+1*STRIDE)*n + i_prev ]  > 0.5;
                            unsigned int prev_mask3 = 0; //hr_sphere_region[ k*n*n + (j+2*STRIDE)*n + i_prev ]  > 0.5;
                            unsigned int prev_mask4 = 0; //hr_sphere_region[ k*n*n + (j+3*STRIDE)*n + i_prev ]  > 0.5;

                            for (int j = jj; j < jj + BLOCK_SIZE; ++j) {

                                /* Load working set */
                                float r1 = hr_sphere_region[ k*n*n + j*n + i+0*STRIDE];  //float r5 = *ptr_hr++;
                                float r2 = hr_sphere_region[ k*n*n + j*n + i+1*STRIDE];  //float r6 = *ptr_hr++;
                                float r3 = hr_sphere_region[ k*n*n + j*n + i+2*STRIDE];  //float r7 = *ptr_hr++;
                                float r4 = hr_sphere_region[ k*n*n + j*n + i+3*STRIDE];  //float r8 = *ptr_hr;

                                /* Perform computation */
                                COMPUTATION

                                /* Update state of prev_mask */
                                prev_mask1 = curr_mask1;
                                prev_mask2 = curr_mask2;
                                prev_mask3 = curr_mask3;
                                prev_mask4 = curr_mask4;
                                count++;
                            }
                        }
                    }

                    acc1 += acc2;   edge_count1 += edge_count2;
                    acc3 += acc4;   edge_count3 += edge_count4;
                    acc5 += acc6;   edge_count5 += edge_count6;
                    acc7 += acc8;   edge_count7 += edge_count8;
                    acc1 += acc3;   edge_count1 += edge_count3;
                    acc5 += acc7;   edge_count5 += edge_count7;
                    directions_vectors_bone_length[v] = acc1 + acc5;
                    directions_vectors_intercepts[v]  = edge_count1 + edge_count5;
                }
            }
        }
    }

    if (!flag) {
        printf("%d\n", count);
        flag = 1;
    }

    for (int v = 0; v < n_vectors; ++v) {
        directions_vectors_mil[v] = directions_vectors_bone_length[v] / directions_vectors_intercepts[v];
    }

}

///
/// Blocking version for vector (0,0,1).
/// Using accumulators and unrolling over dimension x.
///
void mil_block_v3(const float *hr_sphere_region, int n, float *directions_vectors_mil) {

    const int n_vectors = NUM_DIRECTIONS;
    float directions_vectors_bone_length[n_vectors];
    unsigned int directions_vectors_intercepts[n_vectors];
    unsigned int count = 0;
    static unsigned int flag = 0;

    for (int kk = 0; kk < n; kk+=BLOCK_SIZE) {
        for (int jj = 0; jj < n; jj+=BLOCK_SIZE) {
            for (int ii = 0; ii < n; ii+=BLOCK_SIZE) {
                /* for every direction vector */
                for (int v = 0; v < n_vectors; ++v) {

                    /* Init accumulators */
                    float acc1 = 0.0f, acc5 = 0.0f;
                    float acc2 = 0.0f, acc6 = 0.0f;
                    float acc3 = 0.0f, acc7 = 0.0f;
                    float acc4 = 0.0f, acc8 = 0.0f;

                    unsigned int edge_count1 = 1, edge_count5 = 1;
                    unsigned int edge_count2 = 1, edge_count6 = 1;
                    unsigned int edge_count3 = 1, edge_count7 = 1;
                    unsigned int edge_count4 = 1, edge_count8 = 1;

                    for (int j = jj; j < jj + BLOCK_SIZE; j += STRIDE) {
                        for (int i = ii; i < ii + BLOCK_SIZE; i += STRIDE*NUM_ACC) {

                            unsigned int prev_mask1 = 0; //hr_sphere_region[ k*n*n + (j+0*STRIDE)*n + i_prev ]  > 0.5;
                            unsigned int prev_mask2 = 0; //hr_sphere_region[ k*n*n + (j+1*STRIDE)*n + i_prev ]  > 0.5;
                            unsigned int prev_mask3 = 0; //hr_sphere_region[ k*n*n + (j+2*STRIDE)*n + i_prev ]  > 0.5;
                            unsigned int prev_mask4 = 0; //hr_sphere_region[ k*n*n + (j+3*STRIDE)*n + i_prev ]  > 0.5;

                            for (int k = kk; k < kk + BLOCK_SIZE; ++k) {

                                /* Load working set */
                                float r1 = hr_sphere_region[ k*n*n + j*n + i+0*STRIDE];  //float r5 = *ptr_hr++;
                                float r2 = hr_sphere_region[ k*n*n + j*n + i+1*STRIDE];  //float r6 = *ptr_hr++;
                                float r3 = hr_sphere_region[ k*n*n + j*n + i+2*STRIDE];  //float r7 = *ptr_hr++;
                                float r4 = hr_sphere_region[ k*n*n + j*n + i+3*STRIDE];  //float r8 = *ptr_hr;

                                /* Perform computation */
                                COMPUTATION

                                /* Update state of prev_mask */
                                prev_mask1 = curr_mask1;
                                prev_mask2 = curr_mask2;
                                prev_mask3 = curr_mask3;
                                prev_mask4 = curr_mask4;
                                count++;
                            }
                        }
                    }

                    acc1 += acc2;   edge_count1 += edge_count2;
                    acc3 += acc4;   edge_count3 += edge_count4;
                    acc5 += acc6;   edge_count5 += edge_count6;
                    acc7 += acc8;   edge_count7 += edge_count8;
                    acc1 += acc3;   edge_count1 += edge_count3;
                    acc5 += acc7;   edge_count5 += edge_count7;
                    directions_vectors_bone_length[v] = acc1 + acc5;
                    directions_vectors_intercepts[v]  = edge_count1 + edge_count5;
                }
            }
        }
    }

    if (!flag) {
        printf("%d\n", count);
        flag = 1;
    }

    for (int v = 0; v < n_vectors; ++v) {
        directions_vectors_mil[v] = directions_vectors_bone_length[v] / directions_vectors_intercepts[v];
    }

}

///
/// Blocking version for vector (1,1,0)
/// Using accumulators and unrolling over dimensions x,y.
///
void mil_block_v4(const float *hr_sphere_region, int n, float *directions_vectors_mil) {

    const int n_vectors = NUM_DIRECTIONS;
    float directions_vectors_bone_length[n_vectors];
    unsigned int directions_vectors_intercepts[n_vectors];
    unsigned int count = 0;
    static unsigned int flag = 0;

    for (int kk = 0; kk < n; kk+=BLOCK_SIZE) {
        for (int jj = 0; jj < n; jj+=BLOCK_SIZE) {
            for (int ii = 0; ii < n; ii+=BLOCK_SIZE) {
                /* for every direction vector */
                for (int v = 0; v < n_vectors; ++v) {

                    /* Init accumulators */
                    float acc1 = 0.0, acc5 = 0.0;
                    float acc2 = 0.0, acc6 = 0.0;
                    float acc3 = 0.0, acc7 = 0.0;
                    float acc4 = 0.0, acc8 = 0.0;

                    unsigned int edge_count1 = 1, edge_count5 = 1;
                    unsigned int edge_count2 = 1, edge_count6 = 1;
                    unsigned int edge_count3 = 1, edge_count7 = 1;
                    unsigned int edge_count4 = 1, edge_count8 = 1;

                    for (int k = kk; k < kk + BLOCK_SIZE; k += STRIDE) {
                        for (int ij = 0; ij < BLOCK_SIZE; ij += STRIDE*NUM_ACC/2) {

                            unsigned int prev_mask1 = 0; // change this to previous of the initial val
                            unsigned int prev_mask2 = 0; // change this to previous of the initial val
                            unsigned int prev_mask3 = 0; // change this to previous of the initial val
                            unsigned int prev_mask4 = 0; // change this to previous of the initial val
                            int i1 = ii + ij;
                            int j1 = jj;
                            int i2 = ii;
                            int j2 = jj + ij;
                            while (i1 + 1 < ii + BLOCK_SIZE /*&& j2+1 < jj + BLOCK_SIZE - 1*/) {

                                /* Load working set */
                                float r1, r2, r3, r4;
                                r1 = hr_sphere_region[ k*n*n + j1*n + (i1+1+0*STRIDE) ];
                                r2 = hr_sphere_region[ k*n*n + (j2+1+0*STRIDE)*n + i2 ];
                                if (i1+1+1*STRIDE < ii + BLOCK_SIZE) {
                                    r3 = hr_sphere_region[ k*n*n + j1*n + (i1+1+1*STRIDE) ];
                                    r4 = hr_sphere_region[ k*n*n + (j2+1+1*STRIDE)*n + i2 ];
                                } else {
                                    r3 = 0.0f;
                                    r4 = 0.0f;
                                }

                                /* Perform computation */
                                COMPUTATION

                                /* Update state of prev_mask */
                                prev_mask1 = curr_mask1;
                                prev_mask2 = curr_mask2;
                                if (i1+1+1*STRIDE < ii + BLOCK_SIZE) {
                                    prev_mask3 = curr_mask3;
                                    prev_mask4 = curr_mask4;
                                }
                                count++;


                                ++i1;
                                ++j1;
                                ++i2;
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
                    directions_vectors_bone_length[v] = acc1 + acc5;
                    directions_vectors_intercepts[v]  = edge_count1 + edge_count5;

                } /* Vector iteration*/

            } /* Blocking */
        }
    }

    if (!flag) {
        printf("%d\n", count);
        flag = 1;
    }

    for (int v = 0; v < n_vectors; ++v) {
        directions_vectors_mil[v] = directions_vectors_bone_length[v] / directions_vectors_intercepts[v];
    }

}

///
/// Blocking version for vector (-1,1,0)
/// Using accumulators and unrolling over dimensions x,y.
///
void mil_block_v5(const float *hr_sphere_region, int n, float *directions_vectors_mil) {

    const int n_vectors = NUM_DIRECTIONS;
    float directions_vectors_bone_length[n_vectors];
    unsigned int directions_vectors_intercepts[n_vectors];
    unsigned int count = 0;
    static unsigned int flag = 0;

    for (int kk = 0; kk < n; kk+=BLOCK_SIZE) {
        for (int jj = 0; jj < n; jj+=BLOCK_SIZE) {
            for (int ii = 0; ii < n; ii+=BLOCK_SIZE) {
                /* for every direction vector */
                for (int v = 0; v < n_vectors; ++v) {

                    /* Init accumulators */
                    float acc1 = 0.0, acc5 = 0.0;
                    float acc2 = 0.0, acc6 = 0.0;
                    float acc3 = 0.0, acc7 = 0.0;
                    float acc4 = 0.0, acc8 = 0.0;

                    unsigned int edge_count1 = 1, edge_count5 = 1;
                    unsigned int edge_count2 = 1, edge_count6 = 1;
                    unsigned int edge_count3 = 1, edge_count7 = 1;
                    unsigned int edge_count4 = 1, edge_count8 = 1;

                    for (int k = kk; k < kk + BLOCK_SIZE; k += STRIDE) {
                        for (int ij = 0; ij < BLOCK_SIZE; ij += STRIDE*NUM_ACC/2) {

                            unsigned int prev_mask1 = 0; // change this to previous of the initial val
                            unsigned int prev_mask2 = 0;
                            unsigned int prev_mask3 = 0;
                            unsigned int prev_mask4 = 0;
                            int i1 = ii + (BLOCK_SIZE-1) - ij;  /* Start at the end of the row in the block */
                            int j1 = jj;
                            int i2 = ii + (BLOCK_SIZE-1);
                            int j2 = jj + ij;
                            while (j2 + 1 < jj + BLOCK_SIZE) {

                                /* Load working set */
                                float r1, r2, r3, r4;
                                r1 = hr_sphere_region[ k*n*n + j1*n + (i1-1-0*STRIDE) ];
                                r2 = hr_sphere_region[ k*n*n + (j2+1+0*STRIDE)*n + i2 ];
                                if (j2+1+1*STRIDE < jj + BLOCK_SIZE) {
                                    r3 = hr_sphere_region[ k*n*n + j1*n + (i1-1-1*STRIDE) ];
                                    r4 = hr_sphere_region[ k*n*n + (j2+1+1*STRIDE)*n + i2 ];
                                } else {
                                    r3 = 0.0f;
                                    r4 = 0.0f;
                                }

                                /* Perform computation */
                                COMPUTATION

                                /* Update state of prev_mask */
                                prev_mask1 = curr_mask1;
                                prev_mask2 = curr_mask2;
                                if (j2+1+1*STRIDE < jj + BLOCK_SIZE) {
                                    prev_mask3 = curr_mask3;
                                    prev_mask4 = curr_mask4;
                                }
                                count++;


                                ++i1;
                                ++j1;
                                ++i2;
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
                    directions_vectors_bone_length[v] = acc1 + acc5;
                    directions_vectors_intercepts[v]  = edge_count1 + edge_count5;

                } /* Vector iteration*/

            } /* Blocking */
        }
    }

    if (!flag) {
        printf("%d\n", count);
        flag = 1;
    }

    for (int v = 0; v < n_vectors; ++v) {
        directions_vectors_mil[v] = directions_vectors_bone_length[v] / directions_vectors_intercepts[v];
    }

}

///
/// Blocking version for vector (1,0,1)
/// Using accumulators and unrolling over dimensions x,z.
///
void mil_block_v6(const float *hr_sphere_region, int n, float *directions_vectors_mil) {

    const int n_vectors = NUM_DIRECTIONS;
    float directions_vectors_bone_length[n_vectors];
    unsigned int directions_vectors_intercepts[n_vectors];
    unsigned int count = 0;
    static unsigned int flag = 0;

    for (int kk = 0; kk < n; kk+=BLOCK_SIZE) {
        for (int jj = 0; jj < n; jj+=BLOCK_SIZE) {
            for (int ii = 0; ii < n; ii+=BLOCK_SIZE) {
                /* for every direction vector */
                for (int v = 0; v < n_vectors; ++v) {

                    /* Init accumulators */
                    float acc1 = 0.0, acc5 = 0.0;
                    float acc2 = 0.0, acc6 = 0.0;
                    float acc3 = 0.0, acc7 = 0.0;
                    float acc4 = 0.0, acc8 = 0.0;

                    unsigned int edge_count1 = 1, edge_count5 = 1;
                    unsigned int edge_count2 = 1, edge_count6 = 1;
                    unsigned int edge_count3 = 1, edge_count7 = 1;
                    unsigned int edge_count4 = 1, edge_count8 = 1;

                    for (int j = jj; j < jj + BLOCK_SIZE; j += STRIDE) {
                        for (int ik = 0; ik < BLOCK_SIZE; ik += STRIDE*NUM_ACC/2) {

                            unsigned int prev_mask1 = 0; // change this to previous of the initial val
                            unsigned int prev_mask2 = 0;
                            unsigned int prev_mask3 = 0;
                            unsigned int prev_mask4 = 0;
                            int i1 = ii + ik;
                            int k1 = kk;
                            int i2 = ii;
                            int k2 = kk + ik;
                            while (i1 + 1 < ii + BLOCK_SIZE /*&& j2+1 < jj + BLOCK_SIZE - 1*/) {

                                /* Load working set */
                                float r1, r2, r3, r4;
                                r1 = hr_sphere_region[ k1*n*n + j*n + (i1+1+0*STRIDE) ];
                                r2 = hr_sphere_region[ (k2+1+0*STRIDE)*n*n + j*n + i2 ];
                                if (i1+1+1*STRIDE < ii + BLOCK_SIZE) {
                                    r3 = hr_sphere_region[ k1*n*n + j*n + (i1+1+1*STRIDE) ];
                                    r4 = hr_sphere_region[ (k2+1+1*STRIDE)*n*n + j*n + i2 ];
                                } else {
                                    r3 = 0.0f;
                                    r4 = 0.0f;
                                }

                                /* Perform computation */
                                COMPUTATION

                                /* Update state of prev_mask */
                                prev_mask1 = curr_mask1;
                                prev_mask2 = curr_mask2;
                                if (i1+1+1*STRIDE < ii + BLOCK_SIZE) {
                                    prev_mask3 = curr_mask3;
                                    prev_mask4 = curr_mask4;
                                }
                                count++;

                                ++i1;
                                ++k1;
                                ++i2;
                                ++k2;
                            }
                        }
                    } /* End iteration over dimension k */

                    acc1 += acc2;   edge_count1 += edge_count2;
                    acc3 += acc4;   edge_count3 += edge_count4;
                    acc5 += acc6;   edge_count5 += edge_count6;
                    acc7 += acc8;   edge_count7 += edge_count8;
                    acc1 += acc3;   edge_count1 += edge_count3;
                    acc5 += acc7;   edge_count5 += edge_count7;
                    directions_vectors_bone_length[v] = acc1 + acc5;
                    directions_vectors_intercepts[v]  = edge_count1 + edge_count5;

                } /* Vector iteration*/

            } /* Blocking */
        }
    }

    if (!flag) {
        printf("%d\n", count);
        flag = 1;
    }

    for (int v = 0; v < n_vectors; ++v) {
        directions_vectors_mil[v] = directions_vectors_bone_length[v] / directions_vectors_intercepts[v];
    }

}

///
/// Blocking version for vector (-1,0,1)
/// Using accumulators and unrolling over dimensions x,z.
///
void mil_block_v7(const float *hr_sphere_region, int n, float *directions_vectors_mil) {

    const int n_vectors = NUM_DIRECTIONS;
    float directions_vectors_bone_length[n_vectors];
    unsigned int directions_vectors_intercepts[n_vectors];
    unsigned int count = 0;
    static unsigned int flag = 0;

    for (int kk = 0; kk < n; kk+=BLOCK_SIZE) {
        for (int jj = 0; jj < n; jj+=BLOCK_SIZE) {
            for (int ii = 0; ii < n; ii+=BLOCK_SIZE) {
                /* for every direction vector */
                for (int v = 0; v < n_vectors; ++v) {

                    /* Init accumulators */
                    float acc1 = 0.0, acc5 = 0.0;
                    float acc2 = 0.0, acc6 = 0.0;
                    float acc3 = 0.0, acc7 = 0.0;
                    float acc4 = 0.0, acc8 = 0.0;

                    unsigned int edge_count1 = 1, edge_count5 = 1;
                    unsigned int edge_count2 = 1, edge_count6 = 1;
                    unsigned int edge_count3 = 1, edge_count7 = 1;
                    unsigned int edge_count4 = 1, edge_count8 = 1;

                    for (int j = jj; j < jj + BLOCK_SIZE; j += STRIDE) {
                        for (int ik = 0; ik < BLOCK_SIZE; ik += STRIDE*NUM_ACC/2) {

                            unsigned int prev_mask1 = 0; // change this to previous of the initial val
                            unsigned int prev_mask2 = 0;
                            unsigned int prev_mask3 = 0;
                            unsigned int prev_mask4 = 0;
                            int i1 = ii + (BLOCK_SIZE-1) - ik;  /* Start at the end of the row in the block */
                            int k1 = kk;
                            int i2 = ii + (BLOCK_SIZE-1);
                            int k2 = kk + ik;
                            while (k2 + 1 < kk + BLOCK_SIZE) {

                                /* Load working set */
                                float r1, r2, r3, r4;
                                r1 = hr_sphere_region[ k1*n*n + j*n + (i1-1-0*STRIDE) ];
                                r2 = hr_sphere_region[ (k2+1+0*STRIDE)*n*n + j*n + i2 ];
                                if (k2+1+1*STRIDE < kk + BLOCK_SIZE) {
                                    r3 = hr_sphere_region[ k1*n*n + j*n + (i1-1-1*STRIDE) ];
                                    r4 = hr_sphere_region[ (k2+1+1*STRIDE)*n*n + j*n + i2 ];
                                } else {
                                    r3 = 0.0f;
                                    r4 = 0.0f;
                                }

                                /* Perform computation */
                                COMPUTATION

                                /* Update state of prev_mask */
                                prev_mask1 = curr_mask1;
                                prev_mask2 = curr_mask2;
                                if (k2+1+1*STRIDE < kk + BLOCK_SIZE) {
                                    prev_mask3 = curr_mask3;
                                    prev_mask4 = curr_mask4;
                                }
                                count++;

                                ++i1;
                                ++k1;
                                ++i2;
                                ++k2;
                            }
                        }
                    } /* End iteration over dimension k */

                    acc1 += acc2;   edge_count1 += edge_count2;
                    acc3 += acc4;   edge_count3 += edge_count4;
                    acc5 += acc6;   edge_count5 += edge_count6;
                    acc7 += acc8;   edge_count7 += edge_count8;
                    acc1 += acc3;   edge_count1 += edge_count3;
                    acc5 += acc7;   edge_count5 += edge_count7;
                    directions_vectors_bone_length[v] = acc1 + acc5;
                    directions_vectors_intercepts[v]  = edge_count1 + edge_count5;

                } /* Vector iteration*/

            } /* Blocking */
        }
    }

    if (!flag) {
        printf("%d\n", count);
        flag = 1;
    }

    for (int v = 0; v < n_vectors; ++v) {
        directions_vectors_mil[v] = directions_vectors_bone_length[v] / directions_vectors_intercepts[v];
    }

}


///
/// Blocking version for first vector only.
/// Using accumulators and unrolling over dimension x.
///
void mil_block1(const float *hr_sphere_region, int n, float *directions_vectors_mil) {

    const int n_vectors = NUM_DIRECTIONS;
    float directions_vectors_bone_length[n_vectors];
    unsigned int directions_vectors_intercepts[n_vectors];
    unsigned int count = 0;
    static unsigned int flag = 0;

    for (int kk = 0; kk < n; kk+=BLOCK_SIZE) {
        for (int jj = 0; jj < n; jj+=BLOCK_SIZE) {
            for (int ii = 0; ii < n; ii+=BLOCK_SIZE) {
                /* for every direction vector */
                for (int v = 0; v < n_vectors; ++v) {

                    /* Init accumulators */
                    float acc1 = 0.0, acc5 = 0.0;
                    float acc2 = 0.0, acc6 = 0.0;
                    float acc3 = 0.0, acc7 = 0.0;
                    float acc4 = 0.0, acc8 = 0.0;

                    unsigned int edge_count1 = 1, edge_count5 = 1;
                    unsigned int edge_count2 = 1, edge_count6 = 1;
                    unsigned int edge_count3 = 1, edge_count7 = 1;
                    unsigned int edge_count4 = 1, edge_count8 = 1;

                    unsigned int prev_mask = hr_sphere_region[0] > 0.5;

                    for (int k = kk; k < kk + BLOCK_SIZE; k += STRIDE) {
                        for (int j = jj; j < jj + BLOCK_SIZE; j += STRIDE) {
                            for (int i = ii; i < ii + BLOCK_SIZE; i += NUM_ACC) {

                                unsigned int curr_mask1; //, curr_mask5;
                                unsigned int curr_mask2; //, curr_mask6;
                                unsigned int curr_mask3; //, curr_mask7;
                                unsigned int curr_mask4; //, curr_mask8;

                                /* Load working set */
                                const float* ptr_hr = &hr_sphere_region[ k*n*n + j*n + i];
                                float r1 = *ptr_hr++;  //float r5 = *ptr_hr++;
                                float r2 = *ptr_hr++;  //float r6 = *ptr_hr++;
                                float r3 = *ptr_hr++;  //float r7 = *ptr_hr++;
                                float r4 = *ptr_hr++;  //float r8 = *ptr_hr;

                                /* Perform computation */
                                acc1 += r1;  //acc5 += r5;
                                acc2 += r2;  //acc6 += r6;
                                acc3 += r3;  //acc7 += r7;
                                acc4 += r4;  //acc8 += r8;

                                /* Calculate masks */
                                curr_mask1 = r1 > 0.5;  //curr_mask5 = r5 > 0.5;
                                curr_mask2 = r2 > 0.5;  //curr_mask6 = r6 > 0.5;
                                curr_mask3 = r3 > 0.5;  //curr_mask7 = r7 > 0.5;
                                curr_mask4 = r4 > 0.5;  //curr_mask8 = r8 > 0.5;

                                /* Detect edge and add to counter */
                                edge_count1 += curr_mask1 ^ prev_mask;   //edge_count5 += curr_mask5 ^ curr_mask4;
                                edge_count2 += curr_mask2 ^ curr_mask1;  //edge_count6 += curr_mask6 ^ curr_mask5;
                                edge_count3 += curr_mask3 ^ curr_mask2;  //edge_count7 += curr_mask7 ^ curr_mask6;
                                edge_count4 += curr_mask4 ^ curr_mask3;  //edge_count8 += curr_mask8 ^ curr_mask7;

                                /* Update state of prev_mask */
                                prev_mask = curr_mask4;
                                count++;
                            }
                        }
                    }

                    acc1 += acc2;   edge_count1 += edge_count2;
                    acc3 += acc4;   edge_count3 += edge_count4;
                    acc5 += acc6;   edge_count5 += edge_count6;
                    acc7 += acc8;   edge_count7 += edge_count8;
                    acc1 += acc3;   edge_count1 += edge_count3;
                    acc5 += acc7;   edge_count5 += edge_count7;
                    directions_vectors_bone_length[v] = acc1 + acc5;
                    directions_vectors_intercepts[v]  = edge_count1 + edge_count5;
                }
            }
        }
    }

    if (!flag) {
        printf("%d\n", count);
        flag = 1;
    }

    for (int v = 0; v < n_vectors; ++v) {
        directions_vectors_mil[v] = directions_vectors_bone_length[v] / directions_vectors_intercepts[v];
    }


}