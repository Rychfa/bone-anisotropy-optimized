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
        /*  { i,  j,  k} */
            { 0,  1,  1},
            { 1,  0,  1},
            { 1,  1,  0},
        };

typedef struct {
    int k;
    int j;
    int i;
} s_position;

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

int get_start (int ii, int n, int dir, int iter) {
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

int get_start_3D (int ii, int n, int dir, int iter) {
    int i = ii;
    if (dir == -1) {
        if (iter == 0) {
            i = n - 1;
        } else {
            i = ii + 1;
        }
    }
    return i;
}

s_position get_start_positions (s_position pos, int n, const int dir[3], const int iter[3], int num_iter_vecs) {

    s_position start;

    /* 1D and 2D vectors */
    if (num_iter_vecs <= 2) {
        start.k = get_start(pos.k, n, dir[2], iter[2]);
        start.j = get_start(pos.j, n, dir[1], iter[1]);
        start.i = get_start(pos.i, n, dir[0], iter[0]);
    }
    else {
        /* 3D vectors */
        /* (1, 1, 1) */
        if (dir[0] == 1 && dir[1] == 1 && dir[2] == 1) {
            if (iter[0] == 0) {
                start.i = pos.i;
                start.j = pos.j;
                start.k = pos.k + 2;
            }
            else if (iter[1] == 0) {
                start.i = pos.i + 2;
                start.j = pos.j;
                start.k = pos.k + 2;
            }
            else if (iter[2] == 0) {
                start.i = pos.i;
                start.j = pos.j;
                start.k = pos.k;
            }
        }

        /* (-1, 1, 1) */
        if (dir[0] == -1 && dir[1] == 1 && dir[2] == 1) {
            if (iter[0] == 0) {
                start.i = n-1;
                start.j = pos.j + 2;
                start.k = pos.k + 2;
            }
            else if (iter[1] == 0) {
                start.i = pos.i + 1;
                start.j = pos.j;
                start.k = pos.k;
            }
            else if (iter[2] == 0) {
                start.i = pos.i + 1;
                start.j = pos.j + 2;
                start.k = pos.k;
            }
        }

        /* (1, -1, 1) */
        if (dir[0] == 1 && dir[1] == -1 && dir[2] == 1) {
            if (iter[0] == 0) {
                start.i = pos.i;
                start.j = pos.j + 1;
                start.k = pos.k;
            }
            else if (iter[1] == 0) {
                start.i = pos.i + 2;
                start.j = n-1;
                start.k = pos.k + 2;
            }
            else if (iter[2] == 0) {
                start.i = pos.i + 2;
                start.j = pos.j + 1;
                start.k = pos.k;
            }
        }

        /* (1, 1, -1) */
        if (dir[0] == 1 && dir[1] == 1 && dir[2] == -1) {
            if (iter[0] == 0) {
                start.i = pos.i;
                start.j = pos.j;
                start.k = pos.k + 1;
            }
            else if (iter[1] == 0) {
                start.i = pos.i + 2;
                start.j = pos.j;
                start.k = pos.k + 1;
            }
            else if (iter[2] == 0) {
                start.i = pos.i + 2;
                start.j = pos.j + 2;
                start.k = n-1;
            }
        }
    }

    return start;
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
                        s_position face_position = {kk, jj, ii};
                        s_position start = get_start_positions(face_position, n, DIRECTIONS[v], iteratorVectors[f], numIterVecs);
                        int k = start.k;
                        int j = start.j;
                        int i = start.i;

//                        int k = get_start_position(kk, n, DIRECTIONS[v][2], iteratorVectors[f][2]);
//                        int j = get_start_position(jj, n, DIRECTIONS[v][1], iteratorVectors[f][1]);
//                        int i = get_start_position(ii, n, DIRECTIONS[v][0], iteratorVectors[f][0]);

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
    gBone1 = directions_vectors_bone_length[10];
    gInter1 = directions_vectors_intercepts[10];

}

double mil2_3D(const double *hr_sphere_region, int n, const int kb, const int jb, const int ib, const int vec) {

    {
        const int vecID = vec;
        double bone_length_block = 0.0;
        int intercepts_block = 0;

        /* Init accumulators */
        double acc1 = 0.0, acc5 = 0.0, acc9  = 0.0;
        double acc2 = 0.0, acc6 = 0.0, acc10 = 0.0;
        double acc3 = 0.0, acc7 = 0.0, acc11 = 0.0;
        double acc4 = 0.0, acc8 = 0.0, acc12 = 0.0;

        unsigned int edge_count1 = 0, edge_count5 = 0, edge_count9  = 0;
        unsigned int edge_count2 = 0, edge_count6 = 0, edge_count10 = 0;
        unsigned int edge_count3 = 0, edge_count7 = 0, edge_count11 = 0;
        unsigned int edge_count4 = 0, edge_count8 = 0, edge_count12 = 0;

        unsigned int curr_mask1, curr_mask2,  curr_mask3,  curr_mask4;
        unsigned int curr_mask5, curr_mask6,  curr_mask7,  curr_mask8;
        unsigned int curr_mask9, curr_mask10, curr_mask11, curr_mask12;

        unsigned int prev_mask1, prev_mask4, prev_mask7, prev_mask10;
        unsigned int prev_mask2, prev_mask5, prev_mask8, prev_mask11;
        unsigned int prev_mask3, prev_mask6, prev_mask9, prev_mask12;

        int prev1, prev3, prev5;
        int prev2, prev4, prev6;

        double r1, r2,  r3,  r4;
        double r5, r6,  r7,  r8;
        double r9, r10, r11, r12;

        for (int ks = 2; ks < BLOCK_SIZE; ks += 2*STRIDE) {
            for (int js = ks + 2; js < BLOCK_SIZE; js += STRIDE) {
                int k1 = ks;
                int j1 = js;
                int k2 = js;
                int j2 = ks - 2;
                int i = 0;

                LOAD_PREV_3D(vecID)

                while (j1 < BLOCK_SIZE) {

                    LOAD_DATA_3D(vecID)

                    acc1 += r1; acc5 += r5; acc9  += r9 ;
                    acc2 += r2; acc6 += r6; acc10 += r10;
                    acc3 += r3; acc7 += r7; acc11 += r11;
                    acc4 += r4; acc8 += r8; acc12 += r12;

                    /* Calculate masks */
                    curr_mask1 = r1 > 0.5; curr_mask5 = r5 > 0.5; curr_mask9  = r9  > 0.5;
                    curr_mask2 = r2 > 0.5; curr_mask6 = r6 > 0.5; curr_mask10 = r10 > 0.5;
                    curr_mask3 = r3 > 0.5; curr_mask7 = r7 > 0.5; curr_mask11 = r11 > 0.5;
                    curr_mask4 = r4 > 0.5; curr_mask8 = r8 > 0.5; curr_mask12 = r12 > 0.5;

                    /* Detect edge and add to counter */
                    edge_count1 += curr_mask1 ^ prev_mask1;
                    edge_count2 += curr_mask2 ^ prev_mask2;
                    edge_count3 += curr_mask3 ^ prev_mask3;
                    edge_count4 += curr_mask4 ^ prev_mask4;

                    edge_count5 += curr_mask5 ^ prev_mask5;
                    edge_count6 += curr_mask6 ^ prev_mask6;
                    edge_count7 += curr_mask7 ^ prev_mask7;
                    edge_count8 += curr_mask8 ^ prev_mask8;

                    edge_count9  += curr_mask9  ^ prev_mask9 ;
                    edge_count10 += curr_mask10 ^ prev_mask10;
                    edge_count11 += curr_mask11 ^ prev_mask11;
                    edge_count12 += curr_mask12 ^ prev_mask12;

                    prev_mask1 = curr_mask1;  prev_mask5 = curr_mask5;  prev_mask9  = curr_mask9 ;
                    prev_mask2 = curr_mask2;  prev_mask6 = curr_mask6;  prev_mask10 = curr_mask10;
                    prev_mask3 = curr_mask3;  prev_mask7 = curr_mask7;  prev_mask11 = curr_mask11;
                    prev_mask4 = curr_mask4;  prev_mask8 = curr_mask8;  prev_mask12 = curr_mask12;

                    ++k1;
                    ++j1;
                    ++k2;
                    ++j2;
                    ++i;
                }
            }
        }

        /* Calculate remainder */
        for (int jk_s = 2; jk_s < BLOCK_SIZE; jk_s += 2*STRIDE) {
            int k = jk_s;
            int j = jk_s;
            int i = 0;

            /* Start remainder complete vectors */
            LOAD_PREV_REMAINDER_3D(vecID)

            while (j < BLOCK_SIZE) {
                LOAD_DATA_REMAINDER_3D(vecID)

                acc1 += r1; acc5 += r5; acc9  += r9 ;
                acc2 += r2; acc6 += r6; acc10 += r10;

                /* Calculate masks */
                curr_mask1 = r1 > 0.5; curr_mask5 = r5 > 0.5; curr_mask9  = r9  > 0.5;
                curr_mask2 = r2 > 0.5; curr_mask6 = r6 > 0.5; curr_mask10 = r10 > 0.5;

                /* Detect edge and add to counter */
                edge_count1 += curr_mask1 ^ prev_mask1;
                edge_count2 += curr_mask2 ^ prev_mask2;

                edge_count5 += curr_mask5 ^ prev_mask5;
                edge_count6 += curr_mask6 ^ prev_mask6;

                edge_count9  += curr_mask9  ^ prev_mask9 ;
                edge_count10 += curr_mask10 ^ prev_mask10;

                prev_mask1 = curr_mask1;  prev_mask5 = curr_mask5;  prev_mask9  = curr_mask9 ;
                prev_mask2 = curr_mask2;  prev_mask6 = curr_mask6;  prev_mask10 = curr_mask10;

                ++k;
                ++j;
                ++i;
            }
        }

        /* Calculate diagonal */
        DIAGONAL_3D(vecID)

        /* Sum up accumulators */
        acc1 += acc2;
        acc3 += acc4;
        acc5 += acc6;
        acc7 += acc8;
        acc9 += acc10;
        acc11 += acc12;
        edge_count1  += edge_count2;
        edge_count3  += edge_count4;
        edge_count5  += edge_count6;
        edge_count7  += edge_count8;
        edge_count9  += edge_count10;
        edge_count11 += edge_count12;

        bone_length_block += acc1 + acc3 + acc5 + acc7 + acc9 + acc11;
        intercepts_block  += edge_count1 + edge_count3 + edge_count5 + edge_count7 + edge_count9 + edge_count11;

        //bone_length[vecID-1] += bone_length_block;
        //intercepts[vecID-1]  += intercepts_block;

        return bone_length_block / intercepts_block;
    }


}

///
/// Test all vectors.
///
void mil_test_all(const double *hr_sphere_region, int n, double *directions_vectors_mil) {

    double bone_length[13] = {0.0};
    int intercepts[13] = {0};

    for (int kk_b = 0; kk_b < n; kk_b+=BLOCK_SIZE) {
        for (int jj_b = 0; jj_b < n; jj_b+=BLOCK_SIZE) {
            for (int ii_b = 0; ii_b < n; ii_b+=BLOCK_SIZE) {

                BLOCK_KERNEL_1D(1, kk_b, jj_b, ii_b)
                BLOCK_KERNEL_1D(2, kk_b, ii_b, jj_b)
                BLOCK_KERNEL_1D(3, jj_b, ii_b, kk_b)
                BLOCK_KERNEL_2D(4, kk_b, jj_b, ii_b)
                BLOCK_KERNEL_2D(5, jj_b, kk_b, ii_b)
                BLOCK_KERNEL_2D(6, ii_b, jj_b, kk_b)
                BLOCK_KERNEL_2D_NEG(7, kk_b, jj_b, ii_b)
                BLOCK_KERNEL_2D_NEG(8, jj_b, kk_b, ii_b)
                BLOCK_KERNEL_2D_NEG(9, ii_b, jj_b, kk_b)
                BLOCK_KERNEL_3D(10, kk_b, jj_b, ii_b)
//                bone_length[9]  = mil2_3D(hr_sphere_region, n, kk_b, jj_b, ii_b, 10);
//                bone_length[10] = mil2_3D(hr_sphere_region, n, kk_b, jj_b, ii_b, 11);
//                bone_length[11] = mil2_3D(hr_sphere_region, n, kk_b, jj_b, ii_b, 12);
//                bone_length[12] = mil2_3D(hr_sphere_region, n, kk_b, jj_b, ii_b, 13);
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