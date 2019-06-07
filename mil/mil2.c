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

    gBone1 = directions_vectors_bone_length[10];
    gInter1 = directions_vectors_intercepts[10];
#endif

}

void mil2_scalar(const double *hr_sphere_region, int n, double *directions_vectors_mil) {
    double bone_length[13] = {0.0};
    int intercepts[13] = {0};

    for (int kb = 0; kb < n; kb+=BLOCK_SIZE) {
        for (int jb = 0; jb < n; jb+=BLOCK_SIZE) {
            for (int ib = 0; ib < n; ib+=BLOCK_SIZE) {

                BLOCK_KERNEL_1D(1, kb, jb, ib)
                BLOCK_KERNEL_1D(2, kb, ib, jb)
                BLOCK_KERNEL_1D(3, jb, ib, kb)
                BLOCK_KERNEL_2D(4, kb, jb, ib)
                BLOCK_KERNEL_2D(5, jb, kb, ib)
                BLOCK_KERNEL_2D(6, ib, jb, kb)
                BLOCK_KERNEL_2D_NEG(7, kb, jb, ib)
                BLOCK_KERNEL_2D_NEG(8, jb, kb, ib)
                BLOCK_KERNEL_2D_NEG(9, ib, jb, kb)
                BLOCK_KERNEL_3D(10, kb, jb, ib)
                BLOCK_KERNEL_3D(11, kb, jb, ib)
                BLOCK_KERNEL_3D(12, kb, jb, ib)
                BLOCK_KERNEL_3D(13, kb, jb, ib)
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