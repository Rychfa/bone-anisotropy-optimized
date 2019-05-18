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

#define STRIDE 2

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

void mil2_baseline(const double *hr_sphere_region, int n, double *directions_vectors_mil) {

    const int n_vectors = NUM_DIRECTIONS;

    double directions_vectors_bone_length[n_vectors];
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
                        double current_point;
                        double prev_point = hr_sphere_region[ k*n*n + j*n + i];
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
#define NUM_ACC 8
void dummy0(const float *hr_sphere_region, int n, float *directions_vectors_mil) {

    const int n_vectors = NUM_DIRECTIONS;
    int count = 0;

#if 1
    for (int kk = 0; kk < n; kk+=BLOCK_SIZE) {
        for (int jj = 0; jj < n; jj+=BLOCK_SIZE) {
            for (int ii = 0; ii < n; ii+=BLOCK_SIZE) {
                /* for every direction vector */
                for (int v = 0; v < 1; ++v) {
                    float acc1 = 0.0;  float acc5 = 0.0;
                    float acc2 = 0.0;  float acc6 = 0.0;
                    float acc3 = 0.0;  float acc7 = 0.0;
                    float acc4 = 0.0;  float acc8 = 0.0;

                    for (int k = kk; k < kk + BLOCK_SIZE; k += STRIDE) {
                        for (int j = jj; j < jj + BLOCK_SIZE; j += STRIDE) {
                            for (int i = ii; i < ii + BLOCK_SIZE; i += NUM_ACC) {
                                acc1 += hr_sphere_region[ k*n*n + j*n + i];
                                acc2 += hr_sphere_region[ k*n*n + j*n + i + 1];
                                acc3 += hr_sphere_region[ k*n*n + j*n + i + 2];
                                acc4 += hr_sphere_region[ k*n*n + j*n + i + 3];
                                acc5 += hr_sphere_region[ k*n*n + j*n + i + 4];
                                acc6 += hr_sphere_region[ k*n*n + j*n + i + 5];
                                acc7 += hr_sphere_region[ k*n*n + j*n + i + 6];
                                acc8 += hr_sphere_region[ k*n*n + j*n + i + 7];
                            }
                        }
                    }

                    acc1 += acc2;
                    acc3 += acc4;
                    acc5 += acc6;
                    acc7 += acc8;
                    acc1 += acc3;
                    acc5 += acc7;
                    directions_vectors_mil[v] = acc1 + acc5;
                }
            }
        }
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

///
/// Blocking version for first vector only.
/// Using accumulators and unrolling over dimension y.
///
void mil_block2(const float *hr_sphere_region, int n, float *directions_vectors_mil) {

    const int n_vectors = NUM_DIRECTIONS;
    float directions_vectors_bone_length[n_vectors];
    unsigned int directions_vectors_intercepts[n_vectors];

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

                            //unsigned int i_prev = (ii > 0) ? ii - 1 : 0;

                            unsigned int prev_mask1 = 0; //hr_sphere_region[ k*n*n + (j+0*STRIDE)*n + i_prev ]  > 0.5;
                            unsigned int prev_mask2 = 0; //hr_sphere_region[ k*n*n + (j+1*STRIDE)*n + i_prev ]  > 0.5;
                            unsigned int prev_mask3 = 0; //hr_sphere_region[ k*n*n + (j+2*STRIDE)*n + i_prev ]  > 0.5;
                            unsigned int prev_mask4 = 0; //hr_sphere_region[ k*n*n + (j+3*STRIDE)*n + i_prev ]  > 0.5;

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
}


///
/// Blocking version for first vector only.
/// Using accumulators and unrolling over dimension x.
///
void mil_block1(const float *hr_sphere_region, int n, float *directions_vectors_mil) {

    const int n_vectors = NUM_DIRECTIONS;
    float directions_vectors_bone_length[n_vectors];
    unsigned int directions_vectors_intercepts[n_vectors];
    int count = 0;

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


}


void dummy1(const double *hr_sphere_region, int n, double *directions_vectors_mil) {

    const int n_vectors = NUM_DIRECTIONS;
    int count = 0;

    /* for every direction vector */
    for (int v = 0; v < n_vectors; ++v) {
        double acc1 = 0;   double acc2 = 0;
        double acc3 = 0;   double acc4 = 0;
        double acc5 = 0;   double acc6 = 0;
        double acc7 = 0;   double acc8 = 0;
        for (int kk = 0 ; kk < n; kk += STRIDE) {
            for (int jj = 0; jj < n; jj += STRIDE*8) {
                int i1 = 0;   //int i2 = 0;
//                int i3 = 0;   int i4 = 0;
//                int i5 = 0;   int i6 = 0;
//                int i7 = 0;   int i8 = 0;
                while ( i1 < n && i1 >= 0) {
                    acc1 += hr_sphere_region[kk*n*n + jj*n + i1];
                    acc2 += hr_sphere_region[kk*n*n + (jj+1*STRIDE)*n + i1];
                    acc3 += hr_sphere_region[kk*n*n + (jj+2*STRIDE)*n + i1];
                    acc4 += hr_sphere_region[kk*n*n + (jj+3*STRIDE)*n + i1];
                    acc5 += hr_sphere_region[kk*n*n + (jj+4*STRIDE)*n + i1];
                    acc6 += hr_sphere_region[kk*n*n + (jj+5*STRIDE)*n + i1];
                    acc7 += hr_sphere_region[kk*n*n + (jj+6*STRIDE)*n + i1];
                    acc8 += hr_sphere_region[kk*n*n + (jj+7*STRIDE)*n + i1];

                    i1 += 1;   //i2 += 1;
//                    i3 += 1;   i4 += 1;
//                    i5 += 1;   i6 += 1;
//                    i7 += 1;   i8 += 1;
//                    count+=8;
                }
            }
        }

        directions_vectors_mil[v] = acc1 + acc2 + acc3 + acc4 + acc5 + acc6 + acc7 + acc8;
    }

//    printf("count: %d\n", count);

}

void dummy2(const int *hr_sphere_region, int n, double *directions_vectors_mil) {

    const int n_vectors = NUM_DIRECTIONS;

    /* for every direction vector */
    for (int v = 0; v < n_vectors; ++v) {
        unsigned int current_point;
        double sum = 0.0;
        for (int kk = 0 ; kk < n; kk += STRIDE) {
            for (int jj = 0; jj < n; jj += STRIDE) {
                int i = 0;
                int j = jj;
                while ( j < n && j >= 0 ) {
                    current_point = hr_sphere_region[kk*n*n + j*n + i];
                    sum += (double) current_point;
                    i += 1;
                    j += 1;
                }
            }
        }

        for (int kk = 0 ; kk < n; kk += STRIDE) {
            for (int ii = 1; ii < n; ii += STRIDE) {
                int j = 0;
                int i = ii;
                while ( i < n && i >= 0 ) {
                    current_point = hr_sphere_region[kk*n*n + j*n + i];
                    sum += (double) current_point;
                    i += 1;
                    j += 1;
                }
            }
        }
        directions_vectors_mil[v] = sum;
    }

}

void dummy3(const int *hr_sphere_region, int n, double *directions_vectors_mil) {

    const int n_vectors = NUM_DIRECTIONS;

    /* for every direction vector */
    for (int v = 0; v < n_vectors; ++v) {
        unsigned int current_point;
        double sum = 0.0;
        for (int kk = 0 ; kk < n; kk += STRIDE) {
            for (int jj = 0; jj < n; jj += STRIDE) {
                int i = 0;
                int j = jj;
                int k = kk;
                while ( (j < n && j >= 0) && (k < n && k >= 0) ) {
                    current_point = hr_sphere_region[k*n*n + j*n + i];
                    sum += (double) current_point;
                    i += 1;
                    j += 1;
                    k += 1;
                }
            }
        }

        for (int kk = 0 ; kk < n; kk += STRIDE) {
            for (int ii = 1; ii < n; ii += STRIDE) {
                int i = ii;
                int j = 0;
                int k = kk;
                while ( (i < n && i >= 0) && (k < n && k >= 0) ) {
                    current_point = hr_sphere_region[k*n*n + j*n + i];
                    sum += (double) current_point;
                    i += 1;
                    j += 1;
                    k += 1;
                }
            }
        }

        for (int jj = 1 ; jj < n; jj += STRIDE) {
            for (int ii = 1; ii < n; ii += STRIDE) {
                int i = ii;
                int j = jj;
                int k = 0;
                while ( (i < n && i >= 0) && (j < n && j >= 0) ) {
                    current_point = hr_sphere_region[k*n*n + j*n + i];
                    sum += (double) current_point;
                    i += 1;
                    j += 1;
                    k += 1;
                }
            }
        }

        directions_vectors_mil[v] = sum;
    }

}

///
///
///
void mil2_o2(const int *hr_sphere_region, int n, double *directions_vectors_mil) {

    const int n_vectors = NUM_DIRECTIONS;
    int iteratorVectors[3][3];

    /* for every direction vector */
    for (int v = 0; v < n_vectors; ++v) {

        int numIterVecs = 0;
        int d0 = DIRECTIONS[v][0];
        int d1 = DIRECTIONS[v][1];
        int d2 = DIRECTIONS[v][2];
        double sum = 0;
        double intersects = 1;

        /* get array of iterator vectors */
        numIterVecs = get_iterator_vectors(DIRECTIONS[v], iteratorVectors);

        /* Iterate over all initial faces */
        for (int f = 0; f < numIterVecs;  ++f) {

            /* Iterator vector 1 */
            if (iteratorVectors[f][0] == 0 && iteratorVectors[f][1] == 1 && iteratorVectors[f][2] == 1) {
                /* Iterate through initial points in the face (0, 1, 1) */
                int ii = (d0 == -1) ? n - 1 : 0;
                for (int kk = 0 ; kk < n ; kk += STRIDE) {
                    for (int jj = 0 ; jj < n ; jj += STRIDE) {

                        int k = kk;
                        int j = jj;
                        int i = ii;
                        unsigned int current_point;
                        unsigned int prev_point = hr_sphere_region[ k*n*n + j*n + i];

                        while ((k < n && k >= 0) && (j < n && j >= 0) && (i < n && i >= 0)) {

                            current_point = hr_sphere_region[ k*n*n + j*n + i];

                            sum += (double) current_point;
                            intersects += (double) (current_point ^ prev_point);
                            prev_point = current_point;

                            k += d2;
                            j += d1;
                            i += d0;
                        }
                    }
                }
            }

            /* Iterator vector 2 */
            if (iteratorVectors[f][0] == 1 && iteratorVectors[f][1] == 0 && iteratorVectors[f][2] == 1) {
                /* Iterate through initial points in the face (0, 1, 1) */
                int jj = (d1 == -1) ? n - 1 : 0;
                for (int kk = 0 ; kk < n ; kk += STRIDE) {
                    for (int ii = 0 ; ii < n ; ii += STRIDE) {

                        int k = kk;
                        int j = jj;
                        int i = ii;
                        unsigned int current_point;
                        unsigned int prev_point = hr_sphere_region[ k*n*n + j*n + i];

                        while ((k < n && k >= 0) && (j < n && j >= 0) && (i < n && i >= 0)) {

                            current_point = hr_sphere_region[ k*n*n + j*n + i];

                            sum += (double) current_point;
                            intersects  += (double) (current_point ^ prev_point);
                            prev_point = current_point;

                            k += d2;
                            j += d1;
                            i += d0;
                        }
                    }
                }
            }

            /* Iterator vector 3 */
            if (iteratorVectors[f][0] == 1 && iteratorVectors[f][1] == 1 && iteratorVectors[f][2] == 0) {
                /* Iterate through initial points in the face (0, 1, 1) */
                int kk = (d2 == -1) ? n - 1 : 0;
                for (int jj = 0 ; jj < n ; jj += STRIDE) {
                    for (int ii = 0 ; ii < n ; ii += STRIDE) {

                        int k = kk;
                        int j = jj;
                        int i = ii;
                        unsigned int current_point;
                        unsigned int prev_point = hr_sphere_region[ k*n*n + j*n + i];

                        while ((k < n && k >= 0) && (j < n && j >= 0) && (i < n && i >= 0)) {

                            current_point = hr_sphere_region[ k*n*n + j*n + i];

                            sum += (double) current_point;
                            intersects  += (double) (current_point ^ prev_point);
                            prev_point = current_point;

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