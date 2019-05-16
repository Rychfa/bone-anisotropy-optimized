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