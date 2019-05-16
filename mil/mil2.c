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
             /* { k,  j,  i} */
                { 0,  1,  1},
                { 1,  0,  1},
                { 1,  1,  0},
        };

int get_iterator_vectors (int directionVector[3], int iteratorVectors[][3]) {

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