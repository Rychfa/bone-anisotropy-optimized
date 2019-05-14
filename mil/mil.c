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
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <stdio.h>
#include "mil.h"

///
/// Skeleton for MIL
///

void randomly_generate_central_points(int *central_points, int n_central_point, int max_coordinate) {
    for (int i = 0; i < n_central_point*3; i += 3) {
        int z = rand() % max_coordinate; //arc4random_uniform(max_coordinate);
        int y = rand() % max_coordinate; //arc4random_uniform(max_coordinate);
        int x = rand() % max_coordinate; //arc4random_uniform(max_coordinate);

        central_points[i] = z;
        central_points[i + 1] = y;
        central_points[i + 2] = x;
    }
}

int mil(int *hr_sphere_region, int n, double directions_vectors[][3], int n_vectors,
        double *directions_vectors_mil) {

    int flops_counter = 0;
    static int central_points[N_CENTRAL_POINTS*3];

    randomly_generate_central_points(central_points, N_CENTRAL_POINTS, n);

    double directions_vectors_bone_length[n_vectors], directions_vectors_intercepts[n_vectors];

    for (int j = 0; j < n_vectors; ++j) {
        directions_vectors_bone_length[j] = 0.0;
        directions_vectors_intercepts[j] = N_CENTRAL_POINTS;
    }

    for (int k = 0; k < N_CENTRAL_POINTS*3; k += 3) {
        int central_point[3] = {central_points[k], central_points[k + 1], central_points[k + 2]};

        for (int i = 0; i < n_vectors; ++i) {
            double *direction_vector = directions_vectors[i];

            // printf("%f %f %f\n", direction_vector[0], direction_vector[1], direction_vector[2]);

            double z = central_point[0], y = central_point[1], x = central_point[2];
            int z_int = (int) z, y_int = (int) y, x_int = (int) x;
            int vector_state = hr_sphere_region[z_int * n * n + y_int * n + x_int];

            while ((z < n && z >= 0) && (y < n && y >= 0) && (x < n && x >= 0)) {
                z_int = (int) z, y_int = (int) y, x_int = (int) x;

                int current_state = hr_sphere_region[z_int * n * n + y_int * n + x_int];

                if (current_state == 1) {
#if MIL_FLOPS_COUNT > 0
                    ++flops_counter;
#endif
                    ++directions_vectors_bone_length[i];
                }

                if (current_state != vector_state) {
#if MIL_FLOPS_COUNT > 0
                    ++flops_counter;
#endif
                    ++directions_vectors_intercepts[i];
                    vector_state = current_state;
                }
#if MIL_FLOPS_COUNT > 0
                flops_counter += 3;
#endif
                z += direction_vector[2];
                y += direction_vector[1];
                x += direction_vector[0];
            }
        }
    }

    for (int i = 0; i < n_vectors; ++i) {
#if MIL_FLOPS_COUNT > 0
        ++flops_counter;
#endif
        directions_vectors_mil[i] = directions_vectors_bone_length[i] / directions_vectors_intercepts[i];
    }

    return flops_counter;
}