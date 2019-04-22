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
#include "mil.h"

///
/// Skeleton for MIL
///
int find_maximum_absolute_value(int *array, int length) {
    int max = -1;
    for (int i = 0; i < length; ++i) {
        if (max < abs(array[i])) {
            max = abs(array[i]);
        }
    }
    return max;
}

double* init_mil_vector(int n_vectors) {
    double *directions_vectors_mil = (double *) malloc(n_vectors * sizeof(double));
    for (int j = 0; j < n_vectors; ++j) {
        directions_vectors_mil[j] = 0.0;
    }
    return directions_vectors_mil;
}

double *mil(int ***hr_sphere_region, int n, double **directions_vectors, int n_vectors, int dimension) {

    //todo: Extend to n point: loop over n 'central point', and calculate the average MIL over the direction_vectors.
    //todo: Now, for simplicity, only one point.
    int central_point[3] = {n / 2 - 1, n / 2 - 1, n / 2 - 1};

    double directions_vectors_bone_length[n_vectors], directions_vectors_intercepts[n_vectors];
    double *directions_vectors_mil = init_mil_vector(n_vectors);

    // for each direction vector
    for (int i = 0; i < n_vectors; ++i) {
        double *direction_vector = directions_vectors[i]; // of dimension 3
        int final_point[dimension], coordinate_director[dimension];

        for (int dim = 0; dim < dimension; ++dim) {
            if (direction_vector[dim] == 0.0) {
                final_point[dim] = central_point[dim];
            }
            else if (direction_vector[dim] > 0.0) {
                final_point[dim] = (int) floor(direction_vector[dim] * n) - 1;
            }
            else {
                final_point[dim] = (int) floor(direction_vector[dim] * central_point[dim]) + central_point[dim];
            }
        }
        for (int dim = 0; dim < dimension; ++dim)
            coordinate_director[dim] = final_point[dim] - central_point[dim];

        // go through the direction vector, from the central point to final point.
        // get the corresponding index in the grid, and count the nbr of bone cell and intercepts
        int max_step = find_maximum_absolute_value(coordinate_director, dimension);
        double vector_step = 1.0 / max_step;
        int current_point[dimension];
        int nbr_of_bone_cell = 0;
        int nbr_of_intercept = 0;
        int vector_state = hr_sphere_region[central_point[0]][central_point[1]][central_point[2]];
        for (int step = 0; step <= max_step; ++step) {
            double t = step * vector_step;
            for (int dim = 0; dim < dimension; ++dim) {
                current_point[dim] = (int) round(central_point[dim] + t * coordinate_director[dim]);
            }

            int current_state = hr_sphere_region[current_point[0]][current_point[1]][current_point[2]];
            if (current_state == 1) {
                ++nbr_of_bone_cell;
            }
            if (current_state != vector_state) {
                ++nbr_of_intercept;
                vector_state = current_state;
            }
        }

        directions_vectors_bone_length[i] += nbr_of_bone_cell;
        directions_vectors_intercepts[i] += nbr_of_intercept;
    }

    for (int i = 0; i < n_vectors; ++i) {
        directions_vectors_mil[i] = (double) directions_vectors_bone_length[i] / directions_vectors_intercepts[i];
    }

    return directions_vectors_mil;
}