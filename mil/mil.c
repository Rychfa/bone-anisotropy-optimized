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
int find_maximum_absolute_value(int *array, int length) {
    int max = -1;
    for (int i = 0; i < length; ++i) {
        if (max < abs(array[i])) {
            max = abs(array[i]);
        }
    }
    return max;
}

double *init_mil_vector(int n_vectors) {
    double *directions_vectors_mil = (double *) malloc(n_vectors * sizeof(double));
    for (int j = 0; j < n_vectors; ++j) {
        directions_vectors_mil[j] = 0.0;
    }
    return directions_vectors_mil;
}

int **randomly_generate_central_points(int n_central_point, int max_coordinate) {
    int **central_points = (int**)malloc(n_central_point * sizeof(int*));

    for (int i = 0; i < n_central_point; ++i) {
        central_points[i] = (int*)malloc(3 * sizeof(int));


        int z = rand() % max_coordinate; //arc4random_uniform(max_coordinate);
        int y = rand() % max_coordinate; //arc4random_uniform(max_coordinate);
        int x = rand() % max_coordinate; //arc4random_uniform(max_coordinate);

        central_points[i] = (int[3]){z, y, x};
    }

    return central_points;
}

//double l1_norm(double *vector, int size) {
//    double l1_norm = 0;
//    for (int i = 0; i < size; ++i) {
//        l1_norm += fabs(vector[i]);
//    }
//    return l1_norm;
//}
//
//void validate_direction_vectors(double** direction_vectors, int n, int d) {
//    for (int i = 0; i < n; ++i) {
//        double direction_vector_l1_norm = l1_norm(direction_vectors[i], d);
//        assert(direction_vector_l1_norm == 1.0);
//    }
//}

<<<<<<< HEAD
double *mil(int ***hr_sphere_region, int n, double directions_vectors[][3], int n_vectors) {
=======
double *mil(int ***hr_sphere_region, int n, double **directions_vectors, int n_vectors, int dimension) {
>>>>>>> devel_jp
//    validate_direction_vectors(directions_vectors, n_vectors, dimension);

    double *directions_vectors_mil = init_mil_vector(n_vectors);

    double directions_vectors_bone_length[n_vectors], directions_vectors_intercepts[n_vectors];
    for (int j = 0; j < n_vectors; ++j) {
        directions_vectors_bone_length[j] = 0.0;
        directions_vectors_intercepts[j] = 0.0;
    }

    // todo: for now, statically initialized, only one central Point.
    //  Will be replaced by: randomly_generate_central_points(,). Deallocate memory at the end !
    int n_central_point = 1;
    int central_points[1][3] = {{n / 2 - 1, n / 2 - 1, n / 2 - 1}};

    for (int k = 0; k < n_central_point; ++k) {
        int *central_point = central_points[k];

        for (int i = 0; i < n_vectors; ++i) {
            double *direction_vector = directions_vectors[i];

            double z = central_point[0], y = central_point[1], x = central_point[2];
            int z_int = (int) z, y_int = (int) y, x_int = (int) x;
            int vector_state = hr_sphere_region[z_int][y_int][x_int];

            while ((z < n && z >= 0) && (y < n && y >= 0) && (x < n && x >= 0)) {
                z_int = (int) z, y_int = (int) y, x_int = (int) x;

                int current_state = hr_sphere_region[z_int][y_int][x_int];

                if (current_state == 1) {
                    ++directions_vectors_bone_length[i];
                }

                if (current_state != vector_state) {
                    ++directions_vectors_intercepts[i];
                    vector_state = current_state;
                }

                z += direction_vector[0];
                y += direction_vector[1];
                x += direction_vector[2];
            }
        }
    }

    for (int i = 0; i < n_vectors; ++i) {
        directions_vectors_mil[i] =
                directions_vectors_intercepts[i] > 0.0 ? (double) directions_vectors_bone_length[i] /
                                                         directions_vectors_intercepts[i]
                                                       : directions_vectors_bone_length[i];
    }

    return directions_vectors_mil;
}


double *mil2(int ***hr_sphere_region, int n, double **directions_vectors, int n_vectors, int dimension) {

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
            } else if (direction_vector[dim] > 0.0) {
                final_point[dim] = (int) floor(direction_vector[dim] * n) - 1;
            } else {
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
<<<<<<< HEAD
}
=======
}
>>>>>>> devel_jp
