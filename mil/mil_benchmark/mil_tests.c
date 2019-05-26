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

///
/// Test for first vector only (1,0,0).
///
void mil_test_v1(const double *hr_sphere_region, int n, double *directions_vectors_mil) {

    double bone_length[13] = {0.0};
    int intercepts[13] = {0};
    for (int kk = 0; kk < n; kk+=BLOCK_SIZE) {
        for (int jj = 0; jj < n; jj+=BLOCK_SIZE) {
            for (int ii = 0; ii < n; ii+=BLOCK_SIZE) {

                for (int v = 0; v < 13; ++v) {
                    BLOCK_KERNEL_1D(1,kk,kk,ii)
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

///
/// Test for first vector only (0,1,0).
///
void mil_test_v2(const double *hr_sphere_region, int n, double *directions_vectors_mil) {

    double bone_length[13] = {0.0};
    int intercepts[13] = {0};
    for (int kk = 0; kk < n; kk+=BLOCK_SIZE) {
        for (int jj = 0; jj < n; jj+=BLOCK_SIZE) {
            for (int ii = 0; ii < n; ii+=BLOCK_SIZE) {

                for (int v = 0; v < 13; ++v) {
                    BLOCK_KERNEL_1D(2,kk,kk,ii)
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

///
/// Test for first vector only (0,0,1).
///
void mil_test_v3(const double *hr_sphere_region, int n, double *directions_vectors_mil) {

    double bone_length[13] = {0.0};
    int intercepts[13] = {0};
    for (int kk = 0; kk < n; kk+=BLOCK_SIZE) {
        for (int jj = 0; jj < n; jj+=BLOCK_SIZE) {
            for (int ii = 0; ii < n; ii+=BLOCK_SIZE) {

                for (int v = 0; v < 13; ++v) {
                    BLOCK_KERNEL_1D(3,kk,kk,ii)
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

///
/// Test for first vector only (1,1,0).
///
void mil_test_v4(const double *hr_sphere_region, int n, double *directions_vectors_mil) {

    double bone_length[13] = {0.0};
    int intercepts[13] = {0};
    for (int kk = 0; kk < n; kk+=BLOCK_SIZE) {
        for (int jj = 0; jj < n; jj+=BLOCK_SIZE) {
            for (int ii = 0; ii < n; ii+=BLOCK_SIZE) {

                for (int v = 0; v < 13; ++v) {
                    BLOCK_KERNEL_2D(4,kk,kk,ii)
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

///
/// Test for first vector only (1,0,1).
///
void mil_test_v5(const double *hr_sphere_region, int n, double *directions_vectors_mil) {

    double bone_length[13] = {0.0};
    int intercepts[13] = {0};
    for (int kk = 0; kk < n; kk+=BLOCK_SIZE) {
        for (int jj = 0; jj < n; jj+=BLOCK_SIZE) {
            for (int ii = 0; ii < n; ii+=BLOCK_SIZE) {

                for (int v = 0; v < 13; ++v) {
                    BLOCK_KERNEL_2D(5,kk,kk,ii)
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

///
/// Test for first vector only (0,1,1).
///
void mil_test_v6(const double *hr_sphere_region, int n, double *directions_vectors_mil) {

    double bone_length[13] = {0.0};
    int intercepts[13] = {0};
    for (int kk = 0; kk < n; kk+=BLOCK_SIZE) {
        for (int jj = 0; jj < n; jj+=BLOCK_SIZE) {
            for (int ii = 0; ii < n; ii+=BLOCK_SIZE) {

                for (int v = 0; v < 13; ++v) {
                    BLOCK_KERNEL_2D(6,kk,kk,ii)
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