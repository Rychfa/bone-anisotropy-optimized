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
#include "stdio.h"
#include "reader.h"
#include "coord_map.h"
#include "region_extraction.h"
#include "bvtv.h"
#include "mil.h"
#include "ellipsoid.h"
#include "utils.h"

int main () {
    //
    // Call everything
    //
//    reader();
//    coordMap();
//    region_extraction();
//    mil();
//    bvtv();
    int *sphere = malloc( (sizeof (int*)) * SPHERE_ARRAY_SIZE);


    /* calculateRotationMatrix(rotation_matrix, ax, ay, az); */

    /* printf("rotation_matrix: %.3f, %.3f, %.3f \n", rotation_matrix[0], rotation_matrix[1], rotation_matrix[2]); */

    createSphereMask(&sphere);
    //writeVTK(&sphere, HIGH_RES_VOXEL_SIZE, SPHERE_NDIM);


    double* ptrHighRes; // = readHighResImage();
    double* ptrLowRes = readLowResImage();

    double rotation_matrix[9];
    double ax, ay, az, xT, yT, zT, xC, yC, zC;

    ax = -0.000429;
    ay = -0.420749;
    az = -0.028403;
    xT = 4.978221;
    yT = 9.715922;
    zT = 19.025912;
    xC = 40.7950011268;
    yC = 22.7550006285;
    zC = 47.3550013080;

    rotation_matrix[0] = 0.91241511;
    rotation_matrix[1] = 0.02857433;
    rotation_matrix[2] = -0.40826728;
    rotation_matrix[3] = -0.0259223;
    rotation_matrix[4] = 0.99959159;
    rotation_matrix[5] = 0.01202831;
    rotation_matrix[6] = 4.08444242e-01;
    rotation_matrix[7] = -3.91584012e-04;
    rotation_matrix[8] = 9.12783188e-01;
   
    //coordMap(ptrLowRes, ptrHighRes, &sphere, rotation_matrix, xT, yT, zT, xC, yC, zC);
    
    free(sphere);
    destroyImageMatrices();

    printf("\nDone\n");

    return 0;
}
