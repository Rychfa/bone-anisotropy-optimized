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
#include <stdio.h>
#include <stdlib.h>
#include "kernel.h"
#include "reader.h"
#include "region_extraction.h"
#include "mil.h"
#include "test_mil.h"
#include "ellipsoid.h"
#include "utils.h"
#include "eigen.h"

void init (int** sphere, int** ptrHighRes, int** ptrLowRes, double** rotation_matrix, double** ptrEvecOut) {
    //
    // Create sphere mask
    //
    *sphere = malloc( (sizeof (int)) * SPHERE_ARRAY_SIZE);
    createSphereMask(*sphere);
    //writeVTK(*sphere, HIGH_RES_VOXEL_SIZE, SPHERE_NDIM, SPHERE_NDIM, SPHERE_NDIM, "test/sphere.vtk");
    //
    // Read input images
    //
    *ptrHighRes = readHighResImage();
    *ptrLowRes = readLowResImage();
    //writeVTK(*ptrHighRes, HIGH_RES_VOXEL_SIZE, HIGH_RES_D1, HIGH_RES_D2, HIGH_RES_D3, "test/debug_hr_image.vtk");
    //
    // Init rotation matrix
    //
    *rotation_matrix = malloc (sizeof(double) * 9);
    (*rotation_matrix)[0] =  0.91241511;
    (*rotation_matrix)[1] =  0.02857433;
    (*rotation_matrix)[2] = -0.40826728;
    (*rotation_matrix)[3] = -0.0259223;
    (*rotation_matrix)[4] =  0.99959159;
    (*rotation_matrix)[5] =  0.01202831;
    (*rotation_matrix)[6] =  0.40844424;
    (*rotation_matrix)[7] = -3.91584012e-04;
    (*rotation_matrix)[8] =  0.91278319;
    //
    // Allocate space for the output eigen vectors
    //
    *ptrEvecOut = calloc (sizeof(double), 3*3*LOW_RES_SIZE);
}

void deInit (int* sphere, int* ptrHighRes, int* ptrLowRes, double* rotation_matrix, double* ptrEvecOut) {
    free(sphere);
    free(ptrHighRes);
    free(ptrLowRes);
    free(rotation_matrix);
    free(ptrEvecOut);
}

void kernel_basic (int* sphere, int* ptrHighRes, int* ptrLowRes, double* rotation_matrix, double* ptrEvecOut) {
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

    double r00, r01, r02, r10, r11, r12, r20, r21, r22;
    double voxel_size_lr, voxel_size_hr;
    double half_voxel_size_lr, half_voxel_size_hr;
    double xDr00, xDr10, xDr20, yDr01, yDr11, yDr21, zDr02, zDr12, zDr22;
    double tx_hr, ty_hr, tz_hr;
    double nx_hr, ny_hr, nz_hr;
    double nx_sp, ny_sp, nz_sp;
    double x_lr, y_lr, z_lr;
    int ii_lr;
    double zero, half, one;
    int i_hr, j_hr, k_hr, ii_hr;
    int imin, jmin, kmin;
    int radius_i, radius_j, radius_k;
    int *extracted_region = malloc((sizeof (int)) * SPHERE_ARRAY_SIZE);

    FILE *fd;

    zero = 0.0;
    half = 0.5;
    one = 1.0;

    // XCT - High resolution
    // image dimension
    nx_hr = HIGH_RES_D1;
    ny_hr = HIGH_RES_D2;
    nz_hr = HIGH_RES_D3;
    // voxel_size
    voxel_size_hr = 0.082;
    half_voxel_size_hr = 0.041;

    // QCT - Low resolution
    // voxel_size
    voxel_size_lr = 3.0;
    half_voxel_size_lr = 1.5;

    // Rotation matrix
    r00 = rotation_matrix[0];
    r01 = rotation_matrix[1];
    r02 = rotation_matrix[2];
    r10 = rotation_matrix[3];
    r11 = rotation_matrix[4];
    r12 = rotation_matrix[5];
    r20 = rotation_matrix[5];
    r21 = rotation_matrix[6];
    r22 = rotation_matrix[7];

    printf("coordMap: r [%.3f, %.3f, %.3f]\n", r00, r01, r02);
    printf("coordMap: r [%.3f, %.3f, %.3f]\n", r10, r11, r12);
    printf("coordMap: r [%.3f, %.3f, %.3f]\n", r20, r21, r22);
    printf("coordMap: center of rotation: %.3f ,%.3f, %.3f\n", xC, yC, zC);
    printf("coordMap: translation: %.3f ,%.3f, %.3f\n", xT, yT, zT);

     // debug
     region_extraction(366, 341, 915, sphere, extracted_region, ptrHighRes); 
     writeVTK(extracted_region, HIGH_RES_VOXEL_SIZE, SPHERE_NDIM, SPHERE_NDIM, SPHERE_NDIM, "test/region.vtk"); 
    //

    fd = fopen("output.txt","w");

    // loop over all femur voxels
    for (int k_lr=0; k_lr < LOW_RES_D3; k_lr++)
    {
        // calculate vector from the center of the image
        // to the center of this voxel
        z_lr = k_lr*voxel_size_lr - zC + half_voxel_size_lr;
        // multiply the vector with the rotation matrix
        zDr02 = z_lr * r02;
        zDr12 = z_lr * r12;
        zDr22 = z_lr * r22;

        for (int j_lr=0; j_lr < LOW_RES_D2; j_lr++)
        {
            //calculate vector from the center of the image
            //to the center of this voxel
            y_lr = j_lr*voxel_size_lr - yC + half_voxel_size_lr;
            //multiply the vector with the rotation matrix
            yDr01 = y_lr * r01;
            yDr11 = y_lr * r11;
            yDr21 = y_lr * r21;

            for (int i_lr=0; i_lr < LOW_RES_D1; i_lr++)
            {
                // calculate vector from the center of the image
                // to the center of this voxel
                x_lr = i_lr*voxel_size_lr - xC + half_voxel_size_lr;
                // calculate index
                ii_lr = i_lr + j_lr*LOW_RES_D1 + k_lr*LOW_RES_D1*LOW_RES_D2;
                // check if this voxel inside the FE mask
                if (ptrLowRes[ii_lr] > 0)
                {
                    //
                    // This is executed approx. one third of the times.
                    //

                    // multiply the vector with the rotation matrix
                    xDr00 = x_lr * r00;
                    xDr10 = x_lr * r10;
                    xDr20 = x_lr * r20;
                    // Find location in Fabric image
                    // the point in the original domain; tx,ty,tz are float indices,
                    // supposed to point to the SWB corner of that voxel
                    // in the original image.
                    // t = (R*D + original_center - half_voxel)/voxel_size
                    tx_hr = (xDr00 + yDr01 + zDr02 + xC + xT - half_voxel_size_hr)/voxel_size_hr;
                    ty_hr = (xDr10 + yDr11 + zDr12 + yC + yT - half_voxel_size_hr)/voxel_size_hr;
                    tz_hr = (xDr20 + yDr21 + zDr22 + zC + zT - half_voxel_size_hr)/voxel_size_hr;


                    // lower bounds
                    if ((tx_hr < zero) && (tx_hr >= -half)) {
                        tx_hr = zero;
                    }
                    else if ((ty_hr < zero) && (ty_hr >= -half)) {
                        ty_hr = zero;
                    }
                    else if ((tz_hr < zero) && (tz_hr >= -half)) {
                        tz_hr = zero;
                    }
                    // upper bounds
                    if ((tx_hr <= (nx_hr-half)) && (tx_hr > (nx_hr-one))) {
                        tx_hr = nx_hr-one; }
                    else if ((ty_hr <= (ny_hr-half)) && (ty_hr > (ny_hr-one))) {
                        ty_hr = ny_hr-one; }
                    else if ((tz_hr <= (nz_hr-half)) && (tz_hr > (nz_hr-one))) {
                        tz_hr = nz_hr-one; }
                    // get the integer index of the calculated point
                    // which is the corner of that voxel in the original image
                    i_hr = (int) tx_hr;
                    j_hr = (int) ty_hr;
                    k_hr = (int) tz_hr;

                    //printf("coordMap: lr(%d ,%d, %d), hr(%d, %d, %d)\n", i_lr, j_lr, k_lr, i_hr, j_hr, k_hr);

                    // extract a sphere region
                    //region_extraction(i_hr, j_hr, k_hr, sphere, extracted_region, ptrHighRes);
                    // compute fabric
                    double mils[NUM_DIRECTIONS];
                    mil( extracted_region, SPHERE_NDIM, DIRECTIONS, NUM_DIRECTIONS, mils);
                    //print_vector(mils, NUM_DIRECTIONS);

                    double Q[3][3];
                    fit_ellipsoid_mils(mils, Q);

                    double eVecs[3][3];
                    double eVals[3];
                    // eigen3(Q, &ptrEvecOut[ii_lr*9], eVals); //TODO:.. ?
                    eigen3(Q, eVecs, eVals);

                    fprintf(fd, "%d, %d, %d, " /* index into LR image */
                        "%.3f, %.3f, %.3f, "   /* e vals */
                        "%.3f, %.3f, %.3f, "   /* e vecs*/
                        "%.3f, %.3f, %.3f, "
                        "%.3f, %.3f, %.3f\n",
                        i_lr, j_lr, k_lr,
                        eVals[0], eVals[1], eVals[2],
                        eVecs[0][0], eVecs[0][1], eVecs[0][2],
                        eVecs[1][0], eVecs[1][1], eVecs[1][2],
                        eVecs[2][0], eVecs[2][1], eVecs[2][2]);

                }

            }
        }
    }
    fclose(fd);
}
