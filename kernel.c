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
#include "mil2.h"
#include "test_mil.h"
#include "ellipsoid.h"
#include "utils.h"
#include "eigen.h"

static double* ptrHighResGlobal = NULL;
static double* ptrLowResGlobal = NULL;

void init (double** sphere, double** extracted_region, double** ptrHighRes, double** ptrLowRes, double** rotation_matrix, double** ptrEvecOut, double ** ptrEvalsOut) {
    //
    // Create sphere mask
    //
    *extracted_region = malloc( (sizeof (double)) * SPHERE_ARRAY_SIZE);
    *sphere = malloc( (sizeof (double)) * SPHERE_ARRAY_SIZE);
    createSphereMask(*sphere);
    //
    // Read input images
    //
    if (ptrHighResGlobal == NULL) {
        ptrHighResGlobal = readHighResImage();
    }
    if (ptrLowResGlobal == NULL) {
        ptrLowResGlobal = readLowResImage();
    }
    *ptrHighRes = ptrHighResGlobal;
    *ptrLowRes = ptrLowResGlobal;
    //
    // Init rotation matrix
    //
    *rotation_matrix = malloc (sizeof(double) * 9);
    (*rotation_matrix)[0] = 0.93698197;
    (*rotation_matrix)[1] = -0.00553534;
    (*rotation_matrix)[2] = 0.34933386;
    (*rotation_matrix)[3] = 0.00427267;
    (*rotation_matrix)[4] = 0.99998126;
    (*rotation_matrix)[5] = 0.004385;
    (*rotation_matrix)[6] = -0.34935159;
    (*rotation_matrix)[7] = -0.00261608;
    (*rotation_matrix)[8] = 0.93698806;
    //
    // Allocate space for the output eigen vectors
    //
    *ptrEvecOut = calloc (sizeof(double), 3*3*LOW_RES_SIZE);
    *ptrEvalsOut = calloc(sizeof(double), 3*LOW_RES_SIZE);

#ifdef DEBUG
    region_extraction_debug_init();
    //fit_ellipsoid_debug_init();
#endif
}

void deInit (double* sphere, double* extracted_region, double* ptrHighRes, double* ptrLowRes, double* rotation_matrix, double* ptrEvecOut, double *ptrEvalsOut, bool generate_ground_truth) {
    if (generate_ground_truth) {
        FILE *fd = fopen("ground_truth.txt","w");

        for (int k=0; k<LOW_RES_D3; k++) {
            for (int j=0; j<LOW_RES_D2; j++) {
                for (int i=0; i<LOW_RES_D1; i++) {
                    int temp = i + j*LOW_RES_D1 + k*LOW_RES_D1*LOW_RES_D2;
                    fprintf(fd, "%d, %d, %d, " /* index into LR image */
                            "%.20f, %.20f, %.20f, "   /* e vals */
                            "%.20f, %.20f, %.20f, "   /* e vecs*/
                            "%.20f, %.20f, %.20f, "
                            "%.20f, %.20f, %.20f\n",
                            i, j, k, 
                            ptrEvalsOut[3*temp+0], ptrEvalsOut[3*temp+1], ptrEvalsOut[3*temp+2],
                            ptrEvecOut[9*temp+0],ptrEvecOut[9*temp+1],ptrEvecOut[9*temp+2],
                            ptrEvecOut[9*temp+3],ptrEvecOut[9*temp+4],ptrEvecOut[9*temp+5],
                            ptrEvecOut[9*temp+6],ptrEvecOut[9*temp+7],ptrEvecOut[9*temp+8]);
                }
            }
       }
       fclose(fd);
       printf("Generated ground truth file\n");
    }

    free(sphere);
    free(extracted_region);
    //free(ptrHighRes); // we always use the same highres
    //free(ptrLowRes);
    free(rotation_matrix);
    free(ptrEvecOut);
    free(ptrEvalsOut);
 
#ifdef DEBUG
    region_extraction_debug_deinit();
    //fit_ellipsoid_debug_deinit();
#endif
}

void kernel_basic (double* sphere, double* extracted_region, double* ptrHighRes, double* ptrLowRes, double* rotation_matrix, double* ptrEvecOut, double *ptrEvalsOut) {
    double xT, yT, zT, xC, yC, zC;

    xT = 5.386915;
    yT = 5.253300;
    zT = 16.591884;
    xC = 37.0230010226;
    yC = 22.5910006240;
    zC = 45.8790012673;

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
    r20 = rotation_matrix[6];
    r21 = rotation_matrix[7];
    r22 = rotation_matrix[8];
      
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
                if (ptrLowRes[ii_lr] > 0.5)
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
                    if ((tx_hr <= nx_hr-half) && (tx_hr > nx_hr-one)) {
                        tx_hr = nx_hr-one;
                    }
                    else if ((ty_hr <= ny_hr-half) && (ty_hr > ny_hr-one)) {
                        ty_hr = ny_hr-one;
                    }
                    else if ((tz_hr <= nz_hr-half) && (tz_hr > nz_hr-one)) {
                        tz_hr = nz_hr-one;
                    }

                    i_hr = (int) tx_hr;
                    j_hr = (int) ty_hr;
                    k_hr = (int) tz_hr;
 
                    region_extraction(i_hr, j_hr, k_hr, sphere, extracted_region, ptrHighRes);

                    /* // compute fabric */
                    /* double mils[NUM_DIRECTIONS]; */
                    /* mil2(extracted_region, SPHERE_NDIM, mils); */
                    /* //print_vector(mils, NUM_DIRECTIONS); */
                    /* double Q[3][3]; */
                    /* fit_ellipsoid_mils(mils, (double (*)[3][3])Q); */

                    /* eigen3(Q, &ptrEvecOut[ii_lr*9], &ptrEvalsOut[ii_lr*3]); */
                }
            }
        }
    } /* main loop */
}

// Optimized version
void kernel_opt1 (double* sphere, double* extracted_region, double* ptrHighRes, double* ptrLowRes, double* rotation_matrix, double* ptrEvecOut, double *ptrEvalsOut) {
    double xT, yT, zT, xC, yC, zC;

    xT = 5.386915;
    yT = 5.253300;
    zT = 16.591884;
    xC = 37.0230010226;
    yC = 22.5910006240;
    zC = 45.8790012673;

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
    r20 = rotation_matrix[6];
    r21 = rotation_matrix[7];
    r22 = rotation_matrix[8];

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
                if (ptrLowRes[ii_lr] > 0.5)
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
                    if ((tx_hr <= nx_hr-half) && (tx_hr > nx_hr-one)) {
                        tx_hr = nx_hr-one;
                    }
                    else if ((ty_hr <= ny_hr-half) && (ty_hr > ny_hr-one)) {
                        ty_hr = ny_hr-one;
                    }
                    else if ((tz_hr <= nz_hr-half) && (tz_hr > nz_hr-one)) {
                        tz_hr = nz_hr-one;
                    }

                    i_hr = (int) tx_hr;
                    j_hr = (int) ty_hr;
                    k_hr = (int) tz_hr;
                    
                    region_extraction_opt1(i_hr, j_hr, k_hr, sphere, extracted_region, ptrHighRes);

                    /* // compute fabric */
                    /* double mils[NUM_DIRECTIONS]; */
                    /* mil2(extracted_region, SPHERE_NDIM, mils); */
                    /* //print_vector(mils, NUM_DIRECTIONS); */

                    /* double Q[3][3]; */
                    /* fit_ellipsoid_mils(mils, Q); */

                    /* eigen3(Q, &ptrEvecOut[ii_lr*9], &ptrEvalsOut[ii_lr*3]); */
                }
            }
        }
    } /* main loop */
}
