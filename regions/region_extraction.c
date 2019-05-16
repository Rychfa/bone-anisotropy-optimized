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
#include "region_extraction.h"
#include "reader.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#undef DEBUG
//#define DEBUG  
///
/// Write a sphere mask
///
void createSphereMask(int *sphere) {

    double radius = SPHERE_DIAMETER/2.0;
    double distance, xc, yc, zc;
    double x, y, z;
    int ii;

    xc = radius;
    yc = radius;
    zc = radius;

    double halfvoxel = HIGH_RES_VOXEL_SIZE/2.0;

    // loop over all voxels and set 1 in the sphere
    for (int k=0; k < SPHERE_NDIM; k++){
        z = k*HIGH_RES_VOXEL_SIZE + halfvoxel;
        for (int j=0; j < SPHERE_NDIM; j++){
            y = j*HIGH_RES_VOXEL_SIZE + halfvoxel;
            for (int i=0; i < SPHERE_NDIM; i++){
                x = i*HIGH_RES_VOXEL_SIZE + halfvoxel;
 
                distance = sqrt(pow(x-xc, 2) + pow(y-yc, 2) + pow(z-zc, 2));
                ii = i + j*SPHERE_NDIM + k*SPHERE_NDIM*SPHERE_NDIM;
                if (distance <= radius) {
                    sphere[ii] = 1; 
                } else {
                    sphere[ii] = 0; 
                }
            }
        }
    }

}

///
/// Skeleton for region extraction
///
void region_extraction (int i_hr, int j_hr, int k_hr, int *sphere, int *extracted_region, int *ptrHighRes) {

    int ihr_min, jhr_min, khr_min;
    int ii;
    int ihr, jhr, khr, ii_hr;
#ifdef DEBUG
    int flop_count = 0;
#endif

    // find min for sphere
    ihr_min = i_hr - SPHERE_HALF_NDIM;
    jhr_min = j_hr - SPHERE_HALF_NDIM;
    khr_min = k_hr - SPHERE_HALF_NDIM;

#ifdef DEBUG
    flop_count += 3;
#endif

    for (int k=0; k < SPHERE_NDIM; k++) {
        khr = k + khr_min;
#ifdef DEBUG
    flop_count += 1;
#endif
        for (int j=0; j < SPHERE_NDIM; j++) {
            jhr = j + jhr_min;
#ifdef DEBUG
    flop_count += 1;
#endif
            for (int i=0; i < SPHERE_NDIM; i++) {
                ihr = i + ihr_min;
#ifdef DEBUG
    flop_count += 1;
#endif
                // calculate index
                ii = i + j*SPHERE_NDIM + k*SPHERE_NDIM*SPHERE_NDIM;
                ii_hr = ihr + jhr*HIGH_RES_D1 + khr*HIGH_RES_D1*HIGH_RES_D2;
    
#ifdef DEBUG
    flop_count += 10;
#endif
                //
                if (ii_hr < HIGH_RES_SIZE) {
                    extracted_region[ii] = sphere[ii] * ptrHighRes[ii_hr];
#ifdef DEBUG
    flop_count += 1;
#endif
                } else {
                    extracted_region[ii] = 0;
                }
    }}}

#ifdef DEBUG
    printf("region_extraction: flop_count %d flops \n", flop_count);
    long int readwrite_bytes = SPHERE_NDIM*SPHERE_NDIM*SPHERE_NDIM*3*(sizeof(int));
    printf("region_extraction: read & write  %ld bytes \n", readwrite_bytes);
    printf("region_extraction: opt intensity  %f flops/bytes \n", (float) flop_count/readwrite_bytes);
#endif

}


///
/// Skeleton for region extraction
///
void region_extraction_opt1 (int i_hr, int j_hr, int k_hr, int *sphere, int *extracted_region, int *ptrHighRes) {

    int ihr_min, jhr_min, khr_min;
    int isp, jsp, ksp, jksp, ii;
    int ihr, jhr, khr, jkhr, ii_hr;
    int SPHERE_NDIM_SQR = SPHERE_NDIM*SPHERE_NDIM;
    int HIGH_RES_D1D2 = HIGH_RES_D1*HIGH_RES_D2;

    int inb, jnb, knb;
    int inb_start, jnb_start, knb_start;
    int NBLOCKS = 9;
    int BLOCK_SIZE = SPHERE_NDIM/NBLOCKS;
#ifdef DEBUG
    int flop_count = 0;
#endif

    // find min for sphere
    ihr_min = i_hr - SPHERE_HALF_NDIM;
    jhr_min = j_hr - SPHERE_HALF_NDIM;
    khr_min = k_hr - SPHERE_HALF_NDIM;
#ifdef DEBUG
    flop_count += 6;
#endif

    for (int knb=0; knb < NBLOCKS; knb++) {
        knb_start = knb*NBLOCKS;
        for (int jnb=0; jnb < NBLOCKS; jnb++) {
            jnb_start = jnb*NBLOCKS;
            for (int inb=0; inb < NBLOCKS; inb++) {
                inb_start = inb*NBLOCKS;

                for (int k=knb_start; k < knb_start+BLOCK_SIZE; k++) {
                    ksp = k*SPHERE_NDIM_SQR;
                    khr = (k + khr_min)*HIGH_RES_D1D2;
#ifdef DEBUG
    flop_count += 3;
#endif
                    for (int j=jnb_start; j < jnb_start+BLOCK_SIZE; j++) {
                        jksp = j*SPHERE_NDIM + ksp ;
                        jkhr = (j + jhr_min)*HIGH_RES_D1 + khr;

#ifdef DEBUG
    flop_count += 5;
#endif
                        for (int i=inb_start; i < inb_start+BLOCK_SIZE; i+=2) {
                            // calculate index
                            ii = i + jksp;
                            ii_hr = (i + ihr_min) + jkhr;
#ifdef DEBUG
    flop_count += 3;
#endif
                            extracted_region[ii] = 0;
                            extracted_region[ii+1] = 0;
                            // extract region with a sphere mask
                            if (ii_hr < HIGH_RES_SIZE) {
                                extracted_region[ii] = sphere[ii] * ptrHighRes[ii_hr];
#ifdef DEBUG
    flop_count += 1;
#endif
                            }
                            if (ii_hr+1 < HIGH_RES_SIZE) {
                                extracted_region[ii+1] = sphere[ii+1] * ptrHighRes[ii_hr+1];
#ifdef DEBUG
    flop_count += 1;
#endif
                            }
        }}} // mini block loop
    }}} // Block number loop

#ifdef DEBUG
    printf("region_extraction: flop_count %d flops \n", flop_count);
    long int readwrite_bytes = BLOCK_SIZE*BLOCK_SIZE*BLOCK_SIZE*3*(sizeof(int));
    printf("region_extraction: read & write  %ld bytes \n", readwrite_bytes);
    printf("region_extraction: opt intensity  %f flops/bytes \n", (float) flop_count/readwrite_bytes);
#endif
}

///
/// Skeleton for region extraction
///
void region_extraction_opt2 (int i_hr, int j_hr, int k_hr, int *sphere, int *extracted_region, int *ptrHighRes) {

    int ihr_min, jhr_min, khr_min;
    int isp, jsp, ksp, jksp, ii;
    int ihr, jhr, khr, jkhr, ii_hr;
    int SPHERE_NDIM_SQR = SPHERE_NDIM*SPHERE_NDIM;
    int HIGH_RES_D1D2 = HIGH_RES_D1*HIGH_RES_D2;

#ifdef DEBUG
    int flop_count = 0;
#endif

    // find min for sphere
    ihr_min = i_hr - SPHERE_HALF_NDIM;
    jhr_min = j_hr - SPHERE_HALF_NDIM;
    khr_min = k_hr - SPHERE_HALF_NDIM;
#ifdef DEBUG
    flop_count += 5;
#endif

                for (int k=0; k < SPHERE_NDIM; k++) {
                    ksp = k*SPHERE_NDIM_SQR;
                    khr = (k + khr_min)*HIGH_RES_D1D2;
#ifdef DEBUG
    flop_count += 3;
#endif
                    for (int j=0; j < SPHERE_NDIM; j++) {
                        jksp = j*SPHERE_NDIM + ksp ;
                        jkhr = (j + jhr_min)*HIGH_RES_D1 + khr;

#ifdef DEBUG
    flop_count += 5;
#endif
                        for (int i=0; i < SPHERE_NDIM; i+=2) {
                            // calculate index
                            ii = i + jksp;
                            ii_hr = (i + ihr_min) + jkhr;
#ifdef DEBUG
    flop_count += 3;
#endif
                            extracted_region[ii] = 0;
                            extracted_region[ii+1] = 0;
                            // extract region with a sphere mask
                            if (ii_hr < HIGH_RES_SIZE) {
                                extracted_region[ii] = sphere[ii] * ptrHighRes[ii_hr];
#ifdef DEBUG
    flop_count += 1;
#endif
                            }
                            if (ii_hr+1 < HIGH_RES_SIZE) {
                                extracted_region[ii+1] = sphere[ii+1] * ptrHighRes[ii_hr+1];
#ifdef DEBUG
    flop_count += 1;
#endif
                            }
    }}} // Block number loop

#ifdef DEBUG
    printf("region_extraction: flop_count %d flops \n", flop_count);
    long int readwrite_bytes = SPHERE_NDIM*SPHERE_NDIM*SPHERE_NDIM*3*(sizeof(int));
    printf("region_extraction: read & write  %ld bytes \n", readwrite_bytes);
    printf("region_extraction: opt intensity  %f flops/bytes \n", (float) flop_count/readwrite_bytes);
#endif
}
