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

///
/// Write a sphere mask
///
void createSphereMask(double *sphere) {

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
                    sphere[ii] = 1.0; 
                } else {
                    sphere[ii] = 0.0; 
                }
            }
        }
    }

}
#ifdef DEBUG
static long region_extraction_flop_count = 0;
void region_extraction_debug_init(void)
    {
      region_extraction_flop_count = 0;
    }
void region_extraction_debug_deinit(void)
    {
      printf("[regions] flop count = %ld\n", region_extraction_flop_count);
      long int readwrite_bytes = SPHERE_NDIM*SPHERE_NDIM*SPHERE_NDIM*3*(sizeof(int));
      printf("region_extraction: read & write  %ld bytes \n", readwrite_bytes);
      printf("region_extraction: opt intensity  %f flops/bytes \n", (float) region_extraction_flop_count/readwrite_bytes);
    }
#endif

///
/// Skeleton for region extraction
///
void region_extraction (int i_hr, int j_hr, int k_hr, double *sphere, double *extracted_region, double *ptrHighRes) {

    int ihr_min, jhr_min, khr_min;
    int ii;
    int ihr, jhr, khr, ii_hr;

    // find min for sphere
    ihr_min = i_hr - SPHERE_HALF_NDIM;
    jhr_min = j_hr - SPHERE_HALF_NDIM;
    khr_min = k_hr - SPHERE_HALF_NDIM;

    for (int k=0; k < SPHERE_NDIM; k++) {
        khr = k + khr_min;
        for (int j=0; j < SPHERE_NDIM; j++) {
            jhr = j + jhr_min;
            for (int i=0; i < SPHERE_NDIM; i++) {
                ihr = i + ihr_min;
                // calculate index
                ii = i + j*SPHERE_NDIM + k*SPHERE_NDIM*SPHERE_NDIM;
                ii_hr = ihr + jhr*HIGH_RES_D1 + khr*HIGH_RES_D1*HIGH_RES_D2;
    
                //
                // sign(ii_hr - HIGH_RES_SIZE);
                if (ii_hr < HIGH_RES_SIZE) {
                    extracted_region[ii] = sphere[ii] * ptrHighRes[ii_hr];
#ifdef DEBUG
    region_extraction_flop_count += 1;
#endif
                } else {
                    extracted_region[ii] = 0;
                }
    }}}
}



///
/// Skeleton for region extraction
///
void region_extraction_opt1 (int i_hr, int j_hr, int k_hr, double *sphere, double *extracted_region, double *ptrHighRes) {

    int ihr_min, jhr_min, khr_min;
    int ii;
    int ihr, jhr, khr, ii_hr;
    double hr_voxel0, hr_voxel1, hr_voxel2, hr_voxel3;
    double sp_voxel0, sp_voxel1, sp_voxel2, sp_voxel3;
    double ex_voxel0, ex_voxel1, ex_voxel2, ex_voxel3;

    // find min for sphere
    ihr_min = i_hr - SPHERE_HALF_NDIM;
    jhr_min = j_hr - SPHERE_HALF_NDIM;
    khr_min = k_hr - SPHERE_HALF_NDIM;

    for (int k=0; k < SPHERE_NDIM; k++) {
        khr = k + khr_min;
        for (int j=0; j < SPHERE_NDIM; j++) {
            jhr = j + jhr_min;
            for (int i=0; i < SPHERE_NDIM; i+=4) {
                ihr = i + ihr_min;
                // calculate index
                ii = i + j*SPHERE_NDIM + k*SPHERE_NDIM*SPHERE_NDIM;
                ii_hr = ihr + jhr*HIGH_RES_D1 + khr*HIGH_RES_D1*HIGH_RES_D2;
                // 
                sp_voxel0 = sphere[ii];
                sp_voxel1 = sphere[ii+1];
                sp_voxel2 = sphere[ii+2];
                sp_voxel3 = sphere[ii+3];
                //
                hr_voxel0 = 0;
                hr_voxel1 = 0;
                hr_voxel2 = 0;
                hr_voxel3 = 0;
                // 0
                if (ii_hr < HIGH_RES_SIZE) {
                    hr_voxel0 = ptrHighRes[ii_hr];
                } //else { hr_voxel0 = 0; }
                // 1
                if (ii_hr + 1 < HIGH_RES_SIZE) {
                    hr_voxel1 = ptrHighRes[ii_hr + 1];
                } //else { hr_voxel1 = 0; }
                // 2
                if (ii_hr + 2 < HIGH_RES_SIZE) {
                    hr_voxel2 = ptrHighRes[ii_hr + 2];
                } //else { hr_voxel2 = 0; }
                // 3
                if (ii_hr + 3 < HIGH_RES_SIZE) {
                    hr_voxel3 = ptrHighRes[ii_hr + 3];
                } //else { hr_voxel = 0; }

                //
                ex_voxel0 = sp_voxel0 * hr_voxel0;
                ex_voxel1 = sp_voxel1 * hr_voxel1;
                ex_voxel2 = sp_voxel2 * hr_voxel2;
                ex_voxel3 = sp_voxel3 * hr_voxel3;
                //
                extracted_region[ii] = ex_voxel0;
                extracted_region[ii+1] = ex_voxel1;
                extracted_region[ii+2] = ex_voxel2;
                extracted_region[ii+3] = ex_voxel3;
#ifdef DEBUG
    region_extraction_flop_count += 4;
#endif

    }}}
}



///
/// Skeleton for region extraction
///
void region_extraction_opt2 (int i_hr, int j_hr, int k_hr, double *sphere, double *extracted_region, double *ptrHighRes) {

    int ihr_min, jhr_min, khr_min;
    int isp, jsp, ksp, jksp, ii;
    int ihr, jhr, khr, jkhr, ii_hr;

    int inb, jnb, knb;
    int inb_start, jnb_start, knb_start;
    int BLOCK_SIZE = 16;
    int NBLOCKS = SPHERE_NDIM/BLOCK_SIZE;

    // find min for sphere
    ihr_min = i_hr - SPHERE_HALF_NDIM;
    jhr_min = j_hr - SPHERE_HALF_NDIM;
    khr_min = k_hr - SPHERE_HALF_NDIM;

    for (int knb=0; knb < NBLOCKS; knb++) {
        knb_start = knb*BLOCK_SIZE;
        for (int jnb=0; jnb < NBLOCKS; jnb++) {
            jnb_start = jnb*BLOCK_SIZE;
            for (int inb=0; inb < NBLOCKS; inb++) {
                inb_start = inb*BLOCK_SIZE;
                //
                for (int k=knb_start; k < knb_start+BLOCK_SIZE; k++) {
                    khr = k + khr_min;
                    for (int j=jnb_start; j < jnb_start+BLOCK_SIZE; j++) {
                        jhr = j + jhr_min;
                        for (int i=inb_start; i < inb_start+BLOCK_SIZE; i++) {
                            ihr = i + ihr_min;
                            // calculate index
                            ii = i + j*SPHERE_NDIM + k*SPHERE_NDIM*SPHERE_NDIM;
                            ii_hr = ihr + jhr*HIGH_RES_D1 + khr*HIGH_RES_D1*HIGH_RES_D2;
                            //
                            // extract region with a sphere mask
                            if (ii_hr < HIGH_RES_SIZE) {
                                extracted_region[ii] = sphere[ii] * ptrHighRes[ii_hr];
#ifdef DEBUG
    region_extraction_flop_count += 1;
#endif
                            } else {
                                extracted_region[ii] = 0;
                            }
        }}} // mini block loop
    }}} // Block number loop
}



///
/// blocking, loop unrolling and scalar replacement
///
void region_extraction_opt3 (int i_hr, int j_hr, int k_hr, double *sphere, double *extracted_region, double *ptrHighRes) {

    int ihr_min, jhr_min, khr_min;
    int isp, jsp, ksp, jksp, ii;
    int ihr, jhr, khr, jkhr, ii_hr;

    int inb, jnb, knb;
    int inb_start, jnb_start, knb_start;
    double hr_voxel0, hr_voxel1, hr_voxel2, hr_voxel3;
    double sp_voxel0, sp_voxel1, sp_voxel2, sp_voxel3;
    double ex_voxel0, ex_voxel1, ex_voxel2, ex_voxel3;
    int BLOCK_SIZE = 16;
    int NBLOCKS = SPHERE_NDIM/BLOCK_SIZE;

    // find min for sphere
    ihr_min = i_hr - SPHERE_HALF_NDIM;
    jhr_min = j_hr - SPHERE_HALF_NDIM;
    khr_min = k_hr - SPHERE_HALF_NDIM;

    for (int knb=0; knb < NBLOCKS; knb++) {
        knb_start = knb*BLOCK_SIZE;
        for (int jnb=0; jnb < NBLOCKS; jnb++) {
            jnb_start = jnb*BLOCK_SIZE;
            for (int inb=0; inb < NBLOCKS; inb++) {
                inb_start = inb*BLOCK_SIZE;
                //
                for (int k=knb_start; k < knb_start+BLOCK_SIZE; k++) {
                    khr = k + khr_min;
                    for (int j=jnb_start; j < jnb_start+BLOCK_SIZE; j++) {
                        jhr = j + jhr_min;
                        for (int i=inb_start; i < inb_start+BLOCK_SIZE; i+=4) {
                            ihr = i + ihr_min;
                            // calculate index
                            ii = i + j*SPHERE_NDIM + k*SPHERE_NDIM*SPHERE_NDIM;
                            ii_hr = ihr + jhr*HIGH_RES_D1 + khr*HIGH_RES_D1*HIGH_RES_D2;
                            //

                // 
                sp_voxel0 = sphere[ii];
                sp_voxel1 = sphere[ii+1];
                sp_voxel2 = sphere[ii+2];
                sp_voxel3 = sphere[ii+3];
                //
                hr_voxel0 = 0;
                hr_voxel1 = 0;
                hr_voxel2 = 0;
                hr_voxel3 = 0;
                // 0
                if (ii_hr < HIGH_RES_SIZE) {
                    hr_voxel0 = ptrHighRes[ii_hr];
                } 
                // 1
                if (ii_hr + 1 < HIGH_RES_SIZE) {
                    hr_voxel1 = ptrHighRes[ii_hr + 1];
                } 
                // 2
                if (ii_hr + 2 < HIGH_RES_SIZE) {
                    hr_voxel2 = ptrHighRes[ii_hr + 2];
                } 
                // 3
                if (ii_hr + 3 < HIGH_RES_SIZE) {
                    hr_voxel3 = ptrHighRes[ii_hr + 3];
                } 

                //
                ex_voxel0 = sp_voxel0 * hr_voxel0;
                ex_voxel1 = sp_voxel1 * hr_voxel1;
                ex_voxel2 = sp_voxel2 * hr_voxel2;
                ex_voxel3 = sp_voxel3 * hr_voxel3;
                //
                extracted_region[ii] = ex_voxel0;
                extracted_region[ii+1] = ex_voxel1;
                extracted_region[ii+2] = ex_voxel2;
                extracted_region[ii+3] = ex_voxel3;

#ifdef DEBUG
    region_extraction_flop_count += 4;
#endif

        }}} // mini block loop
    }}} // Block number loop
}
