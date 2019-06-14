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
#include <immintrin.h>

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
static long region_extraction_call_count = 0;
void region_extraction_debug_init(void)
    {
      region_extraction_flop_count = 0;
      region_extraction_call_count = 0;
    }
void region_extraction_debug_deinit(void)
    {
      printf("[regions] flop count = %ld\n", region_extraction_flop_count);
      printf("[regions] call count = %ld\n", region_extraction_call_count);
      long int readwrite_bytes = SPHERE_NDIM*SPHERE_NDIM*SPHERE_NDIM*3*(sizeof(double));
      printf("[regions] read & write  %ld bytes \n", readwrite_bytes);
      printf("[regions] opt intensity  %f flops/bytes \n", (float) region_extraction_flop_count/(readwrite_bytes*region_extraction_call_count));
    }
#endif

///
/// Skeleton for region extraction
///
void region_extraction (int i_hr, int j_hr, int k_hr, double *sphere, double *extracted_region, double *ptrHighRes) {

    int ihr_min, jhr_min, khr_min;
    int ii;
    int ihr, jhr, khr, ii_hr;

#ifdef DEBUG
    region_extraction_call_count += 1;
#endif
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
                ii =  i + j*SPHERE_NDIM + k*SPHERE_NDIM*SPHERE_NDIM;
                ii_hr =  ihr + jhr*HIGH_RES_D1 + khr*HIGH_RES_D1*HIGH_RES_D2;
    
                // 
                // if (ii_hr < HIGH_RES_SIZE) {
                if (ihr < HIGH_RES_D1 && jhr < HIGH_RES_D2 && khr < HIGH_RES_D3) {
                    extracted_region[ii] = sphere[ii] * ptrHighRes[ii_hr];
                }

#ifdef DEBUG
    region_extraction_flop_count += 1;
#endif
                else {
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
    double hr_voxel4, hr_voxel5, hr_voxel6, hr_voxel7;
    double sp_voxel0, sp_voxel1, sp_voxel2, sp_voxel3;
    double sp_voxel4, sp_voxel5, sp_voxel6, sp_voxel7;
    double ex_voxel0, ex_voxel1, ex_voxel2, ex_voxel3;
    double ex_voxel4, ex_voxel5, ex_voxel6, ex_voxel7;
    int NUM_ACC = 4;

    // find min for sphere
    ihr_min = i_hr - SPHERE_HALF_NDIM;
    jhr_min = j_hr - SPHERE_HALF_NDIM;
    khr_min = k_hr - SPHERE_HALF_NDIM;

    for (int k=0; k < SPHERE_NDIM; k++) {
        khr = k + khr_min;
        for (int j=0; j < SPHERE_NDIM; j++) {
            jhr = j + jhr_min;
            for (int i=0; i < SPHERE_NDIM; i+=NUM_ACC) {
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
                hr_voxel0 = 0.0;
                hr_voxel1 = 0.0;
                hr_voxel2 = 0.0;
                hr_voxel3 = 0.0;
                //
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

    }}}
}


static void print_simd(__m256d a) {
    printf("(%f, %f, %f, %f)\n", a[0], a[1], a[2], a[3]);
}


///
/// Skeleton for region extraction
///
void region_extraction_simd (int i_hr, int j_hr, int k_hr, double *sphere, double *extracted_region, double *ptrHighRes) {

    int ihr_min, jhr_min, khr_min;
    int ihr_max, jhr_max, khr_max;
    int i, j, k, ii;
    int ihr, jhr, khr, ii_hr;

    __m256d sphere_vec, highres_vec, extract_vec;

    // find min for sphere
    ihr_min = i_hr - SPHERE_HALF_NDIM;
    jhr_min = j_hr - SPHERE_HALF_NDIM;
    khr_min = k_hr - SPHERE_HALF_NDIM;
    // find max for sphere
    ihr_max = i_hr + SPHERE_HALF_NDIM;
    jhr_max = j_hr + SPHERE_HALF_NDIM;
    khr_max = k_hr + SPHERE_HALF_NDIM;
    //
    if (ihr_max > HIGH_RES_D1) {
       ihr_max = HIGH_RES_D1;
    }
    if (jhr_max > HIGH_RES_D2) {
       jhr_max = HIGH_RES_D2;
    }
    if (khr_max > HIGH_RES_D3) {
       khr_max = HIGH_RES_D3;
    }

    extract_vec = _mm256_set1_pd(0.0);
    for (k=0; k < SPHERE_NDIM; k++) {
        for (j=0; j < SPHERE_NDIM; j++) {
            for (i=0; i < SPHERE_NDIM; i+=4) {
                ii = i + j*SPHERE_NDIM + k*SPHERE_NDIM*SPHERE_NDIM;
                _mm256_store_pd(extracted_region + ii, extract_vec);
                
            }
        }
    }
    for (khr=khr_min; khr < khr_max; khr++) {
        k = khr - khr_min;
        for (jhr=jhr_min; jhr < jhr_max; jhr++) {
            j = jhr - jhr_min;
            for (ihr=ihr_min; ihr < ihr_max-4; ihr+=4) {
                i = ihr - ihr_min;
                // calculate index
                ii =  i + j*SPHERE_NDIM + k*SPHERE_NDIM*SPHERE_NDIM;
                ii_hr =  ihr + jhr*HIGH_RES_D1 + khr*HIGH_RES_D1*HIGH_RES_D2;
                // printf("%d, %d, %d\n", i, j, k);
                // 
                sphere_vec = _mm256_load_pd(sphere + ii);
                highres_vec = _mm256_loadu_pd(ptrHighRes + ii_hr);
                extract_vec = _mm256_mul_pd(sphere_vec, highres_vec);

                // store
                _mm256_store_pd(extracted_region + ii, extract_vec);
            }
            for (;ihr<ihr_max; ihr++) { /* need clean-up */
                i = ihr - ihr_min;
                // calculate index
                ii =  i + j*SPHERE_NDIM + k*SPHERE_NDIM*SPHERE_NDIM;
                ii_hr =  ihr + jhr*HIGH_RES_D1 + khr*HIGH_RES_D1*HIGH_RES_D2;
                extracted_region[ii] = sphere[ii] * ptrHighRes[ii_hr];
            }
        }
    }
}
