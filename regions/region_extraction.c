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
                if (ii_hr < HIGH_RES_SIZE){
                    extracted_region[ii] = sphere[ii] * ptrHighRes[ii_hr];
                } else {
                    extracted_region[ii] = 0;
                }

    }}}
}
