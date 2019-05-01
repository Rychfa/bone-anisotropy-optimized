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
#include "coord_map.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

///
/// Write a sphere mask
///
void createSphereMask(double diameter, double voxelsize) {

int ndim = (int) diameter/voxelsize + 1;
int sphere[ndim][ndim][ndim];
double radius = diameter/2.0;
double distance, xc, yc, zc;
double x, y, z;

xc = radius;
yc = radius;
zc = radius;

double halfvoxel = voxelsize/2.0;

// loop over all voxels and set 1 in the sphere
for (int k=0; k < ndim; k++){
    z = k*voxelsize + halfvoxel;
    for (int j=0; j < ndim; j++){
        y = j*voxelsize + halfvoxel;
        for (int i=0; i < ndim; i++){
            x = i*voxelsize + halfvoxel;

            distance = sqrt(pow(x-xc, 2) + pow(y-yc, 2) + pow(z-zc, 2));
            if (distance <= radius){ 
                sphere[k][j][i] = 1; }
        }
    }
}

//return sphere;

}

///
/// Skeleton for coordinate mapping.
///
void coordMap (double ***lr_image, double ***hr_image, double sphere_mask) {

double r00, r01, r02, r10, r11, r12, r20, r21, r22;
double voxel_size_lr, voxel_size_hr;
double half_voxel_size_lr, half_voxel_size_hr;
double xC, yC, zC, xT, yT, zT;
double xDr00, xDr10, xDr20, yDr01, yDr11, yDr21, zDr02, zDr12, zDr22;
double tx_hr, ty_hr, tz_hr;
double nx_hr, ny_hr, nz_hr;
double nx_lr, ny_lr, nz_lr;
double nx_sp, ny_sp, nz_sp;
double x_lr, y_lr, z_lr;
double zero, half, one;
int i_hr0, j_hr0, k_hr0;
int imin, jmin, kmin;
int radius_i, radius_j, radius_k;

zero = 0.0;
half = 0.5;
one = 1.0;

// XCT - High resolution
// image dimension
nx_hr = 0.0;
ny_hr = 0.0;
nz_hr = 0.0;
// voxel_size
voxel_size_hr = 0.082;
half_voxel_size_hr = 0.041;

// QCT - Low resolution
// image dimension
nx_lr = 0.0;
ny_lr = 0.0;
nz_lr = 0.0;
// voxel_size
voxel_size_lr = 3.0;
half_voxel_size_lr = 1.5;

// Sphere mask
// image dimension
nx_sp = 0.0;
ny_sp = 0.0;
nz_sp = 0.0;
// sphere radius
radius_i = 0.0;
radius_j = 0.0;
radius_k = 0.0;


// Rotation matrix
r00 = 1.0; r01 = 0.0; r02 = 0.0;
r10 = 0.0; r11 = 1.0; r12 = 0.0;
r20 = 0.0; r21 = 0.0; r22 = 1.0;
// Center of rotation
xC = 0.0;
yC = 0.0;
zC = 0.0;
// translation
xT = 0.0;
yT = 0.0;
zT = 0.0;


// loop over all femur voxels
for (int k_lr=0; k_lr < nz_lr; k_lr++)
{
    // calculate vector from the center of the image 
    // to the center of this voxel
    z_lr = (float) k_lr*voxel_size_lr - zC + half_voxel_size_lr;
    // multiply the vector with the rotation matrix
    zDr02 = z_lr * r02;
    zDr12 = z_lr * r12; 
    zDr22 = z_lr * r22;

    for (int j_lr=0; j_lr < ny_lr; j_lr++) 
    {
        //calculate vector from the center of the image 
        //to the center of this voxel
        y_lr = (float) j_lr*voxel_size_lr - yC + half_voxel_size_lr;
        //multiply the vector with the rotation matrix
        yDr01 = y_lr * r01;
        yDr11 = y_lr * r11;
        yDr21 = y_lr * r21;

        for (int i_lr=0; i_lr < nx_lr; i_lr++) 
        {
           // calculate vector from the center of the image 
           // to the center of this voxel
           x_lr = (float) i_lr*voxel_size_lr - xC + half_voxel_size_lr;
           // check if this voxel inside the FE mask
           if (lr_image[k_lr][j_lr][i_lr] > 0.5) 
           {
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
             i_hr0 = (int) tx_hr;
             j_hr0 = (int) ty_hr;
             k_hr0 = (int) tz_hr;
             
             // extract a sphere region
             //extracted_sphere_region = region_extraction(i_hr0, j_hr0, k_hr0);
             // compute fabric
             //evec, eval = mil(extracted_sphere_region);
 

          }  

        }
    }
}

}
