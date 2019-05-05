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
#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


///
/// Calculate rotation matrix
///
void calculateRotationMatrix(double r[9], double ax, double ay, double az) {
    double caX, caY, caZ, saX, saY, saZ;
    
    caX = cos(ax);
    caY = cos(ay);
    caZ = cos(az);
    saX = sin(ax);
    saY = sin(ay);
    saZ = sin(az);
    
    r[0] = caY*caZ;
    r[1] = saX*saY*caZ-caX*saZ;
    r[2] = caX*saY*caZ+saX*saZ;
    r[3] = caY*saZ;
    r[4] = saX*saY*saZ+caX*caZ;
    r[5] = caX*saY*saZ-saX*caZ;
    r[6] = -saY;
    r[7] = saX*caY;
    r[8] = caX*caY;

    printf("rotation_matrix: r0 %.3f %.3f %.3f \n", r[0], r[1], r[2]);
    printf("rotation_matrix: r1 %.3f %.3f %.3f \n", r[3], r[4], r[5]);
    printf("rotation_matrix: r2 %.3f %.3f %.3f \n", r[6], r[7], r[8]);
}

void writeVTK(int *image, double voxelsize, int ndim) {
    printf("writeVTK: enter\n");
    // output debug
    FILE *fptr;
    double coord;
    int ndimp1 = ndim+1;
    int ii;

    
    printf("writeVTK: start writing file\n");

    fptr = fopen("sphere.vtk","w");
    fprintf(fptr,"# vtk DataFile Version 2.0\n");
    fprintf(fptr,"Sphere mask\n");
    fprintf(fptr,"ASCII\n");
    fprintf(fptr,"DATASET RECTILINEAR_GRID\n");
    fprintf(fptr,"DIMENSIONS %d %d %d\n", ndimp1, ndimp1, ndimp1);
    
    fprintf(fptr,"X_COORDINATES %d float\n", ndimp1);
    for (int i=0; i < ndimp1; i++){
        coord = i*voxelsize;
        fprintf(fptr, "%.3f ", coord);
    }
    fprintf(fptr,"\n");
    
    fprintf(fptr,"Y_COORDINATES %d float\n", ndimp1);
    for (int j=0; j < ndimp1; j++){
        coord = j*voxelsize;
        fprintf(fptr, "%.3f ", coord);
    }
    fprintf(fptr,"\n");
    
    fprintf(fptr,"Z_COORDINATES %d float\n", ndimp1);
    for (int k=0; k < ndimp1; k++){
        coord = k*voxelsize;
        fprintf(fptr, "%.3f ", coord);
    }
    fprintf(fptr,"\n");
    fprintf(fptr,"CELL_DATA %d\n", ndim*ndim*ndim);
    fprintf(fptr,"FIELD FieldData 1\n");
    fprintf(fptr,"Voxel_value 1 %d float\n", ndim*ndim*ndim);

    printf("writeVTK: writing voxel values\n");
    for (int k=0; k < ndim; k++){
        for (int j=0; j < ndim; j++){
            for (int i=0; i < ndim; i++){
                ii = i + j*ndim + k*ndim*ndim;
                //printf("writeVTK: %d %d\n", ii, image[ii]);
                fprintf(fptr,"%d ", image[ii]);
            }
        }
    }
    fprintf(fptr,"\n");
    fprintf(fptr,"POINT_DATA %d\n",(ndimp1)*(ndimp1)*(ndimp1));

    printf("writeVTK: finish writing file\n");

    fclose(fptr);

}
