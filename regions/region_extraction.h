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
#ifndef BONEMAP_REGION_EXTRACTION_H
#define BONEMAP_REGION_EXTRACTION_H

//#define SPHERE_DIAMETER 4
#define HIGH_RES_VOXEL_SIZE 0.082
#define SPHERE_NDIM 80
#define SPHERE_DIAMETER SPHERE_NDIM*0.082
//#define SPHERE_NDIM (int) (SPHERE_DIAMETER /HIGH_RES_VOXEL_SIZE + 1)
#define SPHERE_HALF_NDIM ((int) (SPHERE_NDIM/2))
#define SPHERE_ARRAY_SIZE (SPHERE_NDIM * SPHERE_NDIM * SPHERE_NDIM)
 
//#define DEBUG 

///
/// Function declarations
///
void createSphereMask(double *sphere);
void region_extraction (int i_hr, int j_hr, int k_hr, double *sphere, double *extracted_region, double *ptrHighRes);
void region_extraction_opt1 (int i_hr, int j_hr, int k_hr, double *sphere, double *extracted_region, double *ptrHighRes);
void region_extraction_simd (int i_hr, int j_hr, int k_hr, double *sphere, double *extracted_region, double *ptrHighRes);

#ifdef DEBUG
	void region_extraction_debug_init(void);
	void region_extraction_debug_deinit(void);
#endif

#endif //BONEMAP_REGION_EXTRACTION_H
