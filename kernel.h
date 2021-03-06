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
#ifndef BONEMAP_KERNEL_H
#define BONEMAP_KERNEL_H

#include <stdbool.h>

#ifdef __cplusplus
extern "C" void init(double** sphere, double** extracted_region, double** ptrHighRes, double** ptrLowRes, double** rotation_matrix, double** ptrEvecOut, double** ptrEvalsOut);
extern "C" void kernel_basic(double* sphere, double* extracted_region, double* ptrHighRes, double* ptrLowRes, double* rotation_matrix, double* ptrEvecOut, double *ptrEvalsOut);
extern "C" void kernel_opt1(double* sphere, double* extracted_region, double* ptrHighRes, double* ptrLowRes, double* rotation_matrix, double* ptrEvecOut, double *ptrEvalsOut);
extern "C" void kernel_opt2(double* sphere, double* extracted_region, double* ptrHighRes, double* ptrLowRes, double* rotation_matrix, double* ptrEvecOut, double *ptrEvalsOut);
extern "C" void deInit(double* sphere, double* extracted_region, double* ptrHighRes, double* ptrLowRes, double* rotation_matrix, double* ptrEvecOut, double *ptrEvalsOut, bool generate_ground_truth);
#else
void init(double** sphere, double** extracted_region, double** ptrHighRes, double** ptrLowRes, double** rotation_matrix, double** ptrEvecOut, double **ptrEvalsOut);
void kernel_basic(double* sphere, double* extracted_region, double* ptrHighRes, double* ptrLowRes, double* rotation_matrix, double* ptrEvecOut, double *ptrEvalsOut);
void kernel_opt1(double* sphere, double* extracted_region, double* ptrHighRes, double* ptrLowRes, double* rotation_matrix, double* ptrEvecOut, double *ptrEvalsOut);
void kernel_opt2(double* sphere, double* extracted_region, double* ptrHighRes, double* ptrLowRes, double* rotation_matrix, double* ptrEvecOut, double *ptrEvalsOut);
void deInit(double* sphere, double* extracted_region, double* ptrHighRes, double* ptrLowRes, double* rotation_matrix, double* ptrEvecOut, double *ptrEvalsOut, bool generate_ground_truth);
#endif

#endif //BONEMAP_KERNEL_H
