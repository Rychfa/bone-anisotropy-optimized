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
#ifndef BONEMAP_READER_H
#define BONEMAP_READER_H

#define LOW_RES_FILE "../images/LowRes_F16_R_stance_3p0_segmented/F16_R_stance_3p0_mask.raw"
#define LOW_RES_D1 26
#define LOW_RES_D2 16
#define LOW_RES_D3 31
#define LOW_RES_SIZE (LOW_RES_D1 * LOW_RES_D2 * LOW_RES_D3)
#define HIGH_RES_FILE "../images/HighRes_F16_R_segmented/626_3776_F16_R_H_SEG_cov_crop_rotate.raw"
#define HIGH_RES_D1 1079
#define HIGH_RES_D2 809
#define HIGH_RES_D3 1238
#define HIGH_RES_SIZE (HIGH_RES_D1 * HIGH_RES_D2 * HIGH_RES_D3)

///
/// Function declarations
///
int* readLowResImage ();
int* readHighResImage ();
void destroyImageMatrices ();

#endif //BONEMAP_READER_H
