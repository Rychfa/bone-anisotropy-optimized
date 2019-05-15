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
#ifndef BONEMAP_ELLIPSOID_H
#define BONEMAP_ELLIPSOID_H

#include <math.h>

///
/// Pre-processor 
///

// #define DEBUG  // TODO move somewhere more centralized?

#ifdef linux
#define M_SQRT3 1.7320508075688772935
#endif

///
/// Public data
///

static int DIRECTIONS[][3] =
	{
//		{ 0,  0, -1},
//		{ 0,  0,  1},
//		{ 0, -1,  0},
//		{ 0, -1, -1},
//		{ 0, -1,  1},
//		{ 0,  1,  0},
//		{ 0,  1, -1},
//		{ 0,  1,  1},
//		{-1,  0,  0},
//		{-1,  0, -1},
//		{-1,  0,  1},
//		{-1, -1,  0},
//		{-1, -1, -1},
//		{-1, -1,  1},
//		{-1,  1,  0},
//		{-1,  1, -1},
//		{-1,  1,  1},
//		{ 1,  0,  0},
//		{ 1,  0, -1},
//		{ 1,  0,  1},
//		{ 1, -1,  0},
//		{ 1, -1, -1},
//		{ 1, -1,  1},
//		{ 1,  1,  0},
//		{ 1,  1, -1},
//		{ 1,  1,  1},

        { 0,  0,  1},
        { 0,  1,  0},
        { 1,  0,  0},

        { 0,  1,  1},
        { 0,  1, -1},
        { 1,  0,  1},
        { 1,  0, -1},
        { 1,  1,  0},
        {-1,  1,  0},

        { 1,  1,  1},
        { 1,  1, -1},
        { 1, -1,  1},
        {-1,  1,  1},
	};

static const double DIRECTIONS_NORMALIZED[][3] = 
	{
//		{ 0,  0, -1},
//		{ 0,  0,  1},
//		{ 0, -1,  0},
//		{ 0, -1/M_SQRT2, -1/M_SQRT2},
//		{ 0, -1/M_SQRT2,  1/M_SQRT2},
//		{ 0,  1,  0},
//		{ 0,  1/M_SQRT2, -1/M_SQRT2},
//		{ 0,  1/M_SQRT2,  1/M_SQRT2},
//		{-1,  0,  0},
//		{-1/M_SQRT2,  0, -1/M_SQRT2},
//		{-1/M_SQRT2,  0,  1/M_SQRT2},
//		{-1/M_SQRT2, -1/M_SQRT2,  0},
//		{-1/M_SQRT3, -1/M_SQRT3, -1/M_SQRT3},
//		{-1/M_SQRT3, -1/M_SQRT3,  1/M_SQRT3},
//		{-1/M_SQRT2,  1/M_SQRT2,  0},
//		{-1/M_SQRT3,  1/M_SQRT3, -1/M_SQRT3},
//		{-1/M_SQRT3,  1/M_SQRT3,  1/M_SQRT3},
//		{ 1,  0,  0},
//		{ 1/M_SQRT2,  0, -1/M_SQRT2},
//		{ 1/M_SQRT2,  0,  1/M_SQRT2},
//		{ 1/M_SQRT2, -1/M_SQRT2,  0},
//		{ 1/M_SQRT3, -1/M_SQRT3, -1/M_SQRT3},
//		{ 1/M_SQRT3, -1/M_SQRT3,  1/M_SQRT3},
//		{ 1/M_SQRT2,  1/M_SQRT2,  0},
//		{ 1/M_SQRT3,  1/M_SQRT3, -1/M_SQRT3},
//		{ 1/M_SQRT3,  1/M_SQRT3,  1/M_SQRT3},

        { 0,          0,          1/M_SQRT2},
        { 0,          1/M_SQRT2,  0},
        { 1/M_SQRT2,  0,          0},

        { 0,          1/M_SQRT2,  1/M_SQRT2},
        { 0,          1/M_SQRT2, -1/M_SQRT2},
        { 1/M_SQRT2,  0,          1/M_SQRT2},
        { 1/M_SQRT2,  0,         -1/M_SQRT2},
        { 1/M_SQRT2,  1/M_SQRT2,  0},
        {-1/M_SQRT2,  1/M_SQRT2,  0},

        { 1/M_SQRT3,  1/M_SQRT3,  1/M_SQRT3},
        { 1/M_SQRT3,  1/M_SQRT3, -1/M_SQRT3},
        { 1/M_SQRT3, -1/M_SQRT3,  1/M_SQRT3},
        {-1/M_SQRT3,  1/M_SQRT3,  1/M_SQRT3},
	};



static const int NUM_DIRECTIONS = sizeof(DIRECTIONS)/sizeof(double)/3;

///
/// Function declarations
///

/*
 * The raw ellipsoid fitting routine
 * @p : pointer to an array of length 3, the points in 3D space on which to fit
 */
void fit_ellipsoid(const double (*p)[3], int n, double Q[3][3]);

/*
 * convenience function, wraps fit_ellipsoid above, given the mils (lengths) along
 * each of the DIRECTIONS defined above
 *
 */
void fit_ellipsoid_mils(const double *mils, double Q[3][3]);

/*
 * debugging tools
 */
#ifdef DEBUG
	void fit_ellipsoid_debug_init(void);
	void fit_ellipsoid_debug_deinit(void);
#endif


#endif //BONEMAP_ELLIPSOID_H
