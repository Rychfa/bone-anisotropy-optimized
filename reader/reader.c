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
#include "reader.h"

///
/// Skeleton for image reader
///
void reader () {

typedef union {
    float f;
    char c[4];
} charFloat_t;

typedef union {
    int i;
    char c[4];
} charInt_t;

static double* ptrLowRes;
static double* ptrHighRes;

double* readHighResImage () {

    FILE *fp;
    charInt_t uCharInt;
    uCharInt.i = 0;
    fp = fopen(HIGH_RES_FILE, "r");
    if (fp == NULL) {
        printf("Error while reading high resolution image.\n");
        return NULL;
    }

    //
    // Allocate memory for low res image
    //
    ptrHighRes = malloc ((sizeof ptrHighRes ) * HIGH_RES_SIZE);

    for (int i = 0; i < HIGH_RES_SIZE; ++i) {
        //
        // Read bytes and convert to float
        // Bytes are in little endian.
        //
        for (int l = 0; l < 2; l++) {
            char c = fgetc((FILE*)fp);
            if (c == EOF) {
                printf("Error while reading high resolution image.\n");
                free(ptrHighRes);
                return NULL;
            }
            uCharInt.c[l] = c;
        }
        ptrHighRes[i] = uCharInt.i;


    }

    if (fgetc((FILE*)fp) == EOF) {
        printf("High res file read correctly!\n");
    }

    return ptrHighRes;
}

double* readLowResImage () {

    FILE *fp;
    FILE *fd;
    charFloat_t uCharFloat;
    fp = fopen(LOW_RES_FILE, "r");
    if (fp == NULL) {
        printf("Error while reading low resolution image.\n");
        return NULL;
    }

    //
    // Allocate memory for low res image
    //
    ptrLowRes = malloc ((sizeof ptrLowRes ) * LOW_RES_SIZE);

    fd = fopen("debug_reader_lowres.txt","w");
    for (int i = 0; i < LOW_RES_SIZE; ++i) {
                //
                // Read bytes and convert to float
                // Bytes are in little endian.
                //
                for (int l = 0; l < 4; l++) {
                    char c = fgetc((FILE*)fp);
                    if (c == EOF) {
                        printf("Error while reading low resolution image.\n");
                        free(ptrLowRes);
                        return NULL;
                    }
                    uCharFloat.c[l] = c;
                }
                ptrLowRes[i] = uCharFloat.f;
                fprintf(fd,"%d, %.3f\n", i, ptrLowRes[i]);
    }
    fclose(fd);

    if (fgetc((FILE*)fp) == EOF) {
        printf("Low res file read correctly!\n");
    }

    return ptrLowRes;
}

void destroyImageMatrices () {
    free(ptrLowRes);
    free(ptrHighRes);
}
