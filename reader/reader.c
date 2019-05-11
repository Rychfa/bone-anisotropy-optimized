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
#include <stdio.h>
#include <stdlib.h>
#include "reader.h"

static char* inputFiles[NUMBER_OF_INPUTS] = {
        "../images/LowRes_F16_R_stance_3p0_segmented/F16_R_stance_3p0_mask.raw",
        "../images/LowRes_F16_R_stance_3p0_segmented/scaled2.raw",
        "../images/LowRes_F16_R_stance_3p0_segmented/scaled4.raw",
        "../images/LowRes_F16_R_stance_3p0_segmented/scaled6.raw",
        "../images/LowRes_F16_R_stance_3p0_segmented/scaled8.raw",
        "../images/LowRes_F16_R_stance_3p0_segmented/scaled10.raw"
};

int scaleFactor[NUMBER_OF_INPUTS] = {
        1, //"../images/LowRes_F16_R_stance_3p0_segmented/F16_R_stance_3p0_mask.raw",
        2, //"../images/LowRes_F16_R_stance_3p0_segmented/scaled2.raw",
        4, //"../images/LowRes_F16_R_stance_3p0_segmented/scaled4.raw",
        6, //"../images/LowRes_F16_R_stance_3p0_segmented/scaled6.raw",
        8, //"../images/LowRes_F16_R_stance_3p0_segmented/scaled8.raw",
        10 //"../images/LowRes_F16_R_stance_3p0_segmented/scaled10.raw"
};


typedef union {
    float f;
    char c[4];
} charFloat_t;

typedef union {
    int i;
    char c[4];
} charInt_t;

static int* ptrLowRes;
static int* ptrHighRes;

int* readHighResImage () {

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
    ptrHighRes = malloc (sizeof(int) * HIGH_RES_SIZE);

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

int* readLowResImage (const int index) {

    FILE *fp;
    charFloat_t uCharFloat;
    fp = fopen(inputFiles[index], "r");
    if (fp == NULL) {
        printf("Error while reading low resolution image %d.\n", index);
        return NULL;
    }

    //
    // Allocate memory for low res image
    //
    ptrLowRes = malloc (sizeof(int) * LOW_RES_SIZE(index));

    for (int i = 0; i < LOW_RES_SIZE(index); ++i) {
                //
                // Read bytes and convert to float
                // Bytes are in little endian.
                //
                for (int l = 0; l < 4; l++) {
                    char c = fgetc((FILE*)fp);
                    if (c == EOF) {
                        printf("Error while reading low resolution image %d.\n", index);
                        free(ptrLowRes);
                        return NULL;
                    }
                    uCharFloat.c[l] = c;
                }
                ptrLowRes[i] = (int) uCharFloat.f;
    }

    if (fgetc((FILE*)fp) == EOF) {
        printf("Low res file read correctly!\n");
    }

    return ptrLowRes;
}

void destroyImageMatrices () {
    free(ptrLowRes);
    free(ptrHighRes);
}
