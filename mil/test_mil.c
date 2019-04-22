//
// Created by ryan on 22/04/2019.
//

#include "mil.h"
#include "stdio.h"
#include <stdlib.h>

void print_sphere(int ***sphere, int n) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                printf("%d ", sphere[i][j][k]);
            }
            printf("\n");
        }
        printf("****************\n");
    }
}

int ***init_sphere_region_hr(int n) {
    int ***sphere_region_hr = (int ***) malloc(n * sizeof(int **));
    for (int i = 0; i < n; ++i) {
        sphere_region_hr[i] = (int **) malloc(n * sizeof(int *));
        for (int j = 0; j < n; ++j) {
            sphere_region_hr[i][j] = (int *) malloc(n * sizeof(int));
        }
    }

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            sphere_region_hr[i][7][j] = 1;
        }
    }

    return sphere_region_hr;
}

double **build_direction_vectors(int n, int d) {
    double **direction_vector = (double **) malloc(n * sizeof(double *));
    for (int i = 0; i < n; ++i) {
        direction_vector[i] = (double *) malloc(d * sizeof(double));
    }
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < d; ++j) {
            direction_vector[i][j] = 0.0;
        }
    }
    direction_vector[0][1] = 1.0; // DO CROSS THE BONE LINE
    direction_vector[1][2] = -1.0; // DO NOT CROSS THE BONE LINE
    direction_vector[2][1] = 1.0; // DO NOT CROSS THE BONE LINE

    return direction_vector;
}

void destroy_sphere(int ***sphere, int n) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            free(sphere[i][j]);
        }
        free(sphere[i]);
    }
    free(sphere);
}

void destroy_matrix(double **matrix, int n) {
    for (int i = 0; i < n; ++i) {
        free(matrix[i]);
    }
    free(matrix);
}

void print_vector(double *vector, int n) {
    for (int i = 0; i < n; ++i) {
        printf("%f ", vector[i]);
    }
    printf("##\n");
}

void print_matrix(double **matrix, int n, int d) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < d; ++j) {
            printf("%f ", matrix[i][j]);
        }
        printf("\n");
    }
    printf("*******\n");
}


//###################

void test_case_1() {
    int n = 10, n_vectors = 3, dimension = 3;
    int ***sphere_region_hr = init_sphere_region_hr(n);
    double **direction_vectors = build_direction_vectors(n_vectors, dimension);
    print_matrix(direction_vectors, n_vectors, dimension);

    double *mean_intercept_length = mil(sphere_region_hr, n, direction_vectors, n_vectors, dimension);
    print_vector(mean_intercept_length, n_vectors);

    destroy_sphere(sphere_region_hr, n);
    destroy_matrix(direction_vectors, n_vectors);
    free(mean_intercept_length);
}
