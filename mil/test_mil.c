//
// Created by ryan on 22/04/2019.
//

#include "stdio.h"
#include "stdlib.h"
#include "mil.h"
#include "test_mil.h"
#include "ellipsoid.h"

void print_sphere(int ***sphere, int n) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                char coma_or_space = k < n - 1 ? ',' : ' ';
                printf("%d%c", sphere[i][j][k], coma_or_space);
            }
            printf("\n");
        }
        printf("\n");
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
    direction_vector[2][1] = -1.0; // DO NOT CROSS THE BONE LINE

    return direction_vector;
}

double **build_direction_vectors2(int n, int d) {
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
    direction_vector[2][1] = -1.0; // DO NOT CROSS THE BONE LINE
    direction_vector[3][2] = 1.0; // DO NOT CROSS THE BONE LINE
    direction_vector[4][0] = -1.0; // DO NOT CROSS THE BONE LINE
    direction_vector[5][0] = 1.0; // DO NOT CROSS THE BONE LINE

    direction_vector[6][1] = 1.0; // DO CROSS THE BONE LINE
    direction_vector[6][2] = 1.0;

    direction_vector[7][1] = 1.0; // DO CROSS THE BONE LINE
    direction_vector[7][2] = -1.0;

    direction_vector[8][1] = 1.0; // DO CROSS THE BONE LINE
    direction_vector[8][0] = -1.0;

    return direction_vector;
}

void build_direction_vectors3(double (*directions)[3], int n) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < 3; ++j) {
            directions[i][j] = 0.0;
        }
    }

    directions[0][1] = 1.0; // DO CROSS THE BONE LINE
    directions[1][2] = -1.0; // DO NOT CROSS THE BONE LINE
    directions[2][1] = -1.0; // DO NOT CROSS THE BONE LINE
    directions[3][2] = 1.0; // DO NOT CROSS THE BONE LINE
    directions[4][0] = -1.0; // DO NOT CROSS THE BONE LINE
    directions[5][0] = 1.0; // DO NOT CROSS THE BONE LINE

    directions[6][1] = 1.0; // DO CROSS THE BONE LINE
    directions[6][2] = 1.0;

    directions[7][1] = 1.0; // DO CROSS THE BONE LINE
    directions[7][2] = -1.0;

    directions[8][1] = 1.0; // DO CROSS THE BONE LINE
    directions[8][0] = -1.0;
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
        char coma_or_space = i < n - 1 ? ',' : ' ';
        printf("%f%c", vector[i], coma_or_space);
    }
    printf("\n*******\n");
}

void print_matrix2(double matrix[][3], int n, int d) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < d; ++j) {
            char coma_or_space = j < d - 1 ? ',' : ' ';
            printf("%f%c", matrix[i][j], coma_or_space);
        }
        printf("\n");
    }
    printf("*******\n");
}

void print_matrix3(double matrix[3][3]) {
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            char coma_or_space = j < 3 - 1 ? ',' : ' ';
            printf("%f%c", matrix[i][j], coma_or_space);
        }
        printf("\n");
    }
    printf("*******\n");
}


//###################

void test_mil() {
//    int n = 10, n_vectors = 3, dimension = 3;
//    int ***sphere_region_hr = init_sphere_region_hr(n);
//    double **direction_vectors = build_direction_vectors(n_vectors, dimension);
//    print_matrix2(direction_vectors, n_vectors, dimension);
//
//    double *mean_intercept_length = mil(sphere_region_hr, n, direction_vectors, n_vectors, dimension);
//    print_vector(mean_intercept_length, n_vectors);
//
//    destroy_sphere(sphere_region_hr, n);
//    destroy_matrix(direction_vectors, n_vectors);
//    free(mean_intercept_length);
}

void test_mil_and_ellipsoid() {
    int n = 10, n_vectors = NUM_DIRECTIONS, dimension = 3;
    int ***sphere_region_hr = init_sphere_region_hr(n);
    print_sphere(sphere_region_hr, n);
    print_matrix2(DIRECTIONS, n_vectors, dimension);
//    printf("%f,%f,%f\n", directions[7][0], directions[7][1], directions[7][2]);
//    printf("%f,%f,%f\n", directions[8][0], directions[8][1], directions[8][2]);

    double *mean_intercept_length = mil(sphere_region_hr, n, DIRECTIONS, n_vectors); // seg faulting here
    printf("test\n");
    print_vector(mean_intercept_length, n_vectors);
    double Q[3][3];
    fit_ellipsoid_mils(mean_intercept_length, Q);
    print_matrix3(Q);

    destroy_sphere(sphere_region_hr, n);
    // destroy_matrix(direction_vectors, n_vectors);
    free(mean_intercept_length);
}


void test_mil_and_ellipsoid2()
{

}
