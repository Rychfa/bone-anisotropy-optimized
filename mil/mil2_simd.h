#ifndef BONEMAP_MIL2_SIMD_H
#define BONEMAP_MIL2_SIMD_H

#include <immintrin.h>
#include <stdio.h>

#define BLOCK_SIZE 16
#define NUM_ACC 4
#define STRIDE 2
/**
 * SIMD IMPLEMNTATION
 */
#define GATHER_SCALE 8
#define ONES _mm256_set1_epi64x(1)
#define threshold _mm256_set1_pd(0.5)


#define SIMD_LOAD_PREV_1D                                                                \
        int offset = STRIDE * n; \
        __m128i bytes_offset_vector1 = _mm_set_epi32(i_prev, offset + i_prev, 2 * offset + i_prev, 3 * offset + i_prev); \
        __m128i bytes_offset_vector2 = _mm_set_epi32(0, STRIDE, 2 * STRIDE, 3 * STRIDE); \
        __m128i bytes_offset_vector3 = _mm_set_epi32(0, STRIDE, 2 * STRIDE, 3 * STRIDE); \
        __m256d tmp_reg1, tmp_reg2, prev_mask_d1, prev_mask_d2;                                                    \
        __m256i tmp_comp1, tmp_comp2; \
        switch (vecID) {                                                                 \
            case 1: /* Vector (1,0,0) */                                                 \
                /*tmp_reg1 = _mm256_i32gather_pd(hr_sphere_region + (knn + j * n), bytes_offset_vector1, GATHER_SCALE); \
                tmp_reg2 = _mm256_i32gather_pd(hr_sphere_region + (knn + (j + STRIDE * NUM_ACC) * n), bytes_offset_vector1, GATHER_SCALE);*/ \
                tmp_reg1 = _mm256_set_pd(hr_sphere_region[ knn + (j+0*STRIDE)*n + i_prev], hr_sphere_region[ knn + (j+1*STRIDE)*n + i_prev], hr_sphere_region[ knn + (j+2*STRIDE)*n + i_prev], hr_sphere_region[ knn + (j+3*STRIDE)*n + i_prev]); \
                tmp_reg2 = _mm256_set_pd(hr_sphere_region[ knn + (j + STRIDE * NUM_ACC) * n + i_prev], hr_sphere_region[ knn + (j + STRIDE* NUM_ACC + STRIDE)*n + i_prev], hr_sphere_region[ knn + (j+STRIDE* NUM_ACC+ 2 * STRIDE)*n + i_prev], hr_sphere_region[ knn + (j+STRIDE* NUM_ACC + 3 * STRIDE)*n + i_prev]);\
                break;                                                                   \
            case 2: /* Vector (0,1,0) */                                                 \
                tmp_reg1 = _mm256_set_pd(hr_sphere_region[ k*n*n + i_prev*n + (j+0*STRIDE)], hr_sphere_region[ k*n*n + i_prev*n + (j+1*STRIDE)], hr_sphere_region[ k*n*n + i_prev*n + (j+2*STRIDE)], hr_sphere_region[ k*n*n + i_prev*n + (j+3*STRIDE)]); \
                tmp_reg2 = _mm256_set_pd(hr_sphere_region[ k*n*n + i_prev*n + (j + STRIDE * NUM_ACC +0*STRIDE)], hr_sphere_region[ k*n*n + i_prev*n + (j + STRIDE * NUM_ACC +1*STRIDE)], hr_sphere_region[ k*n*n + i_prev*n + (j + STRIDE * NUM_ACC +2*STRIDE)], hr_sphere_region[ k*n*n + i_prev*n + (j + STRIDE * NUM_ACC +3*STRIDE)]); \
                break;                                                                   \
            case 3: /* Vector (0,0,1) */                                                 \
                /*tmp_reg1 = _mm256_i32gather_pd(hr_sphere_region + (i_prev*n*n + (k * n) + j), bytes_offset_vector3, GATHER_SCALE); \
                tmp_reg2 = _mm256_i32gather_pd(hr_sphere_region + (i_prev*n*n + k*n + j + STRIDE * NUM_ACC), bytes_offset_vector3, GATHER_SCALE);*/ \
                tmp_reg1 = _mm256_set_pd(hr_sphere_region[ i_prev*n*n + k*n + (j+0*STRIDE)], hr_sphere_region[ i_prev*n*n + k*n + (j+1*STRIDE)], hr_sphere_region[ i_prev*n*n + k*n + (j+2*STRIDE)], hr_sphere_region[ i_prev*n*n + k*n + (j+3*STRIDE)]); \
                tmp_reg2 = _mm256_set_pd(hr_sphere_region[ i_prev*n*n + k*n + (j + STRIDE * NUM_ACC +0*STRIDE)], hr_sphere_region[ i_prev*n*n + k*n + (j + STRIDE * NUM_ACC +1*STRIDE)], hr_sphere_region[ i_prev*n*n + k*n + (j + STRIDE * NUM_ACC +2*STRIDE)], hr_sphere_region[ i_prev*n*n + k*n + (j + STRIDE * NUM_ACC +3*STRIDE)]); \
                break;                                                                   \
            default: ;                                                                   \
        }\
        prev_mask_d1 = _mm256_cmp_pd(tmp_reg1, threshold, _CMP_GT_OQ);\
        prev_mask_d2 = _mm256_cmp_pd(tmp_reg2, threshold, _CMP_GT_OQ);\
        tmp_comp1 = _mm256_castpd_si256(prev_mask_d1); \
        tmp_comp2 = _mm256_castpd_si256(prev_mask_d2); \
        prev_mask = _mm256_and_si256(tmp_comp1, ONES);\
        prev_mask2 = _mm256_and_si256(tmp_comp2, ONES);


#define SIMD_LOAD_DATA_SET_1D                                       \
        __m128i bytes_offset_vector1 = _mm_set_epi32(i, offset + i, 2 * offset + i, 3 * offset + i); \
        __m128i bytes_offset_vector2 = _mm_set_epi32(0, STRIDE, 2 * STRIDE, 3 * STRIDE); \
        __m128i bytes_offset_vector3 = _mm_set_epi32(0, STRIDE, 2 * STRIDE, 3 * STRIDE); \
        switch (vecID) {                                            \
            case 1: /* Vector (1,0,0) */                            \
                /* Unroll over dimension x */                       \
                /*region1 = _mm256_i32gather_pd(hr_sphere_region + (knn + j * n), bytes_offset_vector1, GATHER_SCALE); \
                region2 = _mm256_i32gather_pd(hr_sphere_region + (knn + (j + STRIDE * NUM_ACC) * n), bytes_offset_vector1, GATHER_SCALE);*/ \
                region1 = _mm256_set_pd(hr_sphere_region[ knn + (j+0*STRIDE)*n + i], hr_sphere_region[ knn + (j+1*STRIDE)*n + i], hr_sphere_region[ knn + (j+2*STRIDE)*n + i], hr_sphere_region[ knn + (j+3*STRIDE)*n + i]); \
                region2 = _mm256_set_pd(hr_sphere_region[ knn + (j + STRIDE * NUM_ACC) * n + i], hr_sphere_region[ knn + (j + STRIDE* NUM_ACC + STRIDE)*n + i], hr_sphere_region[ knn + (j+STRIDE* NUM_ACC+ 2 * STRIDE)*n + i], hr_sphere_region[ knn + (j+STRIDE* NUM_ACC + 3 * STRIDE)*n + i]); \
                break;                                              \
            case 2: /* Vector (0,1,0) */                            \
                /* Unroll over dimension y */                       \
                /*region1 = _mm256_i32gather_pd(hr_sphere_region + knn + i * n, bytes_offset_vector2, GATHER_SCALE); \
                region2 = _mm256_i32gather_pd(hr_sphere_region + knn + i * n + j + STRIDE * NUM_ACC, bytes_offset_vector2, GATHER_SCALE);*/ \
                region1 = _mm256_set_pd(hr_sphere_region[ knn + i*n + (j+0*STRIDE)], hr_sphere_region[ knn + i*n + (j+1*STRIDE)], hr_sphere_region[ knn + i*n + (j+2*STRIDE)], hr_sphere_region[ knn + i*n + (j+3*STRIDE)]); \
                region2 = _mm256_set_pd(hr_sphere_region[ knn + i*n + (j + STRIDE * NUM_ACC +0*STRIDE)], hr_sphere_region[ knn + i*n + (j + STRIDE * NUM_ACC +1*STRIDE)], hr_sphere_region[ knn + i*n + (j + STRIDE * NUM_ACC +2*STRIDE)], hr_sphere_region[ knn + i*n + (j + STRIDE * NUM_ACC +3*STRIDE)]); \
                break;                                              \
            case 3: /* Vector (0,0,1) */                            \
                /* Unroll over dimension z */                       \
                /*region1 = _mm256_i32gather_pd(hr_sphere_region + (inn + (k * n) + j), bytes_offset_vector3, GATHER_SCALE); \
                region2 = _mm256_i32gather_pd(hr_sphere_region + (inn + k*n + j + STRIDE * NUM_ACC), bytes_offset_vector3, GATHER_SCALE);*/ \
                region1 = _mm256_set_pd(hr_sphere_region[ inn + k*n + (j+0*STRIDE)], hr_sphere_region[ inn + k*n + (j+1*STRIDE)], hr_sphere_region[ inn + k*n + (j+2*STRIDE)], hr_sphere_region[ inn + k*n + (j+3*STRIDE)]); \
                region2 = _mm256_set_pd(hr_sphere_region[ inn + k*n + (j + STRIDE * NUM_ACC +0*STRIDE)], hr_sphere_region[ inn + k*n + (j + STRIDE * NUM_ACC +1*STRIDE)], hr_sphere_region[ inn + k*n + (j + STRIDE * NUM_ACC +2*STRIDE)], hr_sphere_region[ inn + k*n + (j + STRIDE * NUM_ACC +3*STRIDE)]); \
                break;                                              \
            default: ;                                              \
        }


#define SIMD_LOAD_PREV_2D_POS                                                                          \
    __m256d tmp_reg1, tmp_reg2, prev_mask_d1, prev_mask_d2;                                       \
    __m256i tmp_comp1, tmp_comp2;                                                                 \
    int vector1_id0, vector1_id1, vector1_id2, vector1_id3;                                       \
    int vector2_id0, vector2_id1, vector2_id2, vector2_id3;                                       \
    if (j1 > 0) {                                                                                 \
        i1_prev = i1 - 1;                                                                         \
        j1_prev = j1 - 1;                                                                         \
    }                                                                                             \
    else {                                                                                        \
        i1_prev = i1;                                                                             \
        j1_prev = j1;                                                                             \
    }                                                                                             \
    if (i2 > 0) {                                                                                 \
        i2_prev = i2 - 1;                                                                         \
        j2_prev = j2 - 1;                                                                         \
    }                                                                                             \
    else {                                                                                        \
        i2_prev = i2;                                                                             \
        j2_prev = j2;                                                                             \
    }                                                                                             \
    switch (vecID) {                                                                              \
        case 4: /* Vector (1,1,0) */                                                              \
            /* Unroll over dimensions x,y */                                                      \
            vector1_id0 = k*n*n + j1_prev*n + (i1_prev+1+0*STRIDE);                               \
            vector1_id1 = k*n*n + (j2_prev+1+0*STRIDE)*n + i2_prev;                               \
            vector1_id2 = (k+STRIDE)*n*n + j1_prev*n + (i1_prev+1+0*STRIDE);                      \
            vector1_id3 = (k+STRIDE)*n*n + (j2_prev+1+0*STRIDE)*n + i2_prev;                      \
            \
            vector2_id0 = (k+2*STRIDE)*n*n + j1_prev*n + (i1_prev+1+0*STRIDE);                    \
            vector2_id1 = (k+2*STRIDE)*n*n + (j2_prev+1+0*STRIDE)*n + i2_prev;                    \
            vector2_id2 = (k+3*STRIDE)*n*n + j1_prev*n + (i1_prev+1+0*STRIDE);                    \
            vector2_id3 = (k+3*STRIDE)*n*n + (j2_prev+1+0*STRIDE)*n + i2_prev;                    \
            break;                                                                                \
        case 5: /* Vector (1,0,1) */                                                              \
            /* Unroll over dimensions x,z */                                                      \
            vector1_id0 = j1_prev*n*n + k*n + (i1_prev+1+0*STRIDE);                               \
            vector1_id1 = (j2_prev+1+0*STRIDE)*n*n + k*n + i2_prev;                               \
            vector1_id2 = j1_prev*n*n + (k+STRIDE)*n + (i1_prev+1+0*STRIDE);                      \
            vector1_id3 = (j2_prev+1+0*STRIDE)*n*n + (k+STRIDE)*n + i2_prev;                      \
            \
            vector2_id0 = j1_prev*n*n + (k+2*STRIDE)*n + (i1_prev+1+0*STRIDE);                    \
            vector2_id1 = (j2_prev+1+0*STRIDE)*n*n + (k+2*STRIDE)*n + i2_prev;                    \
            vector2_id2 = j1_prev*n*n + (k+3*STRIDE)*n + (i1_prev+1+0*STRIDE);                    \
            vector2_id3 = (j2_prev+1+0*STRIDE)*n*n + (k+3*STRIDE)*n + i2_prev;                    \
            break;                                                                                \
        case 6: /* Vector (0,1,1) */                                                              \
            /* Unroll over dimensions x,z */                                                      \
            vector1_id0 = (i1_prev+1)*n*n + j1_prev*n + k;                                        \
            vector1_id1 = i2_prev*n*n + (j2_prev+1)*n + k;                                        \
            vector1_id2 = (i1_prev+1)*n*n + j1_prev*n + (k+STRIDE);                               \
            vector1_id3 = i2_prev*n*n + (j2_prev+1)*n + (k+STRIDE);                               \
            \
            vector2_id0 = (i1_prev+1)*n*n + j1_prev*n + (k+2*STRIDE);                             \
            vector2_id1 = i2_prev*n*n + (j2_prev+1)*n + (k+2*STRIDE);                             \
            vector2_id2 = (i1_prev+1)*n*n + j1_prev*n + (k+3*STRIDE);                             \
            vector2_id3 = i2_prev*n*n + (j2_prev+1)*n + (k+3*STRIDE);                             \
            break;                                                                                \
        default: ;                                                                                \
    }                                                                                             \
    double vector1_val0 = hr_sphere_region[vector1_id0];                                          \
    double vector1_val1 = hr_sphere_region[vector1_id1];                                          \
    double vector1_val2 = hr_sphere_region[vector1_id2];                                          \
    double vector1_val3 = hr_sphere_region[vector1_id3];                                          \
    double vector2_val0 = hr_sphere_region[vector2_id0];                                          \
    double vector2_val1 = hr_sphere_region[vector2_id1];                                          \
    double vector2_val2 = hr_sphere_region[vector2_id2];                                          \
    double vector2_val3 = hr_sphere_region[vector2_id3];                                          \
    tmp_reg1 = _mm256_set_pd(vector1_val0, vector1_val1, vector1_val2, vector1_val3);             \
    tmp_reg2 = _mm256_set_pd(vector2_val0, vector2_val1, vector2_val2, vector2_val3);             \
    prev_mask_d1 = _mm256_cmp_pd(tmp_reg1, threshold, _CMP_GT_OQ);                                \
    prev_mask_d2 = _mm256_cmp_pd(tmp_reg2, threshold, _CMP_GT_OQ);                                \
    tmp_comp1 = _mm256_castpd_si256(prev_mask_d1);                                                \
    tmp_comp2 = _mm256_castpd_si256(prev_mask_d2);                                                \
    prev_mask = _mm256_and_si256(tmp_comp1, ONES);                                                \
    prev_mask2 = _mm256_and_si256(tmp_comp2, ONES);

#define SIMD_LOAD_DATA_SET_2D_POS                                                                      \
    switch (vecID) {                                                                              \
        case 4: /* Vector (1,1,0) */                                                              \
            /* Unroll over dimensions x,y */                                                      \
            vector1_id0 = k*n*n + j1*n + (i1+1);                                                  \
            vector1_id1 = k*n*n + (j2+1)*n + i2;                                                  \
            vector1_id2 = (k+STRIDE)*n*n + j1*n + (i1+1);                                         \
            vector1_id3 = (k+STRIDE)*n*n + (j2+1)*n + i2;                                         \
            \
            vector2_id0 = (k+2*STRIDE)*n*n + j1*n + (i1+1);                                       \
            vector2_id1 = (k+2*STRIDE)*n*n + (j2+1)*n + i2;                                       \
            vector2_id2 = (k+3*STRIDE)*n*n + j1*n + (i1+1);                                       \
            vector2_id3 = (k+3*STRIDE)*n*n + (j2+1)*n + i2;                                       \
            break;                                                                                \
        case 5: /* Vector (1,0,1) */                                                              \
            /* Unroll over dimensions x,z */                                                      \
            vector1_id0 = j1*n*n + k*n + (i1+1);                                                  \
            vector1_id1 = (j2+1)*n*n + k*n + i2;                                                  \
            vector1_id2 = j1*n*n + (k+STRIDE)*n + (i1+1);                                         \
            vector1_id3 = (j2+1)*n*n + (k+STRIDE)*n + i2;                                         \
            \
            vector2_id0 = j1*n*n + (k+2*STRIDE)*n + (i1+1);                                       \
            vector2_id1 = (j2+1)*n*n + (k+2*STRIDE)*n + i2;                                       \
            vector2_id2 = j1*n*n + (k+3*STRIDE)*n + (i1+1);                                       \
            vector2_id3 = (j2+1)*n*n + (k+3*STRIDE)*n + i2;                                       \
            break;                                                                                \
        case 6: /* Vector (0,1,1) */                                                              \
            /* Unroll over dimensions y,z */                                                      \
            vector1_id0 = (i1+1)*n*n + j1*n + k;                                                  \
            vector1_id1 = i2*n*n + (j2+1)*n + k;                                                  \
            vector1_id2 = (i1+1)*n*n + j1*n + (k+STRIDE);                                         \
            vector1_id3 = i2*n*n + (j2+1)*n + (k+STRIDE);                                         \
            \
            vector2_id0 = (i1+1)*n*n + j1*n + (k+2*STRIDE);                                       \
            vector2_id1 = i2*n*n + (j2+1)*n + (k+2*STRIDE);                                       \
            vector2_id2 = (i1+1)*n*n + j1*n + (k+3*STRIDE);                                       \
            vector2_id3 = i2*n*n + (j2+1)*n + (k+3*STRIDE);                                       \
            break;                                                                                \
        default: ;                                                                                \
    }                                                                                             \
    double vector1_val0 = hr_sphere_region[vector1_id0];                                          \
    double vector1_val1 = hr_sphere_region[vector1_id1];                                          \
    double vector1_val2 = hr_sphere_region[vector1_id2];                                          \
    double vector1_val3 = hr_sphere_region[vector1_id3];                                          \
    double vector2_val0 = hr_sphere_region[vector2_id0];                                          \
    double vector2_val1 = hr_sphere_region[vector2_id1];                                          \
    double vector2_val2 = hr_sphere_region[vector2_id2];                                          \
    double vector2_val3 = hr_sphere_region[vector2_id3];                                          \
    region1 = _mm256_set_pd(vector1_val0, vector1_val1, vector1_val2, vector1_val3);              \
    region2 = _mm256_set_pd(vector2_val0, vector2_val1, vector2_val2, vector2_val3);              \



#define SIMD_LOAD_PREV_2D_NEG                                                                 \
    __m256d tmp_reg1, tmp_reg2, prev_mask_d1, prev_mask_d2;                              \
    __m256i tmp_comp1, tmp_comp2;                                                        \
    int vector1_id0, vector1_id1, vector1_id2, vector1_id3;                              \
    int vector2_id0, vector2_id1, vector2_id2, vector2_id3;                              \
    if (j1 > 0) {                                                                        \
        i1_prev = i1 + 1;                                                                \
        j1_prev = j1 - 1;                                                                \
    }                                                                                    \
    else {                                                                               \
        i1_prev = i1;                                                                    \
        j1_prev = j1;                                                                    \
    }                                                                                    \
    if (i2 < n-1) {                                                                      \
        i2_prev = i2 + 1;                                                                \
        j2_prev = j2 - 1;                                                                \
    }                                                                                    \
    else {                                                                               \
        i2_prev = i2;                                                                    \
        j2_prev = j2;                                                                    \
    }                                                                                    \
    switch (vecID) {                                                                     \
        case 7: /* Vector (-1,1,0) */                                                    \
            /* Unroll over dimensions x,y */                                             \
            vector1_id0 = k*n*n + j1_prev*n + (i1_prev-1);                               \
            vector1_id1 = k*n*n + (j2_prev+1)*n + i2_prev;                               \
            vector1_id2 = (k+STRIDE)*n*n + j1_prev*n + (i1_prev-1);                      \
            vector1_id3 = (k+STRIDE)*n*n + (j2_prev+1)*n + i2_prev;                      \
            \
            vector2_id0 = (k+2*STRIDE)*n*n + j1_prev*n + (i1_prev-1);                    \
            vector2_id1 = (k+2*STRIDE)*n*n + (j2_prev+1)*n + i2_prev;                    \
            vector2_id2 = (k+3*STRIDE)*n*n + j1_prev*n + (i1_prev-1);                    \
            vector2_id3 = (k+3*STRIDE)*n*n + (j2_prev+1)*n + i2_prev;                    \
            break;                                                                       \
        case 8: /* Vector (-1,0,1) */                                                    \
            /* Unroll over dimensions x,z */                                             \
            vector1_id0 = j1_prev*n*n + k*n + (i1_prev-1);                               \
            vector1_id1 = (j2_prev+1)*n*n + k*n + i2_prev;                               \
            vector1_id2 = j1_prev*n*n + (k+STRIDE)*n + (i1_prev-1);                      \
            vector1_id3 = (j2_prev+1)*n*n + (k+STRIDE)*n + i2_prev;                      \
            \
            vector2_id0 = j1_prev*n*n + (k+2*STRIDE)*n + (i1_prev-1);                    \
            vector2_id1 = (j2_prev+1)*n*n + (k+2*STRIDE)*n + i2_prev;                    \
            vector2_id2 = j1_prev*n*n + (k+3*STRIDE)*n + (i1_prev-1);                    \
            vector2_id3 = (j2_prev+1)*n*n + (k+3*STRIDE)*n + i2_prev;                                                  \
            break;                                                                                \
        case 9: /* Vector (0,1,-1) */                                                             \
            /* Unroll over dimensions x,z */                                                      \
            vector1_id0 = (i1_prev-1)*n*n + j1_prev*n + k;                               \
            vector1_id1 = i2_prev*n*n + (j2_prev+1)*n + k;                               \
            vector1_id2 = (i1_prev-1)*n*n + j1_prev*n + (k+STRIDE);                      \
            vector1_id3 = i2_prev*n*n + (j2_prev+1)*n + (k+STRIDE);                      \
            \
            vector2_id0 = (i1_prev-1)*n*n + j1_prev*n + (k+2*STRIDE);                    \
            vector2_id1 = i2_prev*n*n + (j2_prev+1)*n + (k+2*STRIDE);                    \
            vector2_id2 = (i1_prev-1)*n*n + j1_prev*n + (k+3*STRIDE);                    \
            vector2_id3 = i2_prev*n*n + (j2_prev+1)*n + (k+3*STRIDE);                               \
            break;                                                                                \
        default: ;                                                                                \
    }                                                                                             \
    double vector1_val0 = hr_sphere_region[vector1_id0];                                          \
    double vector1_val1 = hr_sphere_region[vector1_id1];                                          \
    double vector1_val2 = hr_sphere_region[vector1_id2];                                          \
    double vector1_val3 = hr_sphere_region[vector1_id3];                                          \
    double vector2_val0 = hr_sphere_region[vector2_id0];                                          \
    double vector2_val1 = hr_sphere_region[vector2_id1];                                          \
    double vector2_val2 = hr_sphere_region[vector2_id2];                                          \
    double vector2_val3 = hr_sphere_region[vector2_id3];                                          \
    tmp_reg1 = _mm256_set_pd(vector1_val0, vector1_val1, vector1_val2, vector1_val3);             \
    tmp_reg2 = _mm256_set_pd(vector2_val0, vector2_val1, vector2_val2, vector2_val3);             \
    prev_mask_d1 = _mm256_cmp_pd(tmp_reg1, threshold, _CMP_GT_OQ);                                \
    prev_mask_d2 = _mm256_cmp_pd(tmp_reg2, threshold, _CMP_GT_OQ);                                \
    tmp_comp1 = _mm256_castpd_si256(prev_mask_d1);                                                \
    tmp_comp2 = _mm256_castpd_si256(prev_mask_d2);                                                \
    prev_mask = _mm256_and_si256(tmp_comp1, ONES);                                                \
    prev_mask2 = _mm256_and_si256(tmp_comp2, ONES);



#define SIMD_LOAD_DATA_SET_2D_NEG                                                       \
    switch (vecID) {                                                                    \
        case 7: /* Vector (-1,1,0) */                                                   \
            /* Unroll over dimensions x,y */                                            \
            vector1_id0 = k*n*n + j1*n + (i1-1);                                        \
            vector1_id1 = k*n*n + (j2+1)*n + i2;                                        \
            vector1_id2 = (k+STRIDE)*n*n + j1*n + (i1-1);                               \
            vector1_id3 = (k+STRIDE)*n*n + (j2+1)*n + i2;                               \
                                                                                        \
            vector2_id0 = (k+2*STRIDE)*n*n + j1*n + (i1-1);                             \
            vector2_id1 = (k+2*STRIDE)*n*n + (j2+1)*n + i2;                             \
            vector2_id2 = (k+3*STRIDE)*n*n + j1*n + (i1-1);                             \
            vector2_id3 = (k+3*STRIDE)*n*n + (j2+1)*n + i2;                             \
            break;                                                                      \
        case 8: /* Vector (-1,0,1) */                                                   \
            /* Unroll over dimensions x,y */                                            \
            vector1_id0 = j1*n*n + k*n + (i1-1);                                        \
            vector1_id1 = (j2+1)*n*n + k*n + i2;                                        \
            vector1_id2 = j1*n*n + (k+STRIDE)*n + (i1-1);                               \
            vector1_id3 = (j2+1)*n*n + (k+STRIDE)*n + i2;                               \
            \
            vector2_id0 = j1*n*n + (k+2*STRIDE)*n + (i1-1);                             \
            vector2_id1 = (j2+1)*n*n + (k+2*STRIDE)*n + i2;                             \
            vector2_id2 = j1*n*n + (k+3*STRIDE)*n + (i1-1);                             \
            vector2_id3 = (j2+1)*n*n + (k+3*STRIDE)*n + i2;                             \
            break;                                                                      \
        case 9: /* Vector (0,1,-1) */                                                   \
            /* Unroll over dimensions x,y */                                            \
            vector1_id0 = (i1-1)*n*n + j1*n + k;                                        \
            vector1_id1 = i2*n*n + (j2+1)*n + k;                                        \
            vector1_id2 = (i1-1)*n*n + j1*n + (k+STRIDE);                               \
            vector1_id3 = i2*n*n + (j2+1)*n + (k+STRIDE);                               \
            \
            vector2_id0 = (i1-1)*n*n + j1*n + (k+2*STRIDE);                             \
            vector2_id1 = i2*n*n + (j2+1)*n + (k+2*STRIDE);                             \
            vector2_id2 = (i1-1)*n*n + j1*n + (k+3*STRIDE);                             \
            vector2_id3 = i2*n*n + (j2+1)*n + (k+3*STRIDE);                             \
            break;                                                                      \
        default: ;                                                                      \
    }                                                                                   \
    double vector1_val0 = hr_sphere_region[vector1_id0];                                \
    double vector1_val1 = hr_sphere_region[vector1_id1];                                \
    double vector1_val2 = hr_sphere_region[vector1_id2];                                \
    double vector1_val3 = hr_sphere_region[vector1_id3];                                \
    double vector2_val0 = hr_sphere_region[vector2_id0];                                \
    double vector2_val1 = hr_sphere_region[vector2_id1];                                \
    double vector2_val2 = hr_sphere_region[vector2_id2];                                \
    double vector2_val3 = hr_sphere_region[vector2_id3];                                \
    region1 = _mm256_set_pd(vector1_val0, vector1_val1, vector1_val2, vector1_val3);    \
    region2 = _mm256_set_pd(vector2_val0, vector2_val1, vector2_val2, vector2_val3);
// End of LOAD_DATA_SET_2D_NEG
    
#define SIMD_COMPUTATION                                                                      \
    __m256i curr_mask, curr_mask2;                                                            \
    __m256i compi1, compi2;                                                                   \
    bone_count = _mm256_add_pd(bone_count, region1);                                          \
    bone_count2 = _mm256_add_pd(bone_count2, region2);                                        \
    \
    /* Calculate masks */                                                                     \
    compi1 = _mm256_castpd_si256(_mm256_cmp_pd(region1, threshold, _CMP_GT_OQ));              \
    curr_mask = _mm256_and_si256(compi1, ONES);                                               \
    \
    compi2 = _mm256_castpd_si256(_mm256_cmp_pd(region2, threshold, _CMP_GT_OQ));              \
    curr_mask2 = _mm256_and_si256(compi2, ONES);                                              \
    \
    /* Detect edge and add to counter */                                                      \
    __m256i edge1 = _mm256_xor_si256(curr_mask, prev_mask);                                   \
    __m256i edge2 = _mm256_xor_si256(curr_mask2, prev_mask2);                                 \
    edge_count  = _mm256_add_epi64(edge_count, edge1);                                        \
    edge_count2 = _mm256_add_epi64(edge_count2, edge2);

#define BLOCK_KERNEL_1D_SIMD(vec, kk, jj, ii)                                            \
    {                                                                                    \
        const int vecID = vec;                                                           \
        double bone_length_block;                                                        \
        int intercepts_block;                                                            \
        /* Init accumulators */                                                          \
        __m256d bone_count = _mm256_setzero_pd();                                        \
        __m256d bone_count2 = _mm256_setzero_pd();                                       \
                                                                                         \
        __m256i edge_count = _mm256_set1_epi64x(0);                                      \
        __m256i edge_count2 = _mm256_set1_epi64x(0);                                     \
                                                                                         \
        int nn = n * n;                                                                  \
        int i_prev = (ii > 0) ? ii - 1 : 0;                                              \
                                                                                         \
        for (int k = kk + 1; k < kk + BLOCK_SIZE; k += STRIDE) {                         \
            int knn = k * nn;                                                            \
            for (int j = jj + 1; j < jj + BLOCK_SIZE; j += STRIDE*NUM_ACC*2) {           \
                __m256i prev_mask, prev_mask2;                                           \
                                                                                         \
                SIMD_LOAD_PREV_1D                                                        \
                                                                                         \
                for (int i = ii; i < ii + BLOCK_SIZE; ++i) {                             \
                    int inn = i * nn;                                                    \
                    __m256d region1;                                                     \
                    __m256d region2;                                                     \
                                                                                         \
                    /* Load working set */                                               \
                    SIMD_LOAD_DATA_SET_1D                                                \
                                                                                         \
                    /* Perform computation */                                            \
                    SIMD_COMPUTATION                                                     \
                                                                                         \
                    /* Update state of prev_mask */                                      \
                    prev_mask = curr_mask;                                               \
                    prev_mask2 = curr_mask2;                                             \
                }                                                                        \
            }                                                                            \
        }                                                                                \
        bone_count = _mm256_add_pd(bone_count, bone_count2);                             \
        edge_count = _mm256_add_epi64(edge_count, edge_count2);                          \
                                                                                         \
        bone_length_block = horizontal_add(bone_count);                                  \
        intercepts_block  = horizontal_addi(edge_count);                                 \
                                                                                         \
        bone_length[vecID-1] += bone_length_block;                                       \
        intercepts[vecID-1] += intercepts_block;                                         \
    }

#define BLOCK_KERNEL_2D_POS_SIMD(vec, kk, jj, ii)                                        \
    {                                                                                    \
        const int vecID = vec;                                                           \
        double bone_length_block;                                                        \
        int intercepts_block;                                                            \
        /* Init accumulators */                                                          \
        __m256d bone_count = _mm256_setzero_pd();                                        \
        __m256d bone_count2 = _mm256_setzero_pd();                                       \
                                                                                         \
        __m256i edge_count = _mm256_set1_epi64x(0);                                      \
        __m256i edge_count2 = _mm256_set1_epi64x(0);                                     \
                                                                                         \
                                                                                         \
        for (int k = kk + 1; k < kk + BLOCK_SIZE; k += STRIDE * NUM_ACC) {               \
            for (int ij = 0; ij < BLOCK_SIZE; ij += STRIDE) {                            \
                unsigned int i1_prev, i2_prev, j1_prev, j2_prev;                         \
                __m256i prev_mask, prev_mask2;                                           \
                                                                                         \
                int i1 = ii + ij;                                                        \
                int j1 = jj;                                                             \
                int i2 = ii;                                                             \
                int j2 = jj + ij;                                                        \
                                                                                         \
                /* Initialise previous mask */                                           \
                SIMD_LOAD_PREV_2D_POS                                                         \
                                                                                         \
                while (i1 + 1 < ii + BLOCK_SIZE && j2 + 1 < jj + BLOCK_SIZE) {           \
                    __m256d region1;                                                     \
                    __m256d region2;                                                     \
                                                                                         \
                    /* Load working set */                                               \
                    SIMD_LOAD_DATA_SET_2D_POS                                                 \
                                                                                         \
                    /* Perform computation */                                            \
                    SIMD_COMPUTATION                                                     \
                                                                                         \
                    /* Update state of prev_mask */                                      \
                    prev_mask = curr_mask;                                               \
                    prev_mask2 = curr_mask2;                                             \
                                                                                         \
                    ++i1;                                                                \
                    ++j1;                                                                \
                    ++i2;                                                                \
                    ++j2;                                                                \
                }                                                                        \
            }                                                                            \
        } /* End iteration over dimension k */                                           \
                                                                                         \
        bone_count = _mm256_add_pd(bone_count, bone_count2);                             \
        edge_count = _mm256_add_epi64(edge_count, edge_count2);                          \
                                                                                         \
        bone_length_block  = horizontal_add(bone_count);                                 \
        intercepts_block  = horizontal_addi(edge_count);                                 \
                                                                                         \
        bone_length[vecID-1] += bone_length_block;                                       \
        intercepts[vecID-1] += intercepts_block;                                         \
    }

#define BLOCK_KERNEL_2D_NEG_SIMD(vec, kk, jj, ii)                                        \
    {                                                                                    \
        const int vecID = vec;                                                           \
        double bone_length_block = 0;                                                    \
        int intercepts_block = 0;                                                        \
        /* Init accumulators */                                                          \
        __m256d bone_count = _mm256_setzero_pd();                                        \
        __m256d bone_count2 = _mm256_setzero_pd();                                       \
                                                                                         \
        __m256i edge_count = _mm256_set1_epi64x(0);                                      \
        __m256i edge_count2 = _mm256_set1_epi64x(0);                                     \
        for (int k = kk + 1; k < kk + BLOCK_SIZE; k += STRIDE * NUM_ACC) {               \
            for (int ij = 0; ij < BLOCK_SIZE; ij += STRIDE) {                            \
                unsigned int i1_prev, i2_prev, j1_prev, j2_prev;                         \
                __m256i prev_mask, prev_mask2;                                           \
                /* Start at the end of the row in the block */                           \
                int i1 = ii + (BLOCK_SIZE-1) - ij;                                       \
                int j1 = jj;                                                             \
                int i2 = ii + (BLOCK_SIZE-1);                                            \
                int j2 = jj + ij;                                                        \
                                                                                         \
                /* Initialise previous mask */                                           \
                SIMD_LOAD_PREV_2D_NEG                                                    \
                                                                                         \
                while (j2 + 1 < jj + BLOCK_SIZE) {                                       \
                    __m256d region1;                                                     \
                    __m256d region2;                                                     \
                                                                                         \
                    /* Load working set */                                               \
                    SIMD_LOAD_DATA_SET_2D_NEG                                            \
                                                                                         \
                    /* Perform computation */                                            \
                    SIMD_COMPUTATION                                                     \
                                                                                         \
                    /* Update state of prev_mask */                                      \
                    prev_mask = curr_mask;                                               \
                    prev_mask2 = curr_mask2;                                             \
                                                                                         \
                    --i1;                                                                \
                    ++j1;                                                                \
                    --i2;                                                                \
                    ++j2;                                                                \
                }                                                                        \
            }                                                                            \
        } /* End iteration over dimension k */                                           \
                                                                                         \
        bone_count = _mm256_add_pd(bone_count, bone_count2);                             \
        edge_count = _mm256_add_epi64(edge_count, edge_count2);                          \
                                                                                         \
        bone_length_block = horizontal_add(bone_count);                                  \
        intercepts_block  = horizontal_addi(edge_count);                                 \
                                                                                         \
        bone_length[vecID-1] += bone_length_block;                                       \
        intercepts[vecID-1] += intercepts_block;                                         \
    }


#define LOAD_PREV_3D_SIMD(vecID, kb, jb, ib)                                                      \
    __m256d tmp_reg1, tmp_reg2, tmp_reg3, prev_mask_d1, prev_mask_d2, prev_mask_d3;               \
    __m256i tmp_comp1, tmp_comp2, tmp_comp3;                                                      \
    int vector1_id0, vector1_id1, vector1_id2, vector1_id3;                                       \
    int vector2_id0, vector2_id1, vector2_id2, vector2_id3;                                       \
    int vector3_id0, vector3_id1, vector3_id2, vector3_id3;                                       \
    switch (vecID) {                                                                              \
        case 10: /* Vector (1,1,1) */                                                             \
            prev1 = (ib > 0) ? 1 : 0;                                                                                            \
            prev2 = (ib > 0 && jb+j2 > 0) ? 1 : 0;                                                                               \
            vector1_id0 = (kb+k1-prev1)*n*n        + (jb+j1-prev1)*n        + (ib+i-prev1);                                      \
            vector1_id1 = (kb+k1+STRIDE-prev1)*n*n + (jb+j1-prev1)*n        + (ib+i-prev1);                                      \
            vector1_id2 = (kb+k2-prev2)*n*n        + (jb+j2-prev2)*n        + (ib+i-prev2);                                      \
            vector1_id3 = (kb+k2-prev1)*n*n        + (jb+j2+STRIDE-prev1)*n + (ib+i-prev1);                                      \
                                                                                                                                 \
            prev3 = (jb > 0) ? 1 : 0;                                                                                            \
            prev4 = (jb > 0 && kb+j2 > 0) ? 1 : 0;                                                                               \
            vector2_id0 = (kb+j1-prev3)*n*n        + (jb+i-prev3)*n         + (ib+k1)-prev3;                                     \
            vector2_id1 = (kb+j1-prev3)*n*n        + (jb+i-prev3)*n         + (ib+k1+STRIDE)-prev3;                              \
            vector2_id2 = (kb+j2-prev4)*n*n        + (jb+i-prev4)*n         + (ib+k2)-prev4;                                     \
            vector2_id3 = (kb+j2+STRIDE-prev3)*n*n + (jb+i-prev3)*n         + (ib+k2)-prev3;                                     \
                                                                                                                                 \
            prev5 = (kb > 0) ? 1 : 0;                                                                                            \
            prev6 = (kb > 0 && ib+j2 > 0) ? 1 : 0;                                                                               \
            vector3_id0 = (kb+i-prev5)*n*n        + (jb+k1-prev5)*n         + (ib+j1-prev5);                                     \
            vector3_id1 = (kb+i-prev5)*n*n        + (jb+k1+STRIDE-prev5)*n  + (ib+j1-prev5);                                     \
            vector3_id2 = (kb+i-prev6)*n*n        + (jb+k2-prev6)*n         + (ib+j2-prev6);                                     \
            vector3_id3 = (kb+i-prev5)*n*n        + (jb+k2-prev5)*n         + (ib+j2+STRIDE-prev5);                              \
            break;                                                                                                               \
        case 11: /* Vector (1,1,-1) */                                                                                           \
            prev1 = (ib > 0) ? 1 : 0;                                                                                            \
            prev2 = (ib > 0 && jb+j2 > 0) ? 1 : 0;                                                                               \
            vector1_id0 = (kb+BLOCK_SIZE-1-k1+prev1)*n*n        + (jb+j1-prev1)*n        + (ib+i-prev1);                         \
            vector1_id1 = (kb+BLOCK_SIZE-1-k1-STRIDE+prev1)*n*n + (jb+j1-prev1)*n        + (ib+i-prev1);                         \
            vector1_id2 = (kb+BLOCK_SIZE-1-k2+prev2)*n*n        + (jb+j2-prev2)*n        + (ib+i-prev2);                         \
            vector1_id3 = (kb+BLOCK_SIZE-1-k2+prev1)*n*n        + (jb+j2+STRIDE-prev1)*n + (ib+i-prev1);                         \
                                                                                                                                 \
            prev3 = (jb > 0) ? 1 : 0;                                                                                            \
            prev4 = (jb > 0 && kb+BLOCK_SIZE-1-j2 < n-1) ? 1 : 0;                                                                \
            vector2_id0 = (kb+BLOCK_SIZE-1-j1+prev3)*n*n        + (jb+i-prev3)*n         + (ib+k1)-prev3;                        \
            vector2_id1 = (kb+BLOCK_SIZE-1-j1+prev3)*n*n        + (jb+i-prev3)*n         + (ib+k1+STRIDE)-prev3;                 \
            vector2_id2 = (kb+BLOCK_SIZE-1-j2+prev4)*n*n        + (jb+i-prev4)*n         + (ib+k2)-prev4;                        \
            vector2_id3 = (kb+BLOCK_SIZE-1-j2-STRIDE+prev3)*n*n + (jb+i-prev3)*n         + (ib+k2)-prev3;                        \
                                                                                                                                 \
            prev5 = (kb+BLOCK_SIZE-1-i < n-1) ? 1 : 0;                                                                           \
            prev6 = (kb+BLOCK_SIZE-1-i < n-1 && ib+j2 > 0) ? 1 : 0;                                                              \
            vector3_id0 = (kb+BLOCK_SIZE-1-i+prev5)*n*n        + (jb+k1-prev5)*n         + (ib+j1-prev5);                        \
            vector3_id1 = (kb+BLOCK_SIZE-1-i+prev5)*n*n        + (jb+k1+STRIDE-prev5)*n  + (ib+j1-prev5);                        \
            vector3_id2 = (kb+BLOCK_SIZE-1-i+prev6)*n*n        + (jb+k2-prev6)*n         + (ib+j2-prev6);                        \
            vector3_id3 = (kb+BLOCK_SIZE-1-i+prev5)*n*n        + (jb+k2-prev5)*n         + (ib+j2+STRIDE-prev5);                 \
            break;                                                                                                               \
        case 12: /* Vector (1,-1,1) */                                                                                           \
            prev1 = (ib > 0) ? 1 : 0;                                                                                            \
            prev2 = (ib > 0 && jb+BLOCK_SIZE-j2-1 < n-1) ? 1 : 0;                                                                \
            vector1_id0 = (kb+k1-prev1)*n*n        + (jb+BLOCK_SIZE-j1-1+prev1)*n        + (ib+i-prev1);                         \
            vector1_id1 = (kb+k1+STRIDE-prev1)*n*n + (jb+BLOCK_SIZE-j1-1+prev1)*n        + (ib+i-prev1);                         \
            vector1_id2 = (kb+k2-prev2)*n*n        + (jb+BLOCK_SIZE-j2-1+prev2)*n        + (ib+i-prev2);                         \
            vector1_id3 = (kb+k2-prev1)*n*n        + (jb+BLOCK_SIZE-j2-1-STRIDE+prev1)*n + (ib+i-prev1);                         \
                                                                                                                                 \
            prev3 = (jb+BLOCK_SIZE-i-1 < n-1) ? 1 : 0;                                                                           \
            prev4 = (jb+BLOCK_SIZE-i-1 < n-1 && kb+j2 > 0) ? 1 : 0;                                                              \
            vector2_id0 = (kb+j1-prev3)*n*n        + (jb+BLOCK_SIZE-i-1+prev3)*n         + (ib+k1)-prev3;                        \
            vector2_id1 = (kb+j1-prev3)*n*n        + (jb+BLOCK_SIZE-i-1+prev3)*n         + (ib+k1+STRIDE)-prev3;                 \
            vector2_id2 = (kb+j2-prev4)*n*n        + (jb+BLOCK_SIZE-i-1+prev4)*n         + (ib+k2)-prev4;                        \
            vector2_id3 = (kb+j2+STRIDE-prev3)*n*n + (jb+BLOCK_SIZE-i-1+prev3)*n         + (ib+k2)-prev3;                        \
                                                                                                                                 \
            prev5 = (kb > 0) ? 1 : 0;                                                                                            \
            prev6 = (kb > 0 && ib+j2 > 0) ? 1 : 0;                                                                               \
            vector3_id0 = (kb+i-prev5)*n*n        + (jb+BLOCK_SIZE-k1-1+prev5)*n         + (ib+j1-prev5);                        \
            vector3_id1 = (kb+i-prev5)*n*n        + (jb+BLOCK_SIZE-k1-1-STRIDE+prev5)*n  + (ib+j1-prev5);                        \
            vector3_id2 = (kb+i-prev6)*n*n        + (jb+BLOCK_SIZE-k2-1+prev6)*n         + (ib+j2-prev6);                        \
            vector3_id3 = (kb+i-prev5)*n*n        + (jb+BLOCK_SIZE-k2-1+prev5)*n         + (ib+j2+STRIDE-prev5);                 \
            break;                                                                                                               \
        case 13: /* Vector (-1,1,1) */                                                                                           \
            prev1 = (ib+BLOCK_SIZE-1-i < n-1) ? 1 : 0;                                                                           \
            prev2 = (ib+BLOCK_SIZE-1-i < n-1 && jb+j2 > 0) ? 1 : 0;                                                              \
            vector1_id0 = (kb+k1-prev1)*n*n        + (jb+j1-prev1)*n        + (ib+BLOCK_SIZE-1-i+prev1);                         \
            vector1_id1 = (kb+k1+STRIDE-prev1)*n*n + (jb+j1-prev1)*n        + (ib+BLOCK_SIZE-1-i+prev1);                         \
            vector1_id2 = (kb+k2-prev2)*n*n        + (jb+j2-prev2)*n        + (ib+BLOCK_SIZE-1-i+prev2);                         \
            vector1_id3 = (kb+k2-prev1)*n*n        + (jb+j2+STRIDE-prev1)*n + (ib+BLOCK_SIZE-1-i+prev1);                         \
                                                                                                                                 \
            prev3 = (jb > 0) ? 1 : 0;                                                                                            \
            prev4 = (jb > 0 && kb+j2 > 0) ? 1 : 0;                                                                               \
            vector2_id0 = (kb+j1-prev3)*n*n        + (jb+i-prev3)*n         + (ib+BLOCK_SIZE-1-k1+prev3);                        \
            vector2_id1 = (kb+j1-prev3)*n*n        + (jb+i-prev3)*n         + (ib+BLOCK_SIZE-1-k1-STRIDE+prev3);                 \
            vector2_id2 = (kb+j2-prev4)*n*n        + (jb+i-prev4)*n         + (ib+BLOCK_SIZE-1-k2+prev4);                        \
            vector2_id3 = (kb+j2+STRIDE-prev3)*n*n + (jb+i-prev3)*n         + (ib+BLOCK_SIZE-1-k2+prev3);                        \
                                                                                                                                 \
            prev5 = (kb > 0) ? 1 : 0;                                                                                            \
            prev6 = (kb > 0 && ib+BLOCK_SIZE-1-j2 < n-1) ? 1 : 0;                                                                \
            vector3_id0 = (kb+i-prev5)*n*n        + (jb+k1-prev5)*n         + (ib+BLOCK_SIZE-1-j1+prev5);                        \
            vector3_id1 = (kb+i-prev5)*n*n        + (jb+k1+STRIDE-prev5)*n  + (ib+BLOCK_SIZE-1-j1+prev5);                        \
            vector3_id2 = (kb+i-prev6)*n*n        + (jb+k2-prev6)*n         + (ib+BLOCK_SIZE-1-j2+prev6);                        \
            vector3_id3 = (kb+i-prev5)*n*n        + (jb+k2-prev5)*n         + (ib+BLOCK_SIZE-1-j2-STRIDE+prev5);                 \
            break;                                                                                                               \
        default: ;                                                                                                               \
    }                                                                                             \
    double vector1_val0 = hr_sphere_region[vector1_id0];                                          \
    double vector1_val1 = hr_sphere_region[vector1_id1];                                          \
    double vector1_val2 = hr_sphere_region[vector1_id2];                                          \
    double vector1_val3 = hr_sphere_region[vector1_id3];                                          \
    double vector2_val0 = hr_sphere_region[vector2_id0];                                          \
    double vector2_val1 = hr_sphere_region[vector2_id1];                                          \
    double vector2_val2 = hr_sphere_region[vector2_id2];                                          \
    double vector2_val3 = hr_sphere_region[vector2_id3];                                          \
    double vector3_val0 = hr_sphere_region[vector3_id0];                                          \
    double vector3_val1 = hr_sphere_region[vector3_id1];                                          \
    double vector3_val2 = hr_sphere_region[vector3_id2];                                          \
    double vector3_val3 = hr_sphere_region[vector3_id3];                                          \
    tmp_reg1 = _mm256_set_pd(vector1_val0, vector1_val1, vector1_val2, vector1_val3);             \
    tmp_reg2 = _mm256_set_pd(vector2_val0, vector2_val1, vector2_val2, vector2_val3);             \
    tmp_reg3 = _mm256_set_pd(vector3_val0, vector3_val1, vector3_val2, vector3_val3);             \
    prev_mask_d1 = _mm256_cmp_pd(tmp_reg1, threshold, _CMP_GT_OQ);                                \
    prev_mask_d2 = _mm256_cmp_pd(tmp_reg2, threshold, _CMP_GT_OQ);                                \
    prev_mask_d3 = _mm256_cmp_pd(tmp_reg3, threshold, _CMP_GT_OQ);                                \
    tmp_comp1 = _mm256_castpd_si256(prev_mask_d1);                                                \
    tmp_comp2 = _mm256_castpd_si256(prev_mask_d2);                                                \
    tmp_comp3 = _mm256_castpd_si256(prev_mask_d3);                                                \
    prev_mask1 = _mm256_and_si256(tmp_comp1, ONES);                                               \
    prev_mask2 = _mm256_and_si256(tmp_comp2, ONES);                                               \
    prev_mask3 = _mm256_and_si256(tmp_comp3, ONES);

#define LOAD_DATA_3D_SIMD(vecID, kb, jb, ib)                                                       \
    switch (vecID) {                                                                               \
        case 10: /* Vector (1,1,1) */                                                              \
            vector1_id0 = (kb+k1)*n*n        + (jb+j1)*n        + (ib+i);                          \
            vector1_id1 = (kb+k1+STRIDE)*n*n + (jb+j1)*n        + (ib+i);                          \
            vector1_id2 = (kb+k2)*n*n        + (jb+j2)*n        + (ib+i);                          \
            vector1_id3 = (kb+k2)*n*n        + (jb+j2+STRIDE)*n + (ib+i);                          \
                                                                                                   \
            vector2_id0 = (kb+j1)*n*n        + (jb+i)*n         + (ib+k1);                         \
            vector2_id1 = (kb+j1)*n*n        + (jb+i)*n         + (ib+k1+STRIDE);                  \
            vector2_id2 = (kb+j2)*n*n        + (jb+i)*n         + (ib+k2);                         \
            vector2_id3 = (kb+j2+STRIDE)*n*n + (jb+i)*n         + (ib+k2);                         \
                                                                                                   \
            vector3_id0 = (kb+i)*n*n        + (jb+k1)*n         + (ib+j1);                         \
            vector3_id1 = (kb+i)*n*n        + (jb+k1+STRIDE)*n  + (ib+j1);                         \
            vector3_id2 = (kb+i)*n*n        + (jb+k2)*n         + (ib+j2);                         \
            vector3_id3 = (kb+i)*n*n        + (jb+k2)*n         + (ib+j2+STRIDE);                  \
            break;                                                                                 \
        case 11: /* Vector (1,1,-1) */                                                             \
            vector1_id0 = (kb+BLOCK_SIZE-1-k1)*n*n        + (jb+j1)*n        + (ib+i);                            \
            vector1_id1 = (kb+BLOCK_SIZE-1-k1-STRIDE)*n*n + (jb+j1)*n        + (ib+i);                            \
            vector1_id2 = (kb+BLOCK_SIZE-1-k2)*n*n        + (jb+j2)*n        + (ib+i);                            \
            vector1_id3 = (kb+BLOCK_SIZE-1-k2)*n*n        + (jb+j2+STRIDE)*n + (ib+i);                            \
                                                                                                                  \
            vector2_id0 = (kb+BLOCK_SIZE-1-j1)*n*n        + (jb+i)*n         + (ib+k1);                           \
            vector2_id1 = (kb+BLOCK_SIZE-1-j1)*n*n        + (jb+i)*n         + (ib+k1+STRIDE);                    \
            vector2_id2 = (kb+BLOCK_SIZE-1-j2)*n*n        + (jb+i)*n         + (ib+k2);                           \
            vector2_id3 = (kb+BLOCK_SIZE-1-j2-STRIDE)*n*n + (jb+i)*n         + (ib+k2);                           \
                                                                                                                  \
            vector3_id0 = (kb+BLOCK_SIZE-1-i)*n*n        + (jb+k1)*n         + (ib+j1);                           \
            vector3_id1 = (kb+BLOCK_SIZE-1-i)*n*n        + (jb+k1+STRIDE)*n  + (ib+j1);                           \
            vector3_id2 = (kb+BLOCK_SIZE-1-i)*n*n        + (jb+k2)*n         + (ib+j2);                           \
            vector3_id3 = (kb+BLOCK_SIZE-1-i)*n*n        + (jb+k2)*n         + (ib+j2+STRIDE);                    \
            break;                                                                                                \
        case 12: /* Vector (1,-1,1) */                                                                            \
            vector1_id0 = (kb+k1)*n*n        + (jb+BLOCK_SIZE-j1-1)*n        + (ib+i);                            \
            vector1_id1 = (kb+k1+STRIDE)*n*n + (jb+BLOCK_SIZE-j1-1)*n        + (ib+i);                            \
            vector1_id2 = (kb+k2)*n*n        + (jb+BLOCK_SIZE-j2-1)*n        + (ib+i);                            \
            vector1_id3 = (kb+k2)*n*n        + (jb+BLOCK_SIZE-j2-1-STRIDE)*n + (ib+i);                            \
                                                                                                                  \
            vector2_id0 = (kb+j1)*n*n        + (jb+BLOCK_SIZE-i-1)*n         + (ib+k1);                           \
            vector2_id1 = (kb+j1)*n*n        + (jb+BLOCK_SIZE-i-1)*n         + (ib+k1+STRIDE);                    \
            vector2_id2 = (kb+j2)*n*n        + (jb+BLOCK_SIZE-i-1)*n         + (ib+k2);                           \
            vector2_id3 = (kb+j2+STRIDE)*n*n + (jb+BLOCK_SIZE-i-1)*n         + (ib+k2);                          \
                                                                                                                  \
            vector3_id0 = (kb+i)*n*n        + (jb+BLOCK_SIZE-k1-1)*n         + (ib+j1);                           \
            vector3_id1 = (kb+i)*n*n        + (jb+BLOCK_SIZE-k1-1-STRIDE)*n  + (ib+j1);                           \
            vector3_id2 = (kb+i)*n*n        + (jb+BLOCK_SIZE-k2-1)*n         + (ib+j2);                           \
            vector3_id3 = (kb+i)*n*n        + (jb+BLOCK_SIZE-k2-1)*n         + (ib+j2+STRIDE);                    \
    break;                                                                                                        \
        case 13: /* Vector (-1,1,1) */                                                                            \
            vector1_id0 = (kb+k1)*n*n        + (jb+j1)*n        + (ib+BLOCK_SIZE-1-i);                            \
            vector1_id1 = (kb+k1+STRIDE)*n*n + (jb+j1)*n        + (ib+BLOCK_SIZE-1-i);                            \
            vector1_id2 = (kb+k2)*n*n        + (jb+j2)*n        + (ib+BLOCK_SIZE-1-i);                            \
            vector1_id3 = (kb+k2)*n*n        + (jb+j2+STRIDE)*n + (ib+BLOCK_SIZE-1-i);                            \
                                                                                                                  \
            vector2_id0 = (kb+j1)*n*n        + (jb+i)*n         + (ib+BLOCK_SIZE-1-k1);                           \
            vector2_id1 = (kb+j1)*n*n        + (jb+i)*n         + (ib+BLOCK_SIZE-1-k1-STRIDE);                    \
            vector2_id2 = (kb+j2)*n*n        + (jb+i)*n         + (ib+BLOCK_SIZE-1-k2);                           \
            vector2_id3 = (kb+j2+STRIDE)*n*n + (jb+i)*n         + (ib+BLOCK_SIZE-1-k2);                           \
                                                                                                                  \
            vector3_id0 = (kb+i)*n*n        + (jb+k1)*n         + (ib+BLOCK_SIZE-1-j1);                           \
            vector3_id1 = (kb+i)*n*n        + (jb+k1+STRIDE)*n  + (ib+BLOCK_SIZE-1-j1);                           \
            vector3_id2 = (kb+i)*n*n        + (jb+k2)*n         + (ib+BLOCK_SIZE-1-j2);                           \
            vector3_id3 = (kb+i)*n*n        + (jb+k2)*n         + (ib+BLOCK_SIZE-1-j2-STRIDE);                    \
            break;                                                                                \
        default: ;                                                                                \
    }                                                                                             \
    double vector1_val0 = hr_sphere_region[vector1_id0];                                          \
    double vector1_val1 = hr_sphere_region[vector1_id1];                                          \
    double vector1_val2 = hr_sphere_region[vector1_id2];                                          \
    double vector1_val3 = hr_sphere_region[vector1_id3];                                          \
    double vector2_val0 = hr_sphere_region[vector2_id0];                                          \
    double vector2_val1 = hr_sphere_region[vector2_id1];                                          \
    double vector2_val2 = hr_sphere_region[vector2_id2];                                          \
    double vector2_val3 = hr_sphere_region[vector2_id3];                                          \
    double vector3_val0 = hr_sphere_region[vector3_id0];                                          \
    double vector3_val1 = hr_sphere_region[vector3_id1];                                          \
    double vector3_val2 = hr_sphere_region[vector3_id2];                                          \
    double vector3_val3 = hr_sphere_region[vector3_id3];                                          \
    region1 = _mm256_set_pd(vector1_val0, vector1_val1, vector1_val2, vector1_val3);              \
    region2 = _mm256_set_pd(vector2_val0, vector2_val1, vector2_val2, vector2_val3);              \
    region3 = _mm256_set_pd(vector3_val0, vector3_val1, vector3_val2, vector3_val3);              \

#define LOAD_PREV_REMAINDER_3D_SIMD(vecID, kb, jb, ib)                                            \
    __m256d tmp_reg1, tmp_reg2, prev_mask_d1, prev_mask_d2;                                       \
    __m256i tmp_comp1, tmp_comp2;                                                                 \
    int vector1_id0, vector1_id1, vector1_id2;                                                    \
    int vector2_id0, vector2_id1, vector2_id2;                                                    \
    switch (vecID) {                                                                                                       \
        case 10: /* Vector (1,1,1) */                                                                                      \
            prev1 = (ib > 0) ? 1 : 0;                                                                                      \
            prev2 = (ib > 0 && jb+j > STRIDE) ? 1 : 0;                                                                     \
            vector1_id0 = (kb+k-prev1)*n*n        + (jb+j-prev1)*n        + (ib+i-prev1);                                  \
            vector1_id1 = (kb+k-prev2)*n*n        + (jb+j-STRIDE-prev2)*n + (ib+i-prev2);                                  \
                                                                                                                           \
            prev3 = (jb > 0) ? 1 : 0;                                                                                      \
            prev4 = (jb > 0 && kb+k > STRIDE) ? 1 : 0;                                                                     \
            vector1_id2 = (kb+k-prev3)*n*n        + (jb+i-prev3)*n        + (ib+j-prev3);                                  \
            vector2_id0 = (kb+k-STRIDE-prev4)*n*n + (jb+i-prev4)*n        + (ib+j-prev4);                                  \
                                                                                                                           \
            prev5 = (kb > 0) ? 1 : 0;                                                                                      \
            prev6 = (kb > 0 && ib+j > STRIDE) ? 1 : 0;                                                                     \
            vector2_id1 = (kb+i-prev5)*n*n       + (jb+k-prev5)*n        + (ib+j-prev5);                                   \
            vector2_id2 = (kb+i-prev6)*n*n       + (jb+k-prev6)*n        + (ib+j-STRIDE-prev6);                            \
            break;                                                                                                         \
        case 11: /* Vector (1,1,-1) */                                                                                     \
            prev1 = (ib > 0) ? 1 : 0;                                                                                      \
            prev2 = (ib > 0 && jb+j > STRIDE) ? 1 : 0;                                                                     \
            vector1_id0 = (kb+BLOCK_SIZE-1-k+prev1)*n*n        + (jb+j-prev1)*n        + (ib+i-prev1);                     \
            vector1_id1 = (kb+BLOCK_SIZE-1-k+prev2)*n*n        + (jb+j-STRIDE-prev2)*n + (ib+i-prev2);                     \
                                                                                                                           \
            prev3 = (jb > 0) ? 1 : 0;                                                                                      \
            prev4 = (jb > 0 && kb+BLOCK_SIZE-1-k+STRIDE < n-1) ? 1 : 0;                                                    \
            vector1_id2 = (kb+BLOCK_SIZE-1-k+prev3)*n*n        + (jb+i-prev3)*n        + (ib+j-prev3);                     \
            vector2_id0 = (kb+BLOCK_SIZE-1-k+STRIDE+prev4)*n*n + (jb+i-prev4)*n        + (ib+j-prev4);                     \
                                                                                                                           \
            prev5 = (kb+BLOCK_SIZE-1-i < n-1) ? 1 : 0;                                                                     \
            prev6 = ((kb+BLOCK_SIZE-1-i < n-1) && ib+j > STRIDE) ? 1 : 0;                                                  \
            vector2_id1 = (kb+BLOCK_SIZE-1-i+prev5)*n*n       + (jb+k-prev5)*n        + (ib+j-prev5);                      \
            vector2_id2 = (kb+BLOCK_SIZE-1-i+prev6)*n*n       + (jb+k-prev6)*n        + (ib+j-STRIDE-prev6);               \
            break;                                                                                                         \
        case 12: /* Vector (1,-1,1) */                                                                                     \
            prev1 = (ib > 0) ? 1 : 0;                                                                                      \
            prev2 = (ib > 0 && jb+BLOCK_SIZE-k-1+STRIDE < n-1) ? 1 : 0;                                                    \
            vector1_id0 = (kb+k-prev1)*n*n        + (jb+BLOCK_SIZE-k-1+prev1)*n        + (ib+i-prev1);                     \
            vector1_id1 = (kb+k-prev2)*n*n        + (jb+BLOCK_SIZE-k-1+STRIDE+prev2)*n + (ib+i-prev2);                     \
                                                                                                                           \
            prev3 = (jb+BLOCK_SIZE-i-1 < n-1) ? 1 : 0;                                                                     \
            prev4 = (jb+BLOCK_SIZE-i-1 < n-1 && kb+k > STRIDE) ? 1 : 0;                                                    \
            vector1_id2 = (kb+k-prev3)*n*n        + (jb+BLOCK_SIZE-i-1+prev3)*n        + (ib+k-prev3);                     \
            vector2_id0 = (kb+k-STRIDE-prev4)*n*n + (jb+BLOCK_SIZE-i-1+prev4)*n        + (ib+k-prev4);                     \
                                                                                                                           \
            prev5 = (kb > 0) ? 1 : 0;                                                                                      \
            prev6 = (kb > 0 && ib+k-STRIDE > 0) ? 1 : 0;                                                                   \
            vector2_id1 = (kb+i-prev5)*n*n       + (jb+BLOCK_SIZE-k-1+prev5)*n        + (ib+k-prev5);                      \
            vector2_id2 = (kb+i-prev6)*n*n       + (jb+BLOCK_SIZE-k-1+prev6)*n        + (ib+k-STRIDE-prev6);               \
            break;                                                                                                         \
        case 13: /* Vector (-1,1,1) */                                                                                     \
            prev1 = (ib+BLOCK_SIZE-1-i < n-1) ? 1 : 0;                                                                     \
            prev2 = (ib+BLOCK_SIZE-1-i < n-1 && jb+j > STRIDE) ? 1 : 0;                                                    \
            vector1_id0 = (kb+k-prev1)*n*n        + (jb+j-prev1)*n        + (ib+BLOCK_SIZE-1-i+prev1);                     \
            vector1_id1 = (kb+k-prev2)*n*n        + (jb+j-STRIDE-prev2)*n + (ib+BLOCK_SIZE-1-i+prev2);                     \
                                                                                                                           \
            prev3 = (jb > 0) ? 1 : 0;                                                                                      \
            prev4 = (jb > 0 && kb+k > STRIDE) ? 1 : 0;                                                                     \
            vector1_id2 = (kb+k-prev3)*n*n        + (jb+i-prev3)*n        + (ib+BLOCK_SIZE-1-j+prev3);                     \
            vector2_id0 = (kb+k-STRIDE-prev4)*n*n + (jb+i-prev4)*n        + (ib+BLOCK_SIZE-1-j+prev4);                     \
                                                                                                                           \
            prev5 = (kb > 0) ? 1 : 0;                                                                                      \
            prev6 = (kb > 0 && ib+BLOCK_SIZE-1-j+STRIDE < n-1) ? 1 : 0;                                                    \
            vector2_id1 = (kb+i-prev5)*n*n       + (jb+k-prev5)*n        + (ib+BLOCK_SIZE-1-j+prev5);                      \
            vector2_id2 = (kb+i-prev6)*n*n       + (jb+k-prev6)*n        + (ib+BLOCK_SIZE-1-j+STRIDE+prev6);               \
            break;                                                                                                         \
        default: ;                                                                                \
    }                                                                                             \
    double vector1_val0 = hr_sphere_region[vector1_id0];                                          \
    double vector1_val1 = hr_sphere_region[vector1_id1];                                          \
    double vector1_val2 = hr_sphere_region[vector1_id2];                                          \
    double vector1_val3 = 0.0;                                                                    \
    double vector2_val0 = hr_sphere_region[vector2_id0];                                          \
    double vector2_val1 = hr_sphere_region[vector2_id1];                                          \
    double vector2_val2 = hr_sphere_region[vector2_id2];                                          \
    double vector2_val3 = 0.0;                                                                    \
    tmp_reg1 = _mm256_set_pd(vector1_val0, vector1_val1, vector1_val2, vector1_val3);             \
    tmp_reg2 = _mm256_set_pd(vector2_val0, vector2_val1, vector2_val2, vector2_val3);             \
    prev_mask_d1 = _mm256_cmp_pd(tmp_reg1, threshold, _CMP_GT_OQ);                                \
    prev_mask_d2 = _mm256_cmp_pd(tmp_reg2, threshold, _CMP_GT_OQ);                                \
    tmp_comp1 = _mm256_castpd_si256(prev_mask_d1);                                                \
    tmp_comp2 = _mm256_castpd_si256(prev_mask_d2);                                                \
    prev_mask1 = _mm256_and_si256(tmp_comp1, ONES);                                               \
    prev_mask2 = _mm256_and_si256(tmp_comp2, ONES);


#define LOAD_DATA_REMAINDER_3D_SIMD(vecID, kb, jb, ib)                                            \
    switch (vecID) {                                                                              \
        case 10: /* Vector (1,1,1) */                                                             \
            vector1_id0 = (kb+k)*n*n        + (jb+j)*n        + (ib+i);                           \
            vector1_id1 = (kb+k)*n*n        + (jb+j-STRIDE)*n + (ib+i);                           \
                                                                                                  \
            vector1_id2 = (kb+k)*n*n        + (jb+i)*n        + (ib+j);                           \
            vector2_id0 = (kb+k-STRIDE)*n*n + (jb+i)*n        + (ib+j);                           \
                                                                                                  \
            vector2_id1 = (kb+i)*n*n       + (jb+k)*n        + (ib+j);                            \
            vector2_id2 = (kb+i)*n*n       + (jb+k)*n        + (ib+j-STRIDE);                     \
            break;                                                                                \
        case 11: /* Vector (-1,0,1) */                                                            \
            vector1_id0 = (kb+BLOCK_SIZE-1-k)*n*n        + (jb+j)*n        + (ib+i);              \
            vector1_id1 = (kb+BLOCK_SIZE-1-k)*n*n        + (jb+j-STRIDE)*n + (ib+i);              \
                                                                                                  \
            vector1_id2 = (kb+BLOCK_SIZE-1-k)*n*n        + (jb+i)*n        + (ib+j);              \
            vector2_id0 = (kb+BLOCK_SIZE-1-k+STRIDE)*n*n + (jb+i)*n        + (ib+j);              \
                                                                                                  \
            vector2_id1 = (kb+BLOCK_SIZE-1-i)*n*n       + (jb+k)*n        + (ib+j);               \
            vector2_id2 = (kb+BLOCK_SIZE-1-i)*n*n       + (jb+k)*n        + (ib+j-STRIDE);        \
            break;                                                                                \
        case 12: /* Vector (1,-1,1) */                                                            \
            vector1_id0 = (kb+k)*n*n        + (jb+BLOCK_SIZE-k-1)*n        + (ib+i);              \
            vector1_id1 = (kb+k)*n*n        + (jb+BLOCK_SIZE-k-1+STRIDE)*n + (ib+i);              \
                                                                                                  \
            vector1_id2 = (kb+k)*n*n        + (jb+BLOCK_SIZE-i-1)*n        + (ib+k);              \
            vector2_id0 = (kb+k-STRIDE)*n*n + (jb+BLOCK_SIZE-i-1)*n        + (ib+k);              \
                                                                                                  \
            vector2_id1 = (kb+i)*n*n       + (jb+BLOCK_SIZE-k-1)*n        + (ib+k);               \
            vector2_id2 = (kb+i)*n*n       + (jb+BLOCK_SIZE-k-1)*n        + (ib+k-STRIDE);        \
            break;                                                                                \
        case 13: /* Vector (-1,1,1) */                                                            \
            vector1_id0 = (kb+k)*n*n        + (jb+j)*n        + (ib+BLOCK_SIZE-1-i);              \
            vector1_id1 = (kb+k)*n*n        + (jb+j-STRIDE)*n + (ib+BLOCK_SIZE-1-i);              \
                                                                                                  \
            vector1_id2 = (kb+k)*n*n        + (jb+i)*n        + (ib+BLOCK_SIZE-1-j);              \
            vector2_id0 = (kb+k-STRIDE)*n*n + (jb+i)*n        + (ib+BLOCK_SIZE-1-j);              \
                                                                                                  \
            vector2_id1 = (kb+i)*n*n       + (jb+k)*n        + (ib+BLOCK_SIZE-1-j);               \
            vector2_id2 = (kb+i)*n*n       + (jb+k)*n        + (ib+BLOCK_SIZE-1-j+STRIDE);        \
            break;                                                                                \
        default: ;                                                                                \
    }                                                                                             \
    double vector1_val0 = hr_sphere_region[vector1_id0];                                          \
    double vector1_val1 = hr_sphere_region[vector1_id1];                                          \
    double vector1_val2 = hr_sphere_region[vector1_id2];                                          \
    double vector1_val3 = 0.0;                                                                    \
    double vector2_val0 = hr_sphere_region[vector2_id0];                                          \
    double vector2_val1 = hr_sphere_region[vector2_id1];                                          \
    double vector2_val2 = hr_sphere_region[vector2_id2];                                          \
    double vector2_val3 = 0.0;                                                                    \
    region1 = _mm256_set_pd(vector1_val0, vector1_val1, vector1_val2, vector1_val3);              \
    region2 = _mm256_set_pd(vector2_val0, vector2_val1, vector2_val2, vector2_val3);


#define BLOCK_KERNEL_3D_DIAGONALS_SIMD(kb, jb, ib)                                                \
{                                                                                                 \
    __m256d region1;                                                                              \
    __m256d bone_count = _mm256_setzero_pd();                                                     \
    __m256i edge_count = _mm256_set1_epi64x(0);                                                   \
    __m256i compi1;                                                                               \
                                                                                                  \
                                                                                                  \
    /* load prev */                                                                               \
    __m256d tmp_reg1, prev_mask_d1;                                                               \
    __m256i tmp_comp1, curr_mask1, prev_mask1;                                                    \
    int vector1_id0, vector1_id1, vector1_id2, vector1_id3, prev1, prev2, prev3, prev4;           \
    double vector1_val0, vector1_val1, vector1_val2, vector1_val3;                                \
                                                                                                  \
    prev1 = (kb >0 && jb > 0 && ib > 0) ? 1 : 0;                                                  \
    prev2 = (kb+BLOCK_SIZE-1 < n-1 && jb > 0 && ib > 0) ? 1 : 0;                                  \
    prev3 = (kb >0 && jb+BLOCK_SIZE-1 < n-1 && ib > 0) ? 1 : 0;                                   \
    prev4 = (kb >0 && jb > 0 && ib+BLOCK_SIZE-1 < n-1) ? 1 : 0;                                   \
    vector1_id0 = (kb-prev1)*n*n + (jb-prev1)*n + (ib-prev1);                                     \
    vector1_id1 = (kb+BLOCK_SIZE-1+prev2)*n*n + (jb-prev2)*n + (ib-prev2);                        \
    vector1_id2 = (kb-prev3)*n*n + (jb+BLOCK_SIZE-1+prev3)*n + (ib-prev3);                        \
    vector1_id3 = (kb-prev4)*n*n + (jb-prev4)*n + (ib+BLOCK_SIZE-1+prev4);                        \
    vector1_val0 = hr_sphere_region[vector1_id0];                                                 \
    vector1_val1 = hr_sphere_region[vector1_id1];                                                 \
    vector1_val2 = hr_sphere_region[vector1_id2];                                                 \
    vector1_val3 = hr_sphere_region[vector1_id3];                                                 \
    tmp_reg1 = _mm256_set_pd(vector1_val0, vector1_val1, vector1_val2, vector1_val3);             \
    prev_mask_d1 = _mm256_cmp_pd(tmp_reg1, threshold, _CMP_GT_OQ);                                \
    tmp_comp1 = _mm256_castpd_si256(prev_mask_d1);                                                \
    prev_mask1 = _mm256_and_si256(tmp_comp1, ONES);                                               \
                                                                                                  \
                                                                                                  \
    /* load diagonals */                                                                          \
    for (int i = 0; i < BLOCK_SIZE; ++i) {                                                        \
        vector1_id0 = (kb+i)*n*n + (jb+i)*n + (ib+i);                                             \
        vector1_id1 = (kb+BLOCK_SIZE-i-1)*n*n + (jb+i)*n + (ib+i);                                \
        vector1_id2 = (kb+i)*n*n + (jb+BLOCK_SIZE-i-1)*n + (ib+i);                                \
        vector1_id3 = (kb+i)*n*n + (jb+i)*n + (ib+BLOCK_SIZE-i-1);                                \
        vector1_val0 = hr_sphere_region[vector1_id0];                                             \
        vector1_val1 = hr_sphere_region[vector1_id1];                                             \
        vector1_val2 = hr_sphere_region[vector1_id2];                                             \
        vector1_val3 = hr_sphere_region[vector1_id3];                                             \
        region1 = _mm256_set_pd(vector1_val0, vector1_val1, vector1_val2, vector1_val3);          \
                                                                                                  \
                                                                                                  \
        compi1 = _mm256_castpd_si256(_mm256_cmp_pd(region1, threshold, _CMP_GT_OQ));              \
        curr_mask1 = _mm256_and_si256(compi1, ONES);                                              \
        __m256i edge1 = _mm256_xor_si256(curr_mask1, prev_mask1);                                  \
        bone_count = _mm256_add_pd(bone_count, region1);                                          \
        edge_count  = _mm256_add_epi64(edge_count, edge1);                                        \
        prev_mask1 = curr_mask1;                                                                  \
    }                                                                                             \
                                                                                                  \
    double bone_length_block_diagonals[4];                                                        \
    _mm256_store_pd(bone_length_block_diagonals, bone_count);                                     \
    bone_length[9] += bone_length_block_diagonals[3];                                             \
    intercepts[9] += _mm256_extract_epi64(edge_count, 3);                                         \
    bone_length[10] += bone_length_block_diagonals[2];                                            \
    intercepts[10] += _mm256_extract_epi64(edge_count, 2);                                        \
    bone_length[11] += bone_length_block_diagonals[1];                                            \
    intercepts[11] += _mm256_extract_epi64(edge_count, 1);                                        \
    bone_length[12] += bone_length_block_diagonals[0];                                            \
    intercepts[12] += _mm256_extract_epi64(edge_count, 0);                                        \
}



#define BLOCK_KERNEL_3D_SIMD(vec, kb, jb, ib)                                                             \
    {                                                                                                         \
        const int vecID = vec;                                                                                \
        double bone_length_block = 0.0;                                                                       \
        int intercepts_block = 0;                                                                             \
                                                                                                              \
        /* Init accumulators */                                                                               \
        __m256d bone_count1 = _mm256_setzero_pd();                                                            \
        __m256d bone_count2 = _mm256_setzero_pd();                                                            \
        __m256d bone_count3 = _mm256_setzero_pd();                                                            \
                                                                                                              \
        __m256i edge_count1 = _mm256_set1_epi64x(0);                                                          \
        __m256i edge_count2 = _mm256_set1_epi64x(0);                                                          \
        __m256i edge_count3 = _mm256_set1_epi64x(0);                                                          \
                                                                                                              \
        __m256i curr_mask1;                                                                                   \
        __m256i curr_mask2;                                                                                   \
        __m256i curr_mask3;                                                                                   \
                                                                                                              \
        __m256i prev_mask1;                                                                                   \
        __m256i prev_mask2;                                                                                   \
        __m256i prev_mask3;                                                                                   \
                                                                                                              \
        int prev1, prev3, prev5;                                                                              \
        int prev2, prev4, prev6;                                                                              \
                                                                                                              \
        __m256d region1, region2, region3;                                                                    \
        __m256i compi1, compi2, compi3;                                                                       \
                                                                                                              \
        for (int ks = 2; ks < BLOCK_SIZE; ks += 2*STRIDE) {                                                   \
            for (int js = ks + 2; js < BLOCK_SIZE; js += STRIDE) {                                            \
                int k1 = ks;                                                                                  \
                int j1 = js;                                                                                  \
                int k2 = js;                                                                                  \
                int j2 = ks - 2;                                                                              \
                int i = 0;                                                                                    \
                                                                                                              \
                LOAD_PREV_3D_SIMD(vecID, kb, jb, ib)                                                                      \
                                                                                                              \
                while (j1 < BLOCK_SIZE) {                                                                     \
                                                                                                              \
                    LOAD_DATA_3D_SIMD(vecID, kb, jb, ib)                                                                  \
                                                                                                              \
                    bone_count1 = _mm256_add_pd(bone_count1, region1);                                        \
                    bone_count2 = _mm256_add_pd(bone_count2, region2);                                        \
                    bone_count3 = _mm256_add_pd(bone_count3, region3);                                        \
                                                                                                              \
                    /* Calculate masks */                                                                     \
                    compi1 = _mm256_castpd_si256(_mm256_cmp_pd(region1, threshold, _CMP_GT_OQ));              \
                    curr_mask1 = _mm256_and_si256(compi1, ONES);                                              \
                    compi2 = _mm256_castpd_si256(_mm256_cmp_pd(region2, threshold, _CMP_GT_OQ));              \
                    curr_mask2 = _mm256_and_si256(compi2, ONES);                                              \
                    compi3 = _mm256_castpd_si256(_mm256_cmp_pd(region3, threshold, _CMP_GT_OQ));              \
                    curr_mask3 = _mm256_and_si256(compi3, ONES);                                              \
                    /* Detect edge and add to counter */                                                      \
                    __m256i edge1 = _mm256_xor_si256(curr_mask1, prev_mask1);                                 \
                    __m256i edge2 = _mm256_xor_si256(curr_mask2, prev_mask2);                                 \
                    __m256i edge3 = _mm256_xor_si256(curr_mask3, prev_mask3);                                 \
                                                                                                              \
                    edge_count1 = _mm256_add_epi64(edge_count1, edge1);                                       \
                    edge_count2 = _mm256_add_epi64(edge_count2, edge2);                                       \
                    edge_count3 = _mm256_add_epi64(edge_count3, edge3);                                       \
                                                                                                              \
                    prev_mask1 = curr_mask1;                                                                  \
                    prev_mask2 = curr_mask2;                                                                  \
                    prev_mask3 = curr_mask3;                                                                  \
                                                                                                              \
                    ++k1;                                                                                     \
                    ++j1;                                                                                     \
                    ++k2;                                                                                     \
                    ++j2;                                                                                     \
                    ++i;                                                                                      \
                }                                                                                             \
            }                                                                                                 \
        }                                                                                                     \
                                                                                                              \
        /* Calculate remainder */                                                                             \
        for (int jk_s = 2; jk_s < BLOCK_SIZE; jk_s += 2*STRIDE) {                                             \
            int k = jk_s;                                                                                     \
            int j = jk_s;                                                                                     \
            int i = 0;                                                                                        \
                                                                                                              \
            /* Start remainder complete vectors */                                                            \
            LOAD_PREV_REMAINDER_3D_SIMD(vecID, kb, jb, ib)                                                                \
                                                                                                              \
            while (j < BLOCK_SIZE) {                                                                          \
                LOAD_DATA_REMAINDER_3D_SIMD(vecID, kb, jb, ib)                                                            \
                                                                                                              \
                bone_count1 = _mm256_add_pd(bone_count1, region1);                                            \
                bone_count2 = _mm256_add_pd(bone_count2, region2);                                            \
                                                                                                              \
                /* Calculate masks */                                                                         \
                compi1 = _mm256_castpd_si256(_mm256_cmp_pd(region1, threshold, _CMP_GT_OQ));                  \
                curr_mask1 = _mm256_and_si256(compi1, ONES);                                                  \
                compi2 = _mm256_castpd_si256(_mm256_cmp_pd(region2, threshold, _CMP_GT_OQ));                  \
                curr_mask2 = _mm256_and_si256(compi2, ONES);                                                  \
                /* Detect edge and add to counter */                                                          \
                __m256i edge1 = _mm256_xor_si256(curr_mask1, prev_mask1);                                     \
                __m256i edge2 = _mm256_xor_si256(curr_mask2, prev_mask2);                                     \
                edge_count1 = _mm256_add_epi64(edge_count1, edge1);                                           \
                edge_count2 = _mm256_add_epi64(edge_count2, edge2);                                           \
                                                                                                              \
                prev_mask1 = curr_mask1;                                                                      \
                prev_mask2 = curr_mask2;                                                                      \
                                                                                                              \
                ++k;                                                                                          \
                ++j;                                                                                          \
                ++i;                                                                                          \
            }                                                                                                 \
        }                                                                                                     \
                                                                                                              \
        /* Sum up accumulators */                                                                             \
        __m256d bone_count23 = _mm256_add_pd(bone_count2, bone_count3);                                       \
        __m256i edge_count23 = _mm256_add_epi64(edge_count2, edge_count3);                                    \
        bone_count1 = _mm256_add_pd(bone_count1, bone_count23);                                               \
        edge_count1 = _mm256_add_epi64(edge_count1, edge_count23);                                            \
                                                                                                              \
        bone_length_block = horizontal_add(bone_count1);                                                      \
        intercepts_block  = horizontal_addi(edge_count1);                                                     \
                                                                                                              \
        bone_length[vecID-1] += bone_length_block;                                                            \
        intercepts[vecID-1]  += intercepts_block;                                                             \
    }

double simd_mil_1D(const double *hr_sphere_region, int* intercepts, int n, const int kk, const int jj, const int ii,  const int vecID);
double simd_mil_2D_pos(const double *hr_sphere_region, int* intercepts, int n, const int kk, const int jj, const int ii,  const int vecID);
double simd_mil_2D_neg(const double *hr_sphere_region, int* intercepts, int n, const int kk, const int jj, const int ii,  const int vecID);

//void mil2_simd(const double *hr_sphere_region, int n, double *directions_vectors_mil);
void mil2_simd(const double *hr_sphere_region, int n, double *directions_vectors_mil);

#endif //BONEMAP_MIL2_SIMD_H

