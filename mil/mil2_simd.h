#ifndef BONEMAP_MIL2_H
#define BONEMAP_MIL2_H

#include <immintrin.h>
#include <stdio.h>

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


#define LOAD_PREV_2D_POS                                                                          \
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

#define LOAD_DATA_SET_2D_POS                                                                      \
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



#define LOAD_PREV_2D_NEG                                                                 \
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



#define LOAD_DATA_SET_2D_NEG                                                            \
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
    __m256i compi;                                                                            \
    bone_count = _mm256_add_pd(bone_count, region1);                                          \
    bone_count2 = _mm256_add_pd(bone_count2, region2);                                        \
    \
    /* Calculate masks */                                                                     \
    compi = _mm256_castpd_si256(_mm256_cmp_pd(region1, threshold, _CMP_GT_OQ));               \
    curr_mask = _mm256_and_si256(compi, ONES);                                                \
    \
    compi = _mm256_castpd_si256(_mm256_cmp_pd(region2, threshold, _CMP_GT_OQ));               \
    curr_mask2 = _mm256_and_si256(compi, ONES);                                               \
    \
    /* Detect edge and add to counter */                                                      \
    __m256i edge1 = _mm256_xor_si256(curr_mask, prev_mask);                                   \
    __m256i edge2 = _mm256_xor_si256(curr_mask2, prev_mask2);                                 \
    edge_count = _mm256_add_epi64(edge_count, edge1);                                         \
    edge_count2 = _mm256_add_epi64(edge_count2, edge2);

double simd_mil_1D(const double *hr_sphere_region, int* intercepts, int n, const int kk, const int jj, const int ii,  const int vecID);
double simd_mil_2D_pos(const double *hr_sphere_region, int* intercepts, int n, const int kk, const int jj, const int ii,  const int vecID);
double simd_mil_2D_neg(const double *hr_sphere_region, int* intercepts, int n, const int kk, const int jj, const int ii,  const int vecID);
#endif //BONEMAP_MIL2_H

