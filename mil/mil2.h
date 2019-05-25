#ifndef BONEMAP_MIL2_H
#define BONEMAP_MIL2_H

#include <immintrin.h>
#include <stdio.h>

#define LOAD_PREV_1D                                                                     \
        i_prev = (ii > 0) ? ii - 1 : 0;                                                  \
        switch (vecID) {                                                                 \
            case 1: /* Vector (1,0,0) */                                                 \
                prev_mask1 = hr_sphere_region[ k*n*n + (j+0*STRIDE)*n + i_prev ]  > 0.5; \
                prev_mask2 = hr_sphere_region[ k*n*n + (j+1*STRIDE)*n + i_prev ]  > 0.5; \
                prev_mask3 = hr_sphere_region[ k*n*n + (j+2*STRIDE)*n + i_prev ]  > 0.5; \
                prev_mask4 = hr_sphere_region[ k*n*n + (j+3*STRIDE)*n + i_prev ]  > 0.5; \
                break;                                                                   \
            case 2: /* Vector (0,1,0) */                                                 \
                prev_mask1 = hr_sphere_region[ k*n*n + i_prev*n + (j+0*STRIDE) ]  > 0.5; \
                prev_mask2 = hr_sphere_region[ k*n*n + i_prev*n + (j+1*STRIDE) ]  > 0.5; \
                prev_mask3 = hr_sphere_region[ k*n*n + i_prev*n + (j+2*STRIDE) ]  > 0.5; \
                prev_mask4 = hr_sphere_region[ k*n*n + i_prev*n + (j+3*STRIDE) ]  > 0.5; \
                break;                                                                   \
            case 3: /* Vector (0,0,1) */                                                 \
                prev_mask1 = hr_sphere_region[ i_prev*n*n + k*n + (j+0*STRIDE) ]  > 0.5; \
                prev_mask2 = hr_sphere_region[ i_prev*n*n + k*n + (j+1*STRIDE) ]  > 0.5; \
                prev_mask3 = hr_sphere_region[ i_prev*n*n + k*n + (j+2*STRIDE) ]  > 0.5; \
                prev_mask4 = hr_sphere_region[ i_prev*n*n + k*n + (j+3*STRIDE) ]  > 0.5; \
                break;                                                                   \
            default: ;                                                                   \
        }


#define LOAD_DATA_SET_1D                                            \
        switch (vecID) {                                            \
            case 1: /* Vector (1,0,0) */                            \
                /* Unroll over dimension y */                       \
                r1 = hr_sphere_region[ k*n*n + (j+0*STRIDE)*n + i]; \
                r2 = hr_sphere_region[ k*n*n + (j+1*STRIDE)*n + i]; \
                r3 = hr_sphere_region[ k*n*n + (j+2*STRIDE)*n + i]; \
                r4 = hr_sphere_region[ k*n*n + (j+3*STRIDE)*n + i]; \
                break;                                              \
            case 2: /* Vector (0,1,0) */                            \
                /* Unroll over dimension x */                       \
                r1 = hr_sphere_region[ k*n*n + i*n + (j+0*STRIDE)]; \
                r2 = hr_sphere_region[ k*n*n + i*n + (j+1*STRIDE)]; \
                r3 = hr_sphere_region[ k*n*n + i*n + (j+2*STRIDE)]; \
                r4 = hr_sphere_region[ k*n*n + i*n + (j+3*STRIDE)]; \
                break;                                              \
            case 3: /* Vector (0,0,1) */                            \
                /* Unroll over dimension x */                       \
                r1 = hr_sphere_region[ i*n*n + k*n + (j+0*STRIDE)]; \
                r2 = hr_sphere_region[ i*n*n + k*n + (j+1*STRIDE)]; \
                r3 = hr_sphere_region[ i*n*n + k*n + (j+2*STRIDE)]; \
                r4 = hr_sphere_region[ i*n*n + k*n + (j+3*STRIDE)]; \
                break;                                              \
            default: ;                                              \
        }

#define LOAD_PREV_2D_POS                                                                          \
    if (j1 > 0) {           \
        i1_prev = i1 - 1;   \
        j1_prev = j1 - 1;   \
    }                       \
    else {                  \
        i1_prev = i1;       \
        j1_prev = j1;       \
    }                       \
    if (i2 > 0) {           \
        i2_prev = i2 - 1;   \
        j2_prev = j2 - 1;   \
    }                       \
    else {                  \
        i2_prev = i2;       \
        j2_prev = j2;       \
    }                       \
    switch (vecID) {                                                                              \
        case 4: /* Vector (1,1,0) */                                                              \
            /* Unroll over dimensions x,y */                                                      \
            prev_mask1 = hr_sphere_region[ k*n*n + j1_prev*n + (i1_prev+1+0*STRIDE) ]  > 0.5;     \
            prev_mask2 = hr_sphere_region[ k*n*n + (j2_prev+1+0*STRIDE)*n + i2_prev ]  > 0.5;     \
            prev_mask3 = hr_sphere_region[ k*n*n + j1_prev*n + (i1_prev+1+1*STRIDE) ]  > 0.5;     \
            prev_mask4 = hr_sphere_region[ k*n*n + (j2_prev+1+1*STRIDE)*n + i2_prev ]  > 0.5;     \
            break;                                                                                \
        case 5: /* Vector (1,0,1) */                                                              \
            /* Unroll over dimensions x,z */                                                      \
            prev_mask1 = hr_sphere_region[ j1_prev*n*n + k*n + (i1_prev+1+0*STRIDE) ]  > 0.5;     \
            prev_mask2 = hr_sphere_region[ (j2_prev+1+0*STRIDE)*n*n + k*n + i2_prev ]  > 0.5;     \
            prev_mask3 = hr_sphere_region[ j1_prev*n*n + k*n + (i1_prev+1+1*STRIDE) ]  > 0.5;     \
            prev_mask4 = hr_sphere_region[ (j2_prev+1+1*STRIDE)*n*n + k*n + i2_prev ]  > 0.5;     \
            break;                                                                                \
        case 6: /* Vector (0,1,1) */                                                              \
            /* Unroll over dimensions x,z */                                                      \
            prev_mask1 = hr_sphere_region[ (i1_prev+1+0*STRIDE)*n*n + j1_prev*n + k ]  > 0.5;     \
            prev_mask2 = hr_sphere_region[ i2_prev*n*n + (j2_prev+1+0*STRIDE)*n + k ]  > 0.5;     \
            prev_mask3 = hr_sphere_region[ (i1_prev+1+1*STRIDE)*n*n + j1_prev*n + k ]  > 0.5;     \
            prev_mask4 = hr_sphere_region[ i2_prev*n*n + (j2_prev+1+1*STRIDE)*n + k ]  > 0.5;     \
            break;                                                                                \
        default: ;                                                                                \
    }

#define LOAD_DATA_SET_2D_POS                                             \
    switch (vecID) {                                                     \
        case 4: /* Vector (1,1,0) */                                     \
            /* Unroll over dimensions x,y */                             \
            r1 = hr_sphere_region[ k*n*n + j1*n + (i1+1+0*STRIDE) ];     \
            r2 = hr_sphere_region[ k*n*n + (j2+1+0*STRIDE)*n + i2 ];     \
            if (i1+1+1*STRIDE < ii + BLOCK_SIZE) {                       \
                r3 = hr_sphere_region[ k*n*n + j1*n + (i1+1+1*STRIDE) ]; \
                r4 = hr_sphere_region[ k*n*n + (j2+1+1*STRIDE)*n + i2 ]; \
            } else {                                                     \
                r3 = 0.0;                                                \
                r4 = 0.0;                                                \
                prev_mask3 = 0;                                          \
                prev_mask4 = 0;                                          \
            }                                                            \
            break;                                                       \
        case 5: /* Vector (1,0,1) */                                     \
            /* Unroll over dimensions x,z */                             \
            r1 = hr_sphere_region[ j1*n*n + k*n + (i1+1+0*STRIDE) ];     \
            r2 = hr_sphere_region[ (j2+1+0*STRIDE)*n*n + k*n + i2 ];     \
            if (i1+1+1*STRIDE < ii + BLOCK_SIZE) {                       \
                r3 = hr_sphere_region[ j1*n*n + k*n + (i1+1+1*STRIDE) ]; \
                r4 = hr_sphere_region[ (j2+1+1*STRIDE)*n*n + k*n + i2 ]; \
            } else {                                                     \
                r3 = 0.0f;                                               \
                r4 = 0.0f;                                               \
                prev_mask3 = 0;                                          \
                prev_mask4 = 0;                                          \
            }                                                            \
            break;                                                       \
        case 6: /* Vector (0,1,1) */                                     \
            /* Unroll over dimensions y,z */                             \
            r1 = hr_sphere_region[ (i1+1+0*STRIDE)*n*n + j1*n + k ];     \
            r2 = hr_sphere_region[ i2*n*n + (j2+1+0*STRIDE)*n + k ];     \
            if (i1+1+1*STRIDE < ii + BLOCK_SIZE) {                       \
                r3 = hr_sphere_region[ (i1+1+1*STRIDE)*n*n + j1*n + k ]; \
                r4 = hr_sphere_region[ i2*n*n + (j2+1+1*STRIDE)*n + k ]; \
            } else {                                                     \
                r3 = 0.0f;                                               \
                r4 = 0.0f;                                               \
                prev_mask3 = 0;                                          \
                prev_mask4 = 0;                                          \
            }                                                            \
            break;                                                       \
        default: ;                                                       \
    }


#define LOAD_DATA_SET_2D_NEG                                             \
    switch (vecID) {                                                     \
        case 7: /* Vector (-1,1,0) */                                    \
            /* Unroll over dimensions x,y */                             \
            r1 = hr_sphere_region[ k*n*n + j1*n + (i1-1-0*STRIDE) ];     \
            r2 = hr_sphere_region[ k*n*n + (j2+1+0*STRIDE)*n + i2 ];     \
            if (j2+1+1*STRIDE < jj + BLOCK_SIZE) {                       \
                r3 = hr_sphere_region[ k*n*n + j1*n + (i1-1-1*STRIDE) ]; \
                r4 = hr_sphere_region[ k*n*n + (j2+1+1*STRIDE)*n + i2 ]; \
            } else {                                                     \
                r3 = 0.0f;                                               \
                r4 = 0.0f;                                               \
                prev_mask3 = 0;                                          \
                prev_mask4 = 0;                                          \
            }                                                            \
            break;                                                       \
        case 8: /* Vector (-1,0,1) */                                    \
            /* Unroll over dimensions x,y */                             \
            r1 = hr_sphere_region[ j1*n*n + k*n + (i1-1-0*STRIDE) ];     \
            r2 = hr_sphere_region[ (j2+1+0*STRIDE)*n*n + k*n + i2 ];     \
            if (j2+1+1*STRIDE < jj + BLOCK_SIZE) {                       \
                r3 = hr_sphere_region[ j1*n*n + k*n + (i1-1-1*STRIDE) ]; \
                r4 = hr_sphere_region[ (j2+1+1*STRIDE)*n*n + k*n + i2 ]; \
            } else {                                                     \
                r3 = 0.0f;                                               \
                r4 = 0.0f;                                               \
                prev_mask3 = 0;                                          \
                prev_mask4 = 0;                                          \
            }                                                            \
            break;                                                       \
        case 9: /* Vector (0,1,-1) */                                    \
            /* Unroll over dimensions x,y */                             \
            r1 = hr_sphere_region[ (i1-1-0*STRIDE)*n*n + j1*n + k ];     \
            r2 = hr_sphere_region[ i2*n*n + (j2+1+0*STRIDE)*n + k ];     \
            if (j2+1+1*STRIDE < jj + BLOCK_SIZE) {                       \
                r3 = hr_sphere_region[ (i1-1-1*STRIDE)*n*n + j1*n + k ]; \
                r4 = hr_sphere_region[ i2*n*n + (j2+1+1*STRIDE)*n + k ]; \
            } else {                                                     \
                r3 = 0.0f;                                               \
                r4 = 0.0f;                                               \
                prev_mask3 = 0;                                          \
                prev_mask4 = 0;                                          \
            }                                                            \
            break;                                                       \
        default: ;                                                       \
    }
    // End of LOAD_DATA_SET_2D_NEG

#define COMPUTATION \
        unsigned int curr_mask1;                \
        unsigned int curr_mask2;                \
        unsigned int curr_mask3;                \
        unsigned int curr_mask4;                \
        acc1 += r1;                             \
        acc2 += r2;                             \
        acc3 += r3;                             \
        acc4 += r4;                             \
                                                \
        /* Calculate masks */                   \
        curr_mask1 = r1 > 0.5f;                 \
        curr_mask2 = r2 > 0.5f;                 \
        curr_mask3 = r3 > 0.5f;                 \
        curr_mask4 = r4 > 0.5f;                 \
                                                \
        /* Detect edge and add to counter */    \
        edge_count1 += curr_mask1 ^ prev_mask1; \
        edge_count2 += curr_mask2 ^ prev_mask2; \
        edge_count3 += curr_mask3 ^ prev_mask3; \
        edge_count4 += curr_mask4 ^ prev_mask4;
        // End of COMPUTATION


/**
 * SIMD IMPLEMNTATION
 */
//
//#define GATHER_SCALE 8
//#define ONES _mm256_set1_epi64x(1)
//#define threshold _mm256_set1_pd(0.5)
//
//
//#define SIMD_LOAD_PREV_1D                                                                \
//        int offset = STRIDE * n; \
//        __m128i bytes_offset_vector1 = _mm_set_epi32(i_prev, offset + i_prev, 2 * offset + i_prev, 3 * offset + i_prev); \
//        __m128i bytes_offset_vector2 = _mm_set_epi32(0, STRIDE, 2 * STRIDE, 3 * STRIDE); \
//        __m128i bytes_offset_vector3 = _mm_set_epi32(0, STRIDE, 2 * STRIDE, 3 * STRIDE); \
//        __m256d tmp_reg1, tmp_reg2, prev_mask_d1, prev_mask_d2;                                                    \
//        __m256i tmp_comp1, tmp_comp2; \
//        switch (vecID) {                                                                 \
//            case 1: /* Vector (1,0,0) */                                                 \
//                tmp_reg1 = _mm256_i32gather_pd(hr_sphere_region + (knn + j * n), bytes_offset_vector1, GATHER_SCALE); \
//                tmp_reg2 = _mm256_i32gather_pd(hr_sphere_region + (knn + (j + STRIDE * NUM_ACC) * n), bytes_offset_vector1, GATHER_SCALE); \
//                break;                                                                   \
//            case 2: /* Vector (0,1,0) */                                                 \
//                tmp_reg1 = _mm256_i32gather_pd(hr_sphere_region + knn + i_prev * n, bytes_offset_vector2, GATHER_SCALE); \
//                tmp_reg2 = _mm256_i32gather_pd(hr_sphere_region + knn + i_prev * n + j + STRIDE * NUM_ACC, bytes_offset_vector2, GATHER_SCALE); \
//                break;                                                                   \
//            case 3: /* Vector (0,0,1) */                                                 \
//                tmp_reg1 = _mm256_i32gather_pd(hr_sphere_region + (i_prev*n*n + (k * n) + j), bytes_offset_vector3, GATHER_SCALE); \
//                tmp_reg2 = _mm256_i32gather_pd(hr_sphere_region + (i_prev*n*n + k*n + j + STRIDE * NUM_ACC), bytes_offset_vector3, GATHER_SCALE); \
//                break;                                                                   \
//            default: ;                                                                   \
//        }\
//        prev_mask_d1 = _mm256_cmp_pd(tmp_reg1, threshold, _CMP_GT_OQ);\
//        prev_mask_d2 = _mm256_cmp_pd(tmp_reg2, threshold, _CMP_GT_OQ);\
//        tmp_comp1 = _mm256_castpd_si256(prev_mask_d1); \
//        tmp_comp2 = _mm256_castpd_si256(prev_mask_d2); \
//        prev_mask = _mm256_and_si256(tmp_comp1, ONES);\
//        prev_mask2 = _mm256_and_si256(tmp_comp2, ONES);
//
//
//#define SIMD_LOAD_DATA_SET_1D                                       \
//        __m128i bytes_offset_vector1 = _mm_set_epi32(i, offset + i, 2 * offset + i, 3 * offset + i); \
//        __m128i bytes_offset_vector2 = _mm_set_epi32(0, STRIDE, 2 * STRIDE, 3 * STRIDE); \
//        __m128i bytes_offset_vector3 = _mm_set_epi32(0, STRIDE, 2 * STRIDE, 3 * STRIDE); \
//        switch (vecID) {                                            \
//            case 1: /* Vector (1,0,0) */                            \
//                /* Unroll over dimension x */                       \
//                region1 = _mm256_i32gather_pd(hr_sphere_region + (knn + j * n), bytes_offset_vector1, GATHER_SCALE); \
//                region2 = _mm256_i32gather_pd(hr_sphere_region + (knn + (j + STRIDE * NUM_ACC) * n), bytes_offset_vector1, GATHER_SCALE); \
//                /*region1 = _mm256_set_pd(hr_sphere_region[ k*n*n + (j+0*STRIDE)*n + i], hr_sphere_region[ k*n*n + (j+1*STRIDE)*n + i], hr_sphere_region[ k*n*n + (j+2*STRIDE)*n + i], hr_sphere_region[ k*n*n + (j+3*STRIDE)*n + i]); \
//                region2 = _mm256_set_pd(hr_sphere_region[ k*n*n + (j + STRIDE * NUM_ACC) * n + i], hr_sphere_region[ k*n*n + (j + STRIDE* NUM_ACC + STRIDE)*n + i], hr_sphere_region[ k*n*n + (j+STRIDE* NUM_ACC+ 2 * STRIDE)*n + i], hr_sphere_region[ k*n*n + (j+STRIDE* NUM_ACC + 3 * STRIDE)*n + i]);*/ \
//                break;                                              \
//            case 2: /* Vector (0,1,0) */                            \
//                /* Unroll over dimension y */                       \
//                region1 = _mm256_i32gather_pd(hr_sphere_region + knn + i * n, bytes_offset_vector2, GATHER_SCALE); \
//                region2 = _mm256_i32gather_pd(hr_sphere_region + knn + i * n + j + STRIDE * NUM_ACC, bytes_offset_vector2, GATHER_SCALE); \
//                /*region1 = _mm256_set_pd(hr_sphere_region[ k*n*n + i*n + (j+0*STRIDE)], hr_sphere_region[ k*n*n + i*n + (j+1*STRIDE)], hr_sphere_region[ k*n*n + i*n + (j+2*STRIDE)], hr_sphere_region[ k*n*n + i*n + (j+3*STRIDE)]); \
//                region2 = _mm256_set_pd(hr_sphere_region[ k*n*n + i*n + (j + STRIDE * NUM_ACC +0*STRIDE)], hr_sphere_region[ k*n*n + i*n + (j + STRIDE * NUM_ACC +1*STRIDE)], hr_sphere_region[ k*n*n + i*n + (j + STRIDE * NUM_ACC +2*STRIDE)], hr_sphere_region[ k*n*n + i*n + (j + STRIDE * NUM_ACC +3*STRIDE)]);*/ \
//                break;                                              \
//            case 3: /* Vector (0,0,1) */                            \
//                /* Unroll over dimension z */                       \
//                region1 = _mm256_i32gather_pd(hr_sphere_region + (inn + (k * n) + j), bytes_offset_vector3, GATHER_SCALE); \
//                region2 = _mm256_i32gather_pd(hr_sphere_region + (inn + k*n + j + STRIDE * NUM_ACC), bytes_offset_vector3, GATHER_SCALE); \
//                /*region1 = _mm256_set_pd(hr_sphere_region[ i*n*n + k*n + (j+0*STRIDE)], hr_sphere_region[ i*n*n + k*n + (j+1*STRIDE)], hr_sphere_region[ i*n*n + k*n + (j+2*STRIDE)], hr_sphere_region[ i*n*n + k*n + (j+3*STRIDE)]); \
//                region2 = _mm256_set_pd(hr_sphere_region[ i*n*n + k*n + (j + STRIDE * NUM_ACC +0*STRIDE)], hr_sphere_region[ i*n*n + k*n + (j + STRIDE * NUM_ACC +1*STRIDE)], hr_sphere_region[ i*n*n + k*n + (j + STRIDE * NUM_ACC +2*STRIDE)], hr_sphere_region[ i*n*n + k*n + (j + STRIDE * NUM_ACC +3*STRIDE)]);*/ \
//                break;                                              \
//            default: ;                                              \
//        }
//
//
//
//#define SIMD_LOAD_DATA_SET_2D                                            \
//    switch (vecID) {                                                     \
//        case 4: /* Vector (1,1,0) */                                     \
//            /* Unroll over dimensions x,y */                             \
//            r1 = hr_sphere_region[ k*n*n + j1*n + (i1+1+0*STRIDE) ];     \
//            r2 = hr_sphere_region[ k*n*n + (j2+1+0*STRIDE)*n + i2 ];     \
//            if (i1+1+1*STRIDE < ii + BLOCK_SIZE) {                       \
//                r3 = hr_sphere_region[ k*n*n + j1*n + (i1+1+1*STRIDE) ]; \
//                r4 = hr_sphere_region[ k*n*n + (j2+1+1*STRIDE)*n + i2 ]; \
//            } else {                                                     \
//                r3 = 0.0f;                                               \
//                r4 = 0.0f;                                               \
//                prev_mask3 = 0;                                          \
//                prev_mask4 = 0;                                          \
//            }                                                            \
//            break;                                                       \
//        case 5: /* Vector (1,0,1) */                                     \
//            /* Unroll over dimensions x,z */                             \
//            r1 = hr_sphere_region[ j1*n*n + k*n + (i1+1+0*STRIDE) ];     \
//            r2 = hr_sphere_region[ (j2+1+0*STRIDE)*n*n + k*n + i2 ];     \
//            if (i1+1+1*STRIDE < ii + BLOCK_SIZE) {                       \
//                r3 = hr_sphere_region[ j1*n*n + k*n + (i1+1+1*STRIDE) ]; \
//                r4 = hr_sphere_region[ (j2+1+1*STRIDE)*n*n + k*n + i2 ]; \
//            } else {                                                     \
//                r3 = 0.0f;                                               \
//                r4 = 0.0f;                                               \
//                prev_mask3 = 0;                                          \
//                prev_mask4 = 0;                                          \
//            }                                                            \
//            break;                                                       \
//        case 6: /* Vector (0,1,1) */                                     \
//            /* Unroll over dimensions y,z */                             \
//            r1 = hr_sphere_region[ (i1+1+0*STRIDE)*n*n + j1*n + k ];     \
//            r2 = hr_sphere_region[ i2*n*n + (j2+1+0*STRIDE)*n + k ];     \
//            if (i1+1+1*STRIDE < ii + BLOCK_SIZE) {                       \
//                r3 = hr_sphere_region[ (i1+1+1*STRIDE)*n*n + j1*n + k ]; \
//                r4 = hr_sphere_region[ i2*n*n + (j2+1+1*STRIDE)*n + k ]; \
//            } else {                                                     \
//                r3 = 0.0f;                                               \
//                r4 = 0.0f;                                               \
//                prev_mask3 = 0;                                          \
//                prev_mask4 = 0;                                          \
//            }                                                            \
//            break;                                                       \
//        case 7: /* Vector (-1,1,0) */                                    \
//            /* Unroll over dimensions x,y */                             \
//            r1 = hr_sphere_region[ k*n*n + j1*n + (i1-1-0*STRIDE) ];     \
//            r2 = hr_sphere_region[ k*n*n + (j2+1+0*STRIDE)*n + i2 ];     \
//            if (j2+1+1*STRIDE < jj + BLOCK_SIZE) {                       \
//                r3 = hr_sphere_region[ k*n*n + j1*n + (i1-1-1*STRIDE) ]; \
//                r4 = hr_sphere_region[ k*n*n + (j2+1+1*STRIDE)*n + i2 ]; \
//            } else {                                                     \
//                r3 = 0.0f;                                               \
//                r4 = 0.0f;                                               \
//                prev_mask3 = 0;                                          \
//                prev_mask4 = 0;                                          \
//            }                                                            \
//            break;                                                       \
//        case 8: /* Vector (-1,0,1) */                                    \
//            /* Unroll over dimensions x,y */                             \
//            r1 = hr_sphere_region[ j1*n*n + k*n + (i1-1-0*STRIDE) ];     \
//            r2 = hr_sphere_region[ (j2+1+0*STRIDE)*n*n + k*n + i2 ];     \
//            if (j2+1+1*STRIDE < jj + BLOCK_SIZE) {                       \
//                r3 = hr_sphere_region[ j1*n*n + k*n + (i1-1-1*STRIDE) ]; \
//                r4 = hr_sphere_region[ (j2+1+1*STRIDE)*n*n + k*n + i2 ]; \
//            } else {                                                     \
//                r3 = 0.0f;                                               \
//                r4 = 0.0f;                                               \
//                prev_mask3 = 0;                                          \
//                prev_mask4 = 0;                                          \
//            }                                                            \
//            break;                                                       \
//        case 9: /* Vector (0,1,-1) */                                    \
//            /* Unroll over dimensions x,y */                             \
//            r1 = hr_sphere_region[ (i1-1-0*STRIDE)*n*n + j1*n + k ];     \
//            r2 = hr_sphere_region[ i2*n*n + (j2+1+0*STRIDE)*n + k ];     \
//            if (j2+1+1*STRIDE < jj + BLOCK_SIZE) {                       \
//                r3 = hr_sphere_region[ (i1-1-1*STRIDE)*n*n + j1*n + k ]; \
//                r4 = hr_sphere_region[ i2*n*n + (j2+1+1*STRIDE)*n + k ]; \
//            } else {                                                     \
//                r3 = 0.0f;                                               \
//                r4 = 0.0f;                                               \
//                prev_mask3 = 0;                                          \
//                prev_mask4 = 0;                                          \
//            }                                                            \
//            break;                                                       \
//        default: ;                                                       \
//    }
//
//
//#define SIMD_COMPUTATION \
//        __m256i curr_mask, curr_mask2;          \
//        __m256i compi; \
//        bone_count = _mm256_add_pd(bone_count, region1);                             \
//        bone_count2 = _mm256_add_pd(bone_count2, region2);                             \
//                                                \
//        /* Calculate masks */                   \
//        compi = _mm256_castpd_si256(_mm256_cmp_pd(region1, threshold, _CMP_GT_OQ));\
//        curr_mask = _mm256_and_si256(compi, ONES); \
//        \
//        compi = _mm256_castpd_si256(_mm256_cmp_pd(region2, threshold, _CMP_GT_OQ));\
//        curr_mask2 = _mm256_and_si256(compi, ONES); \
//                                                \
//        /* Detect edge and add to counter */    \
//        __m256i edge1 = _mm256_xor_si256(curr_mask, prev_mask); \
//        __m256i edge2 = _mm256_xor_si256(curr_mask2, prev_mask2); \
//        edge_count = _mm256_add_epi64(edge_count, edge1); \
//        edge_count2 = _mm256_add_epi64(edge_count2, edge2);
//
//
////void mil2(const int *hr_sphere_region, int n, double *directions_vectors_mil);
//double simd_mil_1D(const double *hr_sphere_region, int* intercepts, int n, const int kk, const int jj, const int ii,  const int vecID);
#endif //BONEMAP_MIL2_H

