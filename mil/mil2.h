#ifndef BONEMAP_MIL2_H
#define BONEMAP_MIL2_H

#include <immintrin.h>
#include <stdio.h>

#define LOAD_PREV_1D                                                                     \
        switch (vecID) {                                                                 \
            case 1: /* Vector (1,0,0) */                                                 \
                i_prev = (ii > 0) ? ii - 1 : 0;                                          \
                prev_mask1 = hr_sphere_region[ k*n*n + (j+0*STRIDE)*n + i_prev ]  > 0.5; \
                prev_mask2 = hr_sphere_region[ k*n*n + (j+1*STRIDE)*n + i_prev ]  > 0.5; \
                prev_mask3 = hr_sphere_region[ k*n*n + (j+2*STRIDE)*n + i_prev ]  > 0.5; \
                prev_mask4 = hr_sphere_region[ k*n*n + (j+3*STRIDE)*n + i_prev ]  > 0.5; \
                break;                                                                   \
            case 2: /* Vector (0,1,0) */                                                 \
                i_prev = (ii > 0) ? ii - 1 : 0;                                          \
                prev_mask1 = hr_sphere_region[ k*n*n + i_prev*n + (j+0*STRIDE) ]  > 0.5; \
                prev_mask2 = hr_sphere_region[ k*n*n + i_prev*n + (j+1*STRIDE) ]  > 0.5; \
                prev_mask3 = hr_sphere_region[ k*n*n + i_prev*n + (j+2*STRIDE) ]  > 0.5; \
                prev_mask4 = hr_sphere_region[ k*n*n + i_prev*n + (j+3*STRIDE) ]  > 0.5; \
                break;                                                                   \
            case 3: /* Vector (0,0,1) */                                                 \
                i_prev = (ii > 0) ? ii - 1 : 0;                                          \
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

#define LOAD_DATA_SET_2D                                                 \
    switch (vecID) {                                                     \
        case 4: /* Vector (1,1,0) */                                     \
            /* Unroll over dimensions x,y */                             \
            r1 = hr_sphere_region[ k*n*n + j1*n + (i1+1+0*STRIDE) ];     \
            r2 = hr_sphere_region[ k*n*n + (j2+1+0*STRIDE)*n + i2 ];     \
            if (i1+1+1*STRIDE < ii + BLOCK_SIZE) {                       \
                r3 = hr_sphere_region[ k*n*n + j1*n + (i1+1+1*STRIDE) ]; \
                r4 = hr_sphere_region[ k*n*n + (j2+1+1*STRIDE)*n + i2 ]; \
            } else {                                                     \
                r3 = 0.0f;                                               \
                r4 = 0.0f;                                               \
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

#endif //BONEMAP_MIL2_H

