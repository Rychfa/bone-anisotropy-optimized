#ifndef BONEMAP_MIL2_H
#define BONEMAP_MIL2_H

#include <immintrin.h>
#include <stdio.h>

#define BLOCK_SIZE 16
#define NUM_ACC 4
#define STRIDE 2

#define LOAD_PREV_1D \
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

#define LOAD_DATA_SET_1D \
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

#define LOAD_PREV_2D_POS \
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
            prev_mask1 = hr_sphere_region[     k*n*n + j1_prev*n + (i1_prev+1) ]  > 0.5;          \
            prev_mask2 = hr_sphere_region[ (k+STRIDE)*n*n + j1_prev*n + (i1_prev+1) ]  > 0.5;     \
            prev_mask3 = hr_sphere_region[     k*n*n + (j2_prev+1)*n + i2_prev ]  > 0.5;          \
            prev_mask4 = hr_sphere_region[ (k+STRIDE)*n*n + (j2_prev+1)*n + i2_prev ]  > 0.5;     \
            break;                                                                                \
        case 5: /* Vector (1,0,1) */                                                              \
            /* Unroll over dimensions x,z */                                                      \
            prev_mask1 = hr_sphere_region[ j1_prev*n*n +     k*n + (i1_prev+1) ]  > 0.5;          \
            prev_mask2 = hr_sphere_region[ j1_prev*n*n + (k+STRIDE)*n + (i1_prev+1) ]  > 0.5;     \
            prev_mask3 = hr_sphere_region[ (j2_prev+1)*n*n +     k*n + i2_prev ]  > 0.5;          \
            prev_mask4 = hr_sphere_region[ (j2_prev+1)*n*n + (k+STRIDE)*n + i2_prev ]  > 0.5;     \
            break;                                                                                \
        case 6: /* Vector (0,1,1) */                                                              \
            /* Unroll over dimensions x,z */                                                      \
            prev_mask1 = hr_sphere_region[ (i1_prev+1)*n*n + j1_prev*n + k     ]  > 0.5;          \
            prev_mask2 = hr_sphere_region[ (i1_prev+1)*n*n + j1_prev*n + (k+STRIDE) ]  > 0.5;     \
            prev_mask3 = hr_sphere_region[ i2_prev*n*n + (j2_prev+1)*n + k     ]  > 0.5;          \
            prev_mask4 = hr_sphere_region[ i2_prev*n*n + (j2_prev+1)*n + (k+STRIDE) ]  > 0.5;     \
            break;                                                                                \
        default: ;                                                                                \
    }

#define LOAD_DATA_SET_2D_POS \
    switch (vecID) {                                                    \
        case 4: /* Vector (1,1,0) */                                    \
            /* Unroll over dimensions x,y */                            \
            r1 = hr_sphere_region[     k*n*n + j1*n + (i1+1) ];         \
            r2 = hr_sphere_region[ (k+STRIDE)*n*n + j1*n + (i1+1) ];    \
            r3 = hr_sphere_region[     k*n*n + (j2+1)*n + i2 ];         \
            r4 = hr_sphere_region[ (k+STRIDE)*n*n + (j2+1)*n + i2 ];    \
            break;                                                      \
        case 5: /* Vector (1,0,1) */                                    \
            /* Unroll over dimensions x,z */                            \
            r1 = hr_sphere_region[ j1*n*n +     k*n + (i1+1) ];         \
            r2 = hr_sphere_region[ j1*n*n + (k+STRIDE)*n + (i1+1) ];    \
            r3 = hr_sphere_region[ (j2+1)*n*n +     k*n + i2 ];         \
            r4 = hr_sphere_region[ (j2+1)*n*n + (k+STRIDE)*n + i2 ];    \
            break;                                                  \
        case 6: /* Vector (0,1,1) */                                \
            /* Unroll over dimensions y,z */                        \
            r1 = hr_sphere_region[ (i1+1)*n*n + j1*n + k     ];     \
            r2 = hr_sphere_region[ (i1+1)*n*n + j1*n + (k+STRIDE) ];     \
            r3 = hr_sphere_region[ i2*n*n + (j2+1)*n + k     ];     \
            r4 = hr_sphere_region[ i2*n*n + (j2+1)*n + (k+STRIDE) ];     \
            break;                                                  \
        default: ;                                                  \
    }

#define LOAD_PREV_2D_POS_UNEVEN \
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

#define LOAD_DATA_SET_2D_POS_UNEVEN \
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

#define LOAD_PREV_2D_NEG    \
    if (j1 > 0) {           \
        i1_prev = i1 + 1;   \
        j1_prev = j1 - 1;   \
    }                       \
    else {                  \
        i1_prev = i1;       \
        j1_prev = j1;       \
    }                       \
    if (i2 < n-1) {         \
        i2_prev = i2 + 1;   \
        j2_prev = j2 - 1;   \
    }                       \
    else {                  \
        i2_prev = i2;       \
        j2_prev = j2;       \
    }                       \
    switch (vecID) {                                                                              \
        case 7: /* Vector (-1,1,0) */                                                             \
            /* Unroll over dimensions x,y */                                                      \
            prev_mask1 = hr_sphere_region[  k*n*n         +  j1_prev*n    + (i1_prev-1) ] > 0.5;  \
            prev_mask2 = hr_sphere_region[ (k+STRIDE)*n*n +  j1_prev*n    + (i1_prev-1) ] > 0.5;  \
            prev_mask3 = hr_sphere_region[  k*n*n         + (j2_prev+1)*n + i2_prev ] > 0.5;      \
            prev_mask4 = hr_sphere_region[ (k+STRIDE)*n*n + (j2_prev+1)*n + i2_prev ] > 0.5;      \
            break;                                                                                \
        case 8: /* Vector (-1,0,1) */                                                             \
            /* Unroll over dimensions x,z */                                                      \
            prev_mask1 = hr_sphere_region[  j1_prev*n*n    +  k*n         + (i1_prev-1) ]  > 0.5; \
            prev_mask2 = hr_sphere_region[  j1_prev*n*n    + (k+STRIDE)*n + (i1_prev-1) ]  > 0.5; \
            prev_mask3 = hr_sphere_region[ (j2_prev+1)*n*n +  k*n         + i2_prev ]  > 0.5;     \
            prev_mask4 = hr_sphere_region[ (j2_prev+1)*n*n + (k+STRIDE)*n + i2_prev ]  > 0.5;     \
            break;                                                                                \
        case 9: /* Vector (0,1,-1) */                                                             \
            /* Unroll over dimensions x,z */                                                      \
            prev_mask1 = hr_sphere_region[ (i1_prev-1)*n*n + j1_prev*n     +  k ]          > 0.5; \
            prev_mask2 = hr_sphere_region[ (i1_prev-1)*n*n + j1_prev*n     + (k+STRIDE) ]  > 0.5; \
            prev_mask3 = hr_sphere_region[ i2_prev*n*n     + (j2_prev+1)*n +  k ]          > 0.5; \
            prev_mask4 = hr_sphere_region[ i2_prev*n*n     + (j2_prev+1)*n + (k+STRIDE) ]  > 0.5; \
            break;                                                                                \
        default: ;                                                                                \
    }

#define LOAD_DATA_SET_2D_NEG \
    switch (vecID) {                                                      \
        case 7: /* Vector (-1,1,0) */                                     \
            /* Unroll over dimensions x,y */                              \
            r1 = hr_sphere_region[  k*n*n         +  j1*n    + (i1-1) ];  \
            r2 = hr_sphere_region[ (k+STRIDE)*n*n +  j1*n    + (i1-1) ];  \
            r3 = hr_sphere_region[  k*n*n         + (j2+1)*n + i2 ];      \
            r4 = hr_sphere_region[ (k+STRIDE)*n*n + (j2+1)*n + i2 ];      \
            break;                                                        \
        case 8: /* Vector (-1,0,1) */                                     \
            /* Unroll over dimensions x,y */                              \
            r1 = hr_sphere_region[  j1*n*n    +  k*n         + (i1-1) ];  \
            r2 = hr_sphere_region[  j1*n*n    + (k+STRIDE)*n + (i1-1) ];  \
            r3 = hr_sphere_region[ (j2+1)*n*n +  k*n         + i2 ]  ;    \
            r4 = hr_sphere_region[ (j2+1)*n*n + (k+STRIDE)*n + i2 ]  ;    \
            break;                                                        \
        case 9: /* Vector (0,1,-1) */                                     \
            /* Unroll over dimensions x,y */                              \
            r1 = hr_sphere_region[ (i1-1)*n*n + j1*n     +  k ]        ;  \
            r2 = hr_sphere_region[ (i1-1)*n*n + j1*n     + (k+STRIDE) ];  \
            r3 = hr_sphere_region[ i2*n*n     + (j2+1)*n +  k ]        ;  \
            r4 = hr_sphere_region[ i2*n*n     + (j2+1)*n + (k+STRIDE) ];  \
            break;                                                        \
        default: ;                                                        \
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
        curr_mask1 = r1 > 0.5;                  \
        curr_mask2 = r2 > 0.5;                  \
        curr_mask3 = r3 > 0.5;                  \
        curr_mask4 = r4 > 0.5;                  \
                                                \
        /* Detect edge and add to counter */    \
        edge_count1 += curr_mask1 ^ prev_mask1; \
        edge_count2 += curr_mask2 ^ prev_mask2; \
        edge_count3 += curr_mask3 ^ prev_mask3; \
        edge_count4 += curr_mask4 ^ prev_mask4;
        // End of COMPUTATION

#define BLOCK_KERNEL_1D(vec, kk, jj, ii) \
    {                                                                              \
        const int vecID = vec;                                                     \
        double bone_length_block = 0.0;                                            \
        int intercepts_block = 0;                                                  \
                                                                                   \
        /* Init accumulators */                                                    \
        double acc1 = 0.0;                                                         \
        double acc2 = 0.0;                                                         \
        double acc3 = 0.0;                                                         \
        double acc4 = 0.0;                                                         \
                                                                                   \
        unsigned int edge_count1 = 0;                                              \
        unsigned int edge_count2 = 0;                                              \
        unsigned int edge_count3 = 0;                                              \
        unsigned int edge_count4 = 0;                                              \
                                                                                   \
        unsigned int curr_mask1;                                                   \
        unsigned int curr_mask2;                                                   \
        unsigned int curr_mask3;                                                   \
        unsigned int curr_mask4;                                                   \
                                                                                   \
        for (int k = kk + 1; k < kk + BLOCK_SIZE; k += STRIDE) {                   \
            for (int j = jj + 1; j < jj + BLOCK_SIZE; j += STRIDE*NUM_ACC) {           \
                unsigned int i_prev, prev_mask1, prev_mask2, prev_mask3, prev_mask4;   \
                i_prev = (ii > 0) ? ii - 1 : 0;                                        \
                LOAD_PREV_1D                                                           \
                                                                                       \
                for (int i = ii; i < ii + BLOCK_SIZE; ++i) {                           \
                    double r1, r2, r3, r4;                                             \
                                                                                       \
                    /* Load working set */                                             \
                    LOAD_DATA_SET_1D                                                   \
                                                                                       \
                    /* Perform computation */                                          \
                    COMPUTATION                                                        \
                                                                                       \
                    /* Update state of prev_mask */                                    \
                    prev_mask1 = curr_mask1;                                           \
                    prev_mask2 = curr_mask2;                                           \
                    prev_mask3 = curr_mask3;                                           \
                    prev_mask4 = curr_mask4;                                           \
                }                                                                      \
            }                                                                          \
        }                                                                              \
                                                                                       \
        acc1 += acc2;                                                                  \
        acc3 += acc4;                                                                  \
        edge_count1 += edge_count2;                                                    \
        edge_count3 += edge_count4;                                                    \
                                                                                       \
        bone_length_block = acc1 + acc3;                                               \
        intercepts_block  = edge_count1 + edge_count3;                                 \
        bone_length[vecID-1] += bone_length_block;                                     \
        intercepts[vecID-1] += intercepts_block;                                       \
    }


#define BLOCK_KERNEL_2D(vec, kk, jj, ii) \
    {                                                                              \
        const int vecID = vec;                                                     \
        double bone_length_block = 0.0;                                            \
        int intercepts_block = 0;                                                  \
        /* Init accumulators */                                                    \
        double acc1 = 0.0;                                                         \
        double acc2 = 0.0;                                                         \
        double acc3 = 0.0;                                                         \
        double acc4 = 0.0;                                                         \
                                                                                   \
        unsigned int edge_count1 = 0;                                              \
        unsigned int edge_count2 = 0;                                              \
        unsigned int edge_count3 = 0;                                              \
        unsigned int edge_count4 = 0;                                              \
                                                                                   \
        for (int k = kk + 1; k < kk + BLOCK_SIZE; k += STRIDE * NUM_ACC / 2) {     \
            for (int ij = 0; ij < BLOCK_SIZE; ij += STRIDE) {                      \
                unsigned int i1_prev, i2_prev, j1_prev, j2_prev;                   \
                unsigned int prev_mask1, prev_mask2, prev_mask3, prev_mask4;       \
                                                                                   \
                int i1 = ii + ij;                                                  \
                int j1 = jj;                                                       \
                int i2 = ii;                                                       \
                int j2 = jj + ij;                                                  \
                                                                                   \
                /* Initialise previous mask */                                     \
                LOAD_PREV_2D_POS                                                   \
                                                                                   \
                while (i1 + 1 < ii + BLOCK_SIZE && j2 + 1 < jj + BLOCK_SIZE) {     \
                    double r1, r2, r3, r4;                                         \
                                                                                   \
                    /* Load working set */                                         \
                    LOAD_DATA_SET_2D_POS                                           \
                                                                                   \
                    /* Perform computation */                                      \
                    COMPUTATION                                                    \
                                                                                   \
                    /* Update state of prev_mask */                                \
                    prev_mask1 = curr_mask1;                                       \
                    prev_mask2 = curr_mask2;                                       \
                    prev_mask3 = curr_mask3;                                       \
                    prev_mask4 = curr_mask4;                                       \
                                                                                   \
                    ++i1;                                                          \
                    ++j1;                                                          \
                    ++i2;                                                          \
                    ++j2;                                                          \
                }                                                                  \
            }                                                                      \
        } /* End iteration over dimension k */                                     \
                                                                                   \
        /* Sum up accumulators */                                                  \
        acc1 += acc2;                                                              \
        acc3 += acc4;                                                              \
        edge_count1 += edge_count2;                                                \
        edge_count3 += edge_count4;                                                \
                                                                                   \
        bone_length_block += acc1 + acc3;                                          \
        intercepts_block  += edge_count1 + edge_count3;                            \
        bone_length[vecID-1] += bone_length_block;                                 \
        intercepts[vecID-1]  += intercepts_block;                                  \
    }

#define BLOCK_KERNEL_2D_NEG(vec, kk, jj, ii) \
    {                                                                              \
        const int vecID = vec;                                                     \
        double bone_length_block = 0.0;                                            \
        int intercepts_block = 0;                                                  \
        /* Init accumulators */                                                    \
        double acc1 = 0.0;                                                         \
        double acc2 = 0.0;                                                         \
        double acc3 = 0.0;                                                         \
        double acc4 = 0.0;                                                         \
                                                                                   \
        unsigned int edge_count1 = 0;                                              \
        unsigned int edge_count2 = 0;                                              \
        unsigned int edge_count3 = 0;                                              \
        unsigned int edge_count4 = 0;                                              \
                                                                                   \
        for (int k = kk + 1; k < kk + BLOCK_SIZE; k += STRIDE * NUM_ACC / 2) {     \
            for (int ij = 0; ij < BLOCK_SIZE; ij += STRIDE) {                      \
                unsigned int i1_prev, i2_prev, j1_prev, j2_prev;                   \
                unsigned int prev_mask1, prev_mask2, prev_mask3, prev_mask4;       \
                                                                                   \
                int i1 = ii + (BLOCK_SIZE-1) - ij;                                 \
                int j1 = jj;                                                       \
                int i2 = ii + (BLOCK_SIZE-1);                                      \
                int j2 = jj + ij;                                                  \
                /* Initialise previous mask */                                     \
                LOAD_PREV_2D_NEG                                                   \
                                                                                   \
                while (j2 + 1 < jj + BLOCK_SIZE) {                                  \
                    double r1, r2, r3, r4;                                         \
                                                                                   \
                    /* Load working set */                                         \
                    LOAD_DATA_SET_2D_NEG                                           \
                                                                                   \
                    /* Perform computation */                                      \
                    COMPUTATION                                                    \
                                                                                   \
                    /* Update state of prev_mask */                                \
                    prev_mask1 = curr_mask1;                                       \
                    prev_mask2 = curr_mask2;                                       \
                    prev_mask3 = curr_mask3;                                       \
                    prev_mask4 = curr_mask4;                                       \
                                                                                   \
                    --i1;                                                          \
                    ++j1;                                                          \
                    --i2;                                                          \
                    ++j2;                                                          \
                }                                                                  \
            }                                                                      \
        } /* End iteration over dimension k */                                     \
                                                                                   \
        /* Sum up accumulators */                                                  \
        acc1 += acc2;                                                              \
        acc3 += acc4;                                                              \
        edge_count1 += edge_count2;                                                \
        edge_count3 += edge_count4;                                                \
                                                                                   \
        bone_length_block += acc1 + acc3;                                          \
        intercepts_block  += edge_count1 + edge_count3;                            \
        bone_length[vecID-1] += bone_length_block;                                 \
        intercepts[vecID-1]  += intercepts_block;                                  \
    }

#define BLOCK_KERNEL_3D_POS(vec, kb, jb, ib) \
{                                                                                                                               \
    const int vecID = vec;                                                                                                     \
    double bone_length_block = 0.0;                                                                                             \
    int intercepts_block = 0;                                                                                                   \
                                                                                                                                \
    /* Init accumulators */                                                                                                     \
    double acc1 = 0.0, acc5 = 0.0, acc9  = 0.0;                                                                                 \
    double acc2 = 0.0, acc6 = 0.0, acc10 = 0.0;                                                                                 \
    double acc3 = 0.0, acc7 = 0.0, acc11 = 0.0;                                                                                 \
    double acc4 = 0.0, acc8 = 0.0, acc12 = 0.0;                                                                                 \
                                                                                                                                \
    unsigned int edge_count1 = 0, edge_count5 = 0, edge_count9  = 0;                                                            \
    unsigned int edge_count2 = 0, edge_count6 = 0, edge_count10 = 0;                                                            \
    unsigned int edge_count3 = 0, edge_count7 = 0, edge_count11 = 0;                                                            \
    unsigned int edge_count4 = 0, edge_count8 = 0, edge_count12 = 0;                                                            \
                                                                                                                                \
    unsigned int curr_mask1, curr_mask2,  curr_mask3,  curr_mask4;                                                              \
    unsigned int curr_mask5, curr_mask6,  curr_mask7,  curr_mask8;                                                              \
    unsigned int curr_mask9, curr_mask10, curr_mask11, curr_mask12;                                                             \
                                                                                                                                \
    unsigned int prev_mask1;                                                                                                    \
    unsigned int prev_mask2;                                                                                                    \
    unsigned int prev_mask3;                                                                                                    \
    unsigned int prev_mask4;                                                                                                    \
    unsigned int prev_mask5;                                                                                                    \
    unsigned int prev_mask6;                                                                                                    \
    unsigned int prev_mask7;                                                                                                    \
    unsigned int prev_mask8;                                                                                                    \
    unsigned int prev_mask9;                                                                                                    \
    unsigned int prev_mask10;                                                                                                   \
    unsigned int prev_mask11;                                                                                                   \
    unsigned int prev_mask12;                                                                                                   \
    int prev1;                                                                                                                  \
    int prev2;                                                                                                                  \
    int prev3;                                                                                                                  \
    int prev4;                                                                                                                  \
    int prev5;                                                                                                                  \
    int prev6;                                                                                                                  \
                                                                                                                                \
    double r1, r2,  r3,  r4;                                                                                                    \
    double r5, r6,  r7,  r8;                                                                                                    \
    double r9, r10, r11, r12;                                                                                                   \
                                                                                                                                \
    for (int ks = 2; ks < BLOCK_SIZE; ks += 2*STRIDE) {                                                                         \
        for (int js = ks + 2; js < BLOCK_SIZE; js += STRIDE) {                                                                  \
            int k1 = ks;                                                                                                        \
            int j1 = js;                                                                                                        \
            int k2 = js;                                                                                                        \
            int j2 = ks - 2;                                                                                                    \
            int i = 0;                                                                                                          \
                                                                                                                                \
            prev1 = (ib > 0) ? 1 : 0;                                                                                           \
            prev2 = (ib > 0 && jb+j2 > 0) ? 1 : 0;                                                                              \
            prev_mask1 = hr_sphere_region[ (kb+k1-prev1)*n*n        + (jb+j1-prev1)*n        + (ib+i-prev1)  ]        > 0.5;    \
            prev_mask2 = hr_sphere_region[ (kb+k1+STRIDE-prev1)*n*n + (jb+j1-prev1)*n        + (ib+i-prev1)  ]        > 0.5;    \
            prev_mask3 = hr_sphere_region[ (kb+k2-prev2)*n*n        + (jb+j2-prev2)*n        + (ib+i-prev2)  ]        > 0.5;    \
            prev_mask4 = hr_sphere_region[ (kb+k2-prev1)*n*n        + (jb+j2+STRIDE-prev1)*n + (ib+i-prev1)  ]        > 0.5;    \
                                                                                                                                \
            prev3 = (jb > 0) ? 1 : 0;                                                                                           \
            prev4 = (jb > 0 && kb+j2 > 0) ? 1 : 0;                                                                              \
            prev_mask5 = hr_sphere_region[ (kb+j1-prev3)*n*n        + (jb+i-prev3)*n         + (ib+k1)-prev3 ]        > 0.5;    \
            prev_mask6 = hr_sphere_region[ (kb+j1-prev3)*n*n        + (jb+i-prev3)*n         + (ib+k1+STRIDE)-prev3 ] > 0.5;    \
            prev_mask7 = hr_sphere_region[ (kb+j2-prev4)*n*n        + (jb+i-prev4)*n         + (ib+k2)-prev4 ]        > 0.5;    \
            prev_mask8 = hr_sphere_region[ (kb+j2+STRIDE-prev3)*n*n + (jb+i-prev3)*n         + (ib+k2)-prev3 ]        > 0.5;    \
                                                                                                                                \
            prev5 = (kb > 0) ? 1 : 0;                                                                                           \
            prev6 = (kb > 0 && ib+j2 > 0) ? 1 : 0;                                                                              \
            prev_mask9  = hr_sphere_region[ (kb+i-prev5)*n*n        + (jb+k1-prev5)*n         + (ib+j1-prev5) ]       > 0.5;    \
            prev_mask10 = hr_sphere_region[ (kb+i-prev5)*n*n        + (jb+k1+STRIDE-prev5)*n  + (ib+j1-prev5) ]       > 0.5;    \
            prev_mask11 = hr_sphere_region[ (kb+i-prev6)*n*n        + (jb+k2-prev6)*n         + (ib+j2-prev6) ]       > 0.5;    \
            prev_mask12 = hr_sphere_region[ (kb+i-prev5)*n*n        + (jb+k2-prev5)*n         + (ib+j2+STRIDE-prev5)] > 0.5;    \
                                                                                                                                \
            while (j1 < BLOCK_SIZE) {                                                                                           \
                                                                                                                                \
                r1 = hr_sphere_region[ (kb+k1)*n*n        + (jb+j1)*n        + (ib+i)  ];                                       \
                r2 = hr_sphere_region[ (kb+k1+STRIDE)*n*n + (jb+j1)*n        + (ib+i)  ];                                       \
                r3 = hr_sphere_region[ (kb+k2)*n*n        + (jb+j2)*n        + (ib+i)  ];                                       \
                r4 = hr_sphere_region[ (kb+k2)*n*n        + (jb+j2+STRIDE)*n + (ib+i)  ];                                       \
                                                                                                                                \
                r5 = hr_sphere_region[ (kb+j1)*n*n        + (jb+i)*n         + (ib+k1) ];                                       \
                r6 = hr_sphere_region[ (kb+j1)*n*n        + (jb+i)*n         + (ib+k1+STRIDE) ];                                \
                r7 = hr_sphere_region[ (kb+j2)*n*n        + (jb+i)*n         + (ib+k2) ];                                       \
                r8 = hr_sphere_region[ (kb+j2+STRIDE)*n*n + (jb+i)*n         + (ib+k2)];                                        \
                                                                                                                                \
                r9  = hr_sphere_region[ (kb+i)*n*n        + (jb+k1)*n         + (ib+j1) ]      ;                                \
                r10 = hr_sphere_region[ (kb+i)*n*n        + (jb+k1+STRIDE)*n  + (ib+j1) ]      ;                                \
                r11 = hr_sphere_region[ (kb+i)*n*n        + (jb+k2)*n         + (ib+j2) ]      ;                                \
                r12 = hr_sphere_region[ (kb+i)*n*n        + (jb+k2)*n         + (ib+j2+STRIDE)];                                \
                                                                                                                                \
                acc1 += r1; acc5 += r5; acc9  += r9 ;                                                                           \
                acc2 += r2; acc6 += r6; acc10 += r10;                                                                           \
                acc3 += r3; acc7 += r7; acc11 += r11;                                                                           \
                acc4 += r4; acc8 += r8; acc12 += r12;                                                                           \
                                                                                                                                \
                /* Calculate masks */                                                                                           \
                curr_mask1 = r1 > 0.5; curr_mask5 = r5 > 0.5; curr_mask9  = r9  > 0.5;                                          \
                curr_mask2 = r2 > 0.5; curr_mask6 = r6 > 0.5; curr_mask10 = r10 > 0.5;                                          \
                curr_mask3 = r3 > 0.5; curr_mask7 = r7 > 0.5; curr_mask11 = r11 > 0.5;                                          \
                curr_mask4 = r4 > 0.5; curr_mask8 = r8 > 0.5; curr_mask12 = r12 > 0.5;                                          \
                                                                                                                                \
                /* Detect edge and add to counter */                                                                            \
                edge_count1 += curr_mask1 ^ prev_mask1;                                                                         \
                edge_count2 += curr_mask2 ^ prev_mask2;                                                                         \
                edge_count3 += curr_mask3 ^ prev_mask3;                                                                         \
                edge_count4 += curr_mask4 ^ prev_mask4;                                                                         \
                                                                                                                                \
                edge_count5 += curr_mask5 ^ prev_mask5;                                                                         \
                edge_count6 += curr_mask6 ^ prev_mask6;                                                                         \
                edge_count7 += curr_mask7 ^ prev_mask7;                                                                         \
                edge_count8 += curr_mask8 ^ prev_mask8;                                                                         \
                                                                                                                                \
                edge_count9  += curr_mask9  ^ prev_mask9 ;                                                                      \
                edge_count10 += curr_mask10 ^ prev_mask10;                                                                      \
                edge_count11 += curr_mask11 ^ prev_mask11;                                                                      \
                edge_count12 += curr_mask12 ^ prev_mask12;                                                                      \
                                                                                                                                \
                prev_mask1 = curr_mask1;  prev_mask5 = curr_mask5;  prev_mask9  = curr_mask9 ;                                  \
                prev_mask2 = curr_mask2;  prev_mask6 = curr_mask6;  prev_mask10 = curr_mask10;                                  \
                prev_mask3 = curr_mask3;  prev_mask7 = curr_mask7;  prev_mask11 = curr_mask11;                                  \
                prev_mask4 = curr_mask4;  prev_mask8 = curr_mask8;  prev_mask12 = curr_mask12;                                  \
                                                                                                                                \
                ++k1;                                                                                                           \
                ++j1;                                                                                                           \
                ++k2;                                                                                                           \
                ++j2;                                                                                                           \
                ++i;                                                                                                            \
            }                                                                                                                   \
        }                                                                                                                       \
    }                                                                                                                           \
                                                                                                                                \
    /* Calculate remainder */                                                                                                   \
    for (int jk_s = 2; jk_s < BLOCK_SIZE; jk_s += 2*STRIDE) {                                                                   \
        int k = jk_s;                                                                                                           \
        int j = jk_s;                                                                                                           \
        int i = 0;                                                                                                              \
                                                                                                                                \
        /* Start remainder complete vectors */                                                                                  \
        prev1 = (ib > 0) ? 1 : 0;                                                                                               \
        prev2 = (ib > 0 && jb+j > STRIDE) ? 1 : 0;                                                                              \
        prev_mask1 = hr_sphere_region[ (kb+k-prev1)*n*n        + (jb+j-prev1)*n        + (ib+i-prev1)        ] > 0.5;           \
        prev_mask2 = hr_sphere_region[ (kb+k-prev2)*n*n        + (jb+j-STRIDE-prev2)*n + (ib+i-prev2)        ] > 0.5;           \
                                                                                                                                \
        prev3 = (jb > 0) ? 1 : 0;                                                                                               \
        prev4 = (jb > 0 && kb+k > STRIDE) ? 1 : 0;                                                                              \
        prev_mask5 = hr_sphere_region[ (kb+k-prev3)*n*n        + (jb+i-prev3)*n        + (ib+j-prev3)        ] > 0.5;           \
        prev_mask6 = hr_sphere_region[ (kb+k-STRIDE-prev4)*n*n + (jb+i-prev4)*n        + (ib+j-prev4)        ] > 0.5;           \
                                                                                                                                \
        prev5 = (kb > 0) ? 1 : 0;                                                                                               \
        prev6 = (kb > 0 && ib+j > STRIDE) ? 1 : 0;                                                                              \
        prev_mask9  = hr_sphere_region[ (kb+i-prev5)*n*n       + (jb+k-prev5)*n        + (ib+j-prev5)        ] > 0.5;           \
        prev_mask10 = hr_sphere_region[ (kb+i-prev6)*n*n       + (jb+k-prev6)*n        + (ib+j-STRIDE-prev6) ] > 0.5;           \
                                                                                                                                \
        while (j < BLOCK_SIZE) {                                                                                                \
            r1 = hr_sphere_region[ (kb+k)*n*n        + (jb+j)*n        + (ib+i)        ];                                       \
            r2 = hr_sphere_region[ (kb+k)*n*n        + (jb+j-STRIDE)*n + (ib+i)        ];                                       \
                                                                                                                                \
            r5 = hr_sphere_region[ (kb+k)*n*n        + (jb+i)*n        + (ib+j)        ];                                       \
            r6 = hr_sphere_region[ (kb+k-STRIDE)*n*n + (jb+i)*n        + (ib+j)        ];                                       \
                                                                                                                                \
            r9  = hr_sphere_region[ (kb+i)*n*n       + (jb+k)*n        + (ib+j)        ];                                       \
            r10 = hr_sphere_region[ (kb+i)*n*n       + (jb+k)*n        + (ib+j-STRIDE) ];                                       \
                                                                                                                                \
            acc1 += r1; acc5 += r5; acc9  += r9 ;                                                                               \
            acc2 += r2; acc6 += r6; acc10 += r10;                                                                               \
                                                                                                                                \
            /* Calculate masks */                                                                                               \
            curr_mask1 = r1 > 0.5; curr_mask5 = r5 > 0.5; curr_mask9  = r9  > 0.5;                                              \
            curr_mask2 = r2 > 0.5; curr_mask6 = r6 > 0.5; curr_mask10 = r10 > 0.5;                                              \
                                                                                                                                \
            /* Detect edge and add to counter */                                                                                \
            edge_count1 += curr_mask1 ^ prev_mask1;                                                                             \
            edge_count2 += curr_mask2 ^ prev_mask2;                                                                             \
                                                                                                                                \
            edge_count5 += curr_mask5 ^ prev_mask5;                                                                             \
            edge_count6 += curr_mask6 ^ prev_mask6;                                                                             \
                                                                                                                                \
            edge_count9  += curr_mask9  ^ prev_mask9 ;                                                                          \
            edge_count10 += curr_mask10 ^ prev_mask10;                                                                          \
                                                                                                                                \
            prev_mask1 = curr_mask1;  prev_mask5 = curr_mask5;  prev_mask9  = curr_mask9 ;                                      \
            prev_mask2 = curr_mask2;  prev_mask6 = curr_mask6;  prev_mask10 = curr_mask10;                                      \
                                                                                                                                \
            ++k;                                                                                                                \
            ++j;                                                                                                                \
            ++i;                                                                                                                \
        }                                                                                                                       \
                                                                                                                                \
    }                                                                                                                           \
                                                                                                                                \
    /* Calculate diagonal */                                                                                                    \
    prev1 = (kb >0 && jb > 0 && ib > 0) ? 1 : 0;                                                                                \
    prev_mask1 = hr_sphere_region[ (kb-prev1)*n*n + (jb-prev1)*n + (ib-prev1) ] > 0.5;                                          \
    for (int i = 0; i < BLOCK_SIZE; ++i) {                                                                                      \
        r1 = hr_sphere_region[ (kb+i)*n*n + (jb+i)*n + (ib+i)];                                                                 \
        acc1 += r1;                                                                                                             \
        curr_mask1 = r1 > 0.5;                                                                                                  \
        edge_count1 += curr_mask1 ^ prev_mask1;                                                                                 \
        prev_mask1 = curr_mask1;                                                                                                \
    }                                                                                                                           \
                                                                                                                                \
    /* Sum up accumulators */                                                                                                   \
    acc1 += acc2;                                                                                                               \
    acc3 += acc4;                                                                                                               \
    acc5 += acc6;                                                                                                               \
    acc7 += acc8;                                                                                                               \
    acc9 += acc10;                                                                                                              \
    acc11 += acc12;                                                                                                             \
    edge_count1  += edge_count2;                                                                                                \
    edge_count3  += edge_count4;                                                                                                \
    edge_count5  += edge_count6;                                                                                                \
    edge_count7  += edge_count8;                                                                                                \
    edge_count9  += edge_count10;                                                                                               \
    edge_count11 += edge_count12;                                                                                               \
                                                                                                                                \
    bone_length_block += acc1 + acc3 + acc5 + acc7 + acc9 + acc11;                                                              \
    intercepts_block  += edge_count1 + edge_count3 + edge_count5 + edge_count7 + edge_count9 + edge_count11;                    \
    bone_length[vecID-1] += bone_length_block;                                                                                  \
    intercepts[vecID-1]  += intercepts_block;                                                                                   \
}


/* Function declarations*/
double mil_1D(const double *hr_sphere_region, int* intercepts, int n, const int kk, const int jj, const int ii,  const int vecID);
double mil_2D_pos(const double *hr_sphere_region, int* intercepts, int n, const int kk, const int jj, const int ii,  const int vecID);
double mil_2D_neg(const double *hr_sphere_region, int* intercepts, int n, const int kk, const int jj, const int ii,  const int vecID);

#endif //BONEMAP_MIL2_H

