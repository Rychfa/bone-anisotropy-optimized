#ifndef BONEMAP_MIL2_H
#define BONEMAP_MIL2_H

#include <immintrin.h>
#include <stdio.h>

#define BLOCK_SIZE 8
#define NUM_ACC 4
#define STRIDE 2

#define LOAD_PREV_1D                                                                     \
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

#define LOAD_DATA_SET_2D_POS                                            \
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

#define LOAD_DATA_SET_2D_POS_UNEVEN                                             \
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

#define LOAD_DATA_SET_2D_NEG                                              \
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


#define BLOCK_KERNEL_1D(vec, kk, jj, ii)                                           \
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

#define BLOCK_KERNEL_3D_POS(vec, kk, jj, ii)                                       \
    {                                                                              \
        const int vecID = vec;                                                     \
        double bone_length_block = 0.0;                                            \
        int intercepts_block = 0;                                                  \
        /* Init accumulators */                                                    \
        double acc1 = 0.0, acc5 = 0.0, acc9  = 0.0;                                \
        double acc2 = 0.0, acc6 = 0.0, acc10 = 0.0;                                \
        double acc3 = 0.0, acc7 = 0.0, acc11 = 0.0;                                \
        double acc4 = 0.0, acc8 = 0.0, acc12 = 0.0;                                \
                                                                                   \
        unsigned int edge_count1 = 0;                                              \
        unsigned int edge_count2 = 0;                                              \
        unsigned int edge_count3 = 0;                                              \
        unsigned int edge_count4 = 0;                                              \
                                                                                   \
        for (int k_start = kk + 1; k_start < kk + BLOCK_SIZE; k_start += STRIDE * 2) {               \
            for (int j_start = jj + 1; j_start < jj + BLOCK_SIZE; j_start += STRIDE * 2) {           \
                unsigned int i1_prev, i2_prev, j1_prev, j2_prev;                   \
                unsigned int prev_mask1, prev_mask2, prev_mask3, prev_mask4;       \
                                                                                   \
                int i = ii;                                                        \
                int j = j_start;                                                   \
                int k = k_start;                                                   \
                LOAD_PREV_3D_POS                                                   \
                                                                                   \
                while (j + 1 < j + BLOCK_SIZE && k + 1 < k + BLOCK_SIZE) {         \
                    double r1, r2, r3, r4;                                         \
                                                                                   \
                    /* Load working set */                                         \
                    LOAD_DATA_SET_3D                                               \
                                                                                   \
                    /* Perform computation */                                      \
                    COMPUTATION_3D                                                 \
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


/* Function declarations*/
double mil_1D(const double *hr_sphere_region, int* intercepts, int n, const int kk, const int jj, const int ii,  const int vecID);
double mil_2D_pos(const double *hr_sphere_region, int* intercepts, int n, const int kk, const int jj, const int ii,  const int vecID);
double mil_2D_neg(const double *hr_sphere_region, int* intercepts, int n, const int kk, const int jj, const int ii,  const int vecID);

#endif //BONEMAP_MIL2_H

