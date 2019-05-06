#include <stdio.h>
#include <math.h>

#define TOLERANCE 1e-5

void print_vec (double* v, int n) {
    for (int i = 0; i < n; ++i) {
        printf("%f\n", v[i]);
    }
}

void print_matrix (double m[3][3], int n1, int n2) {
    for (int i = 0; i < n1; ++i) {
        for (int j = 0; j < n2; ++j) {
            printf("%f ", m[i][j]);
        }
        printf("\n");
    }
}

///
/// \param v
/// \return 0: if zero vector found. 1: otherwise
/// flop count: 3 adds + 3 mults + 1 sqrt + 3 divs (if non-zero vector)
int normalize(double v[3]) {
    //
    // Determine magnitud
    //
    double m = 0.0;
    for (int i = 0; i < 3; ++i) {
        m += v[i] * v[i];
    }
    m = sqrt(m);

    if (m < TOLERANCE) {
        for (int i = 0; i < 3; ++i) {
            v[i] = 0.0;
        }
        return 0;
    }
    else {
        for (int i = 0; i < 3; ++i) {
            v[i] /= m;
        }
        return 1;
    }
}

// Smith, Oliver K. (April 1961),
// "Eigenvalues of a symmetric 3 × 3 matrix.",
// Communications of the ACM, 4 (4): 168,
// doi:10.1145/355578.366316
///
/// \param M: Input matrix
/// \param eVec: eigen vectors
/// \param eVal: eigen values
/// flop count (best-case) : 12  mults + 14  adds (assumming diagonal matrix)
/// flop count (worst-case): 132 mults + 127 adds + 16 divs + 3 sqrts + 1*Cost_acos + 2*Cost_cos
void eigen3(double M[3][3], double eVec[3][3], double eVal[3]) {
    //
    // Clean eVal and create 3x3 identity matrix
    //
    double I[3][3];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            I[i][j] = 0.0;
            eVec[i][j] = 0.0;
        }
        I[i][i] = 1.0;
    }

    double p1 = M[0][1]*M[0][1] + M[0][2]*M[0][2] + M[1][2]*M[1][2];
    if (p1 == 0) {
        //
        // M is a diagonal matrix
        //
        eVal[0] = M[0][0];
        eVal[1] = M[1][1];
        eVal[2] = M[2][2];
    }
    else {
        double diagMean = ( M[0][0] + M[1][1] + M[2][2] ) / 3;
        double a1 = M[0][0] - diagMean;
        double a2 = M[1][1] - diagMean;
        double a3 = M[2][2] - diagMean;
        double p2 = a1*a1 + a2*a2 + a3*a3 + 2*p1;
        double p = sqrt(p2 / 6);
        //
        // Create matrix B
        //
        double B[3][3];
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                B[i][j] = (1 / p) * (M[i][j] - diagMean * I[i][j]);
            }
        }
        //
        // Determinant of matrix B
        //
        double r = B[0][0]*B[1][1]*B[2][2] + B[0][1]*B[1][2]*B[2][0] + B[0][2]*B[1][0]*B[2][1] \
                 - B[0][2]*B[1][1]*B[2][0] - B[0][1]*B[1][0]*B[2][2] - B[0][0]*B[1][2]*B[2][1];
        r /= 2;
        double phi;
        if (r <= -1)
            phi = M_PI / 3;
        else if (r >= 1)
            phi = 0;
        else
            phi = acos(r) / 3; // Determine flop count of acos function
        //
        // The eigenvalues satisfy eig3 <= eig2 <= eig1
        //
        eVal[0] = diagMean + 2 * p * cos(phi);
        eVal[2] = diagMean + 2 * p * cos(phi + (2*M_PI/3));
        eVal[1] = 3 * diagMean - eVal[0] - eVal[2];     // since trace(A) = eig1 + eig2 + eig3
    }

    //
    // Find eigenvectors given the eigen values
    //
    double A[3][3][3];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            A[0][i][j] = M[i][j];
            A[1][i][j] = M[i][j];
            A[2][i][j] = M[i][j];
        }
        A[0][i][i] -= eVal[0];
        A[1][i][i] -= eVal[1];
        A[2][i][i] -= eVal[2];
    }

    //
    // flop count (worst case): 27*(3 adds + 3 mults) + 3*(3 adds + 3 mults + 1 sqrt) + 3 divs
    //
    for (int e = 0; e < 3; ++e) {
        //
        // Matrix matrix multiplication till non zero column is found.
        // The non zero column is the eigen-vector
        //
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                for (int k = 0; k < 3; ++k) {
                    int m1 = (e + 1) % 3;
                    int m2 = (e + 2) % 3;
                    eVec[e][j] += A[m1][j][k] * A[m2][k][i];
                }
            }

            if ( normalize(eVec[e]) ) {
                break;
            }
        }
    }

}

void eigen_test () {
    double M[3][3] = { {1, 3, 6} ,
                       {3,-5,-6} ,
                       {6, -6, 4} };
    double eVal[3];
    double eVec[3][3];

    eigen3(M, eVec, eVal);

    print_vec(eVal, 3);
    print_matrix( eVec, 3, 3);
}