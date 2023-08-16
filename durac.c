//
// Created by Corey Richardson on 8/15/23.
//

#include "durac.h"
#include <stdlib.h>
#include <complex.h>
#include <math.h>

#define N 12 // Number of spatial points in each dimension
#define dx 0.12 // Spatial step size
#define dt 0.06 // Time step size
#define m 17.6 // Mass of the particle
#define M 7.7 // Mass of the singularity

// Simplified Dirac equation in 4+1 dimensions with Schwarzschild-like metric
double long solve_dirac_5D() {
    double complex psi[N][N][N][N][4]; // 4-component spinor in 5D
    double g_tt[N][N][N][N], g_xx[N][N][N][N]; // Metric components

    // Initialize wave function with random values
    for (int a = 0; a < N; a++)
        for (int b = 0; b < N; b++)
            for (int c = 0; c < N; c++)
                for (int d = 0; d < N; d++) {
                    psi[a][b][c][d][0] = rand() / (double)RAND_MAX;
                    psi[a][b][c][d][1] = rand() / (double)RAND_MAX;
                    psi[a][b][c][d][2] = rand() / (double)RAND_MAX;
                    psi[a][b][c][d][3] = rand() / (double)RAND_MAX;
                }

    // Simplified Schwarzschild-like metric in 5D
    for (int a = 0; a < N; a++)
        for (int b = 0; b < N; b++)
            for (int c = 0; c < N; c++)
                for (int d = 0; d < N; d++) {
                    double r = sqrt(a * a + b * b + c * c + d * d) * dx;
                    g_tt[a][b][c][d] = -1.0 + 2.0 * M / r;
                    g_xx[a][b][c][d] = 1.0 / (1.0 - 2.0 * M / r);
                }

    // Time evolution
    for (int t = 0; t < 1000; t++) {
        // Spatial loops
        for (int a = 1; a < N - 1; a++)
            for (int b = 1; b < N - 1; b++)
                for (int c = 1; c < N - 1; c++)
                    for (int d = 1; d < N - 1; d++) {
                        // Dirac's equation with metric components
                        for (int k = 0; k < 4; k++) {
                            double complex d_psi_dx = (psi[a + 1][b][c][d][k] - psi[a - 1][b][c][d][k]) / (2 * dx);
                            psi[a][b][c][d][k] += dt * (-I * g_tt[a][b][c][d] * d_psi_dx + m * g_xx[a][b][c][d] * psi[a][b][c][d][k]);
                        }
                    }
    }

    // Output some results
    for (int a = 0; a < N; a++)
        for (int b = 0; b < N; b++)
            for (int c = 0; c < N; c++)
                for (int d = 0; d < N; d++) {
                    if(a == 7 && b == 7 && c == 7 && d == 7) {
                        return fabsl(sqrtl(creal(psi[0][0][0][0][0])));
                    }
                }
}
