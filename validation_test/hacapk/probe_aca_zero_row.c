/* Call the production ACA function on a matrix with a hidden residual row.
 * Link cHACApK_base.c and cHACApK_lib.c with function-level dead stripping.
 * This diagnostic exits nonzero when an inaccurate early stop is reproduced.
 */
#include "cHACApK_base.h"
#include <math.h>
#include <stdio.h>

static const double matrix[3][3] = {
    {1.0, 1.0, 1.0},
    {1.0, 1.0, 2.0},
    {1.0, 1.0, 1.0}
};

double cHACApK_entry_ij(int i, int j, int problem)
{
    (void)problem;
    return matrix[i - 1][j - 1];
}

int main(void)
{
    double left[9] = {0}, right[9] = {0}, param[100] = {0};
    int lod[4] = {0, 1, 2, 3};
    double error2 = 0.0, norm2 = 0.0;
    int i, j, k, rank;
    param[61] = 1;
    param[64] = 1;
    rank = cHACApK_aca(left, right, param, 3, 3, 1, 1, lod,
                       0, 3, 1e-10, 1.0, 1e-13);
    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            double value = 0.0, residual;
            for (k = 0; k < rank; ++k)
                value += left[i + 3 * k] * right[j + 3 * k];
            residual = value - matrix[i][j];
            error2 += residual * residual;
            norm2 += matrix[i][j] * matrix[i][j];
        }
    }
    printf("rank=%d relative_frobenius_error=%.17g\n", rank, sqrt(error2 / norm2));
    return sqrt(error2 / norm2) < 1e-10 ? 0 : 1;
}
