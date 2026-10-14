/* Test-only bridge to the production C ACA function; not a solver ABI.
 * Calls are serial because the HACApK callback has no context argument.
 */
#include "cHACApK_base.h"
#include <stdlib.h>

static const double *input;
static int rows, columns, calls;

double cHACApK_entry_ij(int i, int j, int problem)
{
    (void)problem;
    ++calls;
    return input[(i - 1) * columns + j - rows - 1];
}

__declspec(dllexport) int aca_probe(const double *matrix, int m, int n,
                                  double eps, double *output, int *entries)
{
    int kmax = m < n ? m : n;
    double param[100] = {0};
    double *left, *right;
    int *lod;
    int rank, i, j, k;
    if (m < 1 || n < 1 || !matrix || !output || !entries) return -1;
    left = (double *)calloc((size_t)m * kmax, sizeof(double));
    right = (double *)calloc((size_t)n * kmax, sizeof(double));
    lod = (int *)calloc((size_t)m + n + 1, sizeof(int));
    if (!left || !right || !lod) {
        free(left); free(right); free(lod);
        return -1;
    }
    for (i = 1; i <= m + n; ++i) lod[i] = i;
    input = matrix; rows = m; columns = n; calls = 0;
    param[61] = 1; param[64] = 1;
    rank = cHACApK_aca(left, right, param, m, n, 1, m + 1, lod,
                      0, kmax, eps, 1.0, eps * 1e-3);
    for (i = 0; i < m; ++i) {
        for (j = 0; j < n; ++j) {
            double value = 0.0;
            for (k = 0; k < rank; ++k)
                value += left[i + m * k] * right[j + n * k];
            output[i * n + j] = value;
        }
    }
    *entries = calls;
    free(left); free(right); free(lod);
    input = NULL;
    return rank;
}
