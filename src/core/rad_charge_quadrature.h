// Exact charge-Gram quadrature shared by CPU and mexcuda host setup.
#pragma once
#include <vector>
static void rad_gl4_duffy_tet(std::vector<double>& pts, std::vector<double>& w)
{
    static const double GL4x[4] = {0.06943184420297371, 0.33000947820757187,
                                   0.66999052179242813, 0.93056815579702629};   // 4-pt Gauss-Legendre on [0,1]
    static const double GL4w[4] = {0.17392742256872693, 0.32607257743127307,
                                   0.32607257743127307, 0.17392742256872693};
    pts.clear(); w.clear();
    pts.reserve(64 * 3); w.reserve(64);
    for (int ia = 0; ia < 4; ++ia)
        for (int ib = 0; ib < 4; ++ib)
            for (int ic = 0; ic < 4; ++ic) {
                const double aa = GL4x[ia], bb = GL4x[ib], cc = GL4x[ic];
                pts.push_back(aa);
                pts.push_back(bb * (1.0 - aa));
                pts.push_back(cc * (1.0 - aa) * (1.0 - bb));
                w.push_back(GL4w[ia] * GL4w[ib] * GL4w[ic] * (1.0 - aa) * (1.0 - aa) * (1.0 - bb));
            }
}

// Built-in Dunavant degree-5 symmetric triangle rule (7 nodes; bary (l1,l2,l3) + weight, weights sum to 1).
static const double RAD_DUN5[7][4] = {
    {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0, 0.225},
    {0.0597158717, 0.4701420641, 0.4701420641, 0.1323941527},
    {0.4701420641, 0.0597158717, 0.4701420641, 0.1323941527},
    {0.4701420641, 0.4701420641, 0.0597158717, 0.1323941527},
    {0.7974269853, 0.1012865073, 0.1012865073, 0.1259391805},
    {0.1012865073, 0.7974269853, 0.1012865073, 0.1259391805},
    {0.1012865073, 0.1012865073, 0.7974269853, 0.1259391805},
};
