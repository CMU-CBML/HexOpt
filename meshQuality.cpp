#include "meshQuality.h"

void lapSmth(std::vector<double>& x, const std::vector<std::vector<int>>& adjPts, const int& pNum, const std::vector<int>& pNumM3) {
#pragma omp parallel for
    for (int i = 0; i < pNum; ++i) {
        double dist, maxDist = 0, minDist = DBL_MAX;
        int minDistIdx;
        for (const auto& j : adjPts[i]) {
            dist = (x[pNumM3[i]] - x[j]) * (x[pNumM3[i]] - x[j]) + (x[pNumM3[i] + 1] - x[j + 1]) * (x[pNumM3[i] + 1] - x[j + 1]) + (x[pNumM3[i] + 2] - x[j + 2]) * (x[pNumM3[i] + 2] - x[j + 2]);
            maxDist = std::max(maxDist, dist);
            if (dist < minDist) {
                minDist = dist;
                minDistIdx = j;
            }
        }
        if (minDist < 0.01 * maxDist && minDist > 0) {
            x[pNumM3[i]] += sqrt(maxDist / (1000000000 * minDist)) * (x[pNumM3[i]] - x[minDistIdx]);
            x[pNumM3[i] + 1] += sqrt(maxDist / (1000000000 * minDist)) * (x[pNumM3[i] + 1] - x[minDistIdx + 1]);
            x[pNumM3[i] + 2] += sqrt(maxDist / (1000000000 * minDist)) * (x[pNumM3[i] + 2] - x[minDistIdx + 2]);
        }
    }
}

void sJGrad(const double& x0, const double& y0, const double& z0,
    const double& x1, const double& y1, const double& z1,
    const double& x2, const double& y2, const double& z2,
    const double& x3, const double& y3, const double& z3,
    const double& x4, const double& y4, const double& z4,
    const double& x5, const double& y5, const double& z5,
    const double& x6, const double& y6, const double& z6,
    const double& x7, const double& y7, const double& z7,
    double& dx0, double& dy0, double& dz0,
    double& dx1, double& dy1, double& dz1,
    double& dx2, double& dy2, double& dz2,
    double& dx3, double& dy3, double& dz3,
    double& dx4, double& dy4, double& dz4,
    double& dx5, double& dy5, double& dz5,
    double& dx6, double& dy6, double& dz6,
    double& dx7, double& dy7, double& dz7,
    const double& sJThres, bool& eFSJ) {
    // faceCenter
    double vec[3][3] = {
        {x1 + x2 + x5 + x6 - x0 - x3 - x4 - x7, y1 + y2 + y5 + y6 - y0 - y3 - y4 - y7, z1 + z2 + z5 + z6 - z0 - z3 - z4 - z7},
        {x2 + x3 + x6 + x7 - x0 - x1 - x4 - x5, y2 + y3 + y6 + y7 - y0 - y1 - y4 - y5, z2 + z3 + z6 + z7 - z0 - z1 - z4 - z5},
        {x4 + x5 + x6 + x7 - x0 - x1 - x2 - x3, y4 + y5 + y6 + y7 - y0 - y1 - y2 - y3, z4 + z5 + z6 + z7 - z0 - z1 - z2 - z3}
    };
    double vol = vec[0][0] * (vec[2][1] * vec[1][2] - vec[1][1] * vec[2][2]) +
        vec[0][1] * (vec[2][2] * vec[1][0] - vec[1][2] * vec[2][0]) +
        vec[0][2] * (vec[2][0] * vec[1][1] - vec[1][0] * vec[2][1]);
    double len3 = sqrt(DOT(vec[0], vec[0]) * DOT(vec[1], vec[1]) * DOT(vec[2], vec[2]));
    double sz = (sqrt((x0 - x1) * (x0 - x1) + (y0 - y1) * (y0 - y1) + (z0 - z1) * (z0 - z1)) +
                 sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2)) +
                 sqrt((x2 - x3) * (x2 - x3) + (y2 - y3) * (y2 - y3) + (z2 - z3) * (z2 - z3)) +
                 sqrt((x3 - x0) * (x3 - x0) + (y3 - y0) * (y3 - y0) + (z3 - z0) * (z3 - z0)) +
                 sqrt((x0 - x4) * (x0 - x4) + (y0 - y4) * (y0 - y4) + (z0 - z4) * (z0 - z4)) +
                 sqrt((x1 - x5) * (x1 - x5) + (y1 - y5) * (y1 - y5) + (z1 - z5) * (z1 - z5)) +
                 sqrt((x2 - x6) * (x2 - x6) + (y2 - y6) * (y2 - y6) + (z2 - z6) * (z2 - z6)) +
                 sqrt((x3 - x7) * (x3 - x7) + (y3 - y7) * (y3 - y7) + (z3 - z7) * (z3 - z7)) +
                 sqrt((x4 - x5) * (x4 - x5) + (y4 - y5) * (y4 - y5) + (z4 - z5) * (z4 - z5)) +
                 sqrt((x5 - x6) * (x5 - x6) + (y5 - y6) * (y5 - y6) + (z5 - z6) * (z5 - z6)) +
                 sqrt((x6 - x7) * (x6 - x7) + (y6 - y7) * (y6 - y7) + (z6 - z7) * (z6 - z7)) +
                 sqrt((x7 - x4) * (x7 - x4) + (y7 - y4) * (y7 - y4) + (z7 - z4) * (z7 - z4))) / 12.0, sz2 = 1.0 / sz; sz *= sz;
    if (vol >= 0) {
        eFSJ = false;
        sz2 /= 16.0;
        dx0 = sz2 * (y4 * z1 + y5 * z1 + y1 * z2 + y5 * z2 - y7 * z2 + y1 * z3 - y4 * z3 - y7 * z3 -
            y1 * z4 - y5 * z4 + y7 * z4 - y1 * z5 + y4 * z5 + y7 * z5 - y4 * z7 - y5 * z7 +
            y3 * (-z1 - z2 + z4 + z7) + y2 * (-z1 + z3 - z5 + z7));
        dy0 = -sz2 * (x4 * z1 + x5 * z1 + x1 * z2 + x5 * z2 - x7 * z2 + x1 * z3 - x4 * z3 - x7 * z3 -
            x1 * z4 - x5 * z4 + x7 * z4 - x1 * z5 + x4 * z5 + x7 * z5 - x4 * z7 - x5 * z7 +
            x3 * (-z1 - z2 + z4 + z7) + x2 * (-z1 + z3 - z5 + z7));
        dz0 = sz2 * (x4 * y1 + x5 * y1 + x1 * y2 + x5 * y2 - x7 * y2 + x1 * y3 - x4 * y3 - x7 * y3 -
            x1 * y4 - x5 * y4 + x7 * y4 - x1 * y5 + x4 * y5 + x7 * y5 - x4 * y7 - x5 * y7 +
            x3 * (-y1 - y2 + y4 + y7) + x2 * (-y1 + y3 - y5 + y7));
        dx1 = -sz2 * (y4 * z0 + y5 * z0 + y0 * z2 - y5 * z2 - y6 * z2 + y0 * z3 + y4 * z3 - y6 * z3 -
            y0 * z4 + y5 * z4 + y6 * z4 - y0 * z5 - y4 * z5 + y6 * z5 - y4 * z6 - y5 * z6 +
            y3 * (-z0 + z2 - z4 + z6) + y2 * (-z0 - z3 + z5 + z6));
        dy1 = sz2 * (x4 * z0 + x5 * z0 + x0 * z2 - x5 * z2 - x6 * z2 + x0 * z3 + x4 * z3 - x6 * z3 -
            x0 * z4 + x5 * z4 + x6 * z4 - x0 * z5 - x4 * z5 + x6 * z5 - x4 * z6 - x5 * z6 +
            x3 * (-z0 + z2 - z4 + z6) + x2 * (-z0 - z3 + z5 + z6));
        dz1 = -sz2 * (x4 * y0 + x5 * y0 + x0 * y2 - x5 * y2 - x6 * y2 + x0 * y3 + x4 * y3 - x6 * y3 -
            x0 * y4 + x5 * y4 + x6 * y4 - x0 * y5 - x4 * y5 + x6 * y5 - x4 * y6 - x5 * y6 +
            x3 * (-y0 + y2 - y4 + y6) + x2 * (-y0 - y3 + y5 + y6));
        dx2 = sz2 * (-y5 * z0 + y7 * z0 + y0 * z1 - y5 * z1 - y6 * z1 - y0 * z3 + y6 * z3 + y7 * z3 +
            y0 * z5 - y6 * z5 - y7 * z5 + y5 * z6 - y7 * z6 + y1 * (-z0 - z3 + z5 + z6) +
            y3 * (z0 + z1 - z6 - z7) - y0 * z7 + y5 * z7 + y6 * z7);
        dy2 = -sz2 * (-x5 * z0 + x7 * z0 + x0 * z1 - x5 * z1 - x6 * z1 - x0 * z3 + x6 * z3 + x7 * z3 +
            x0 * z5 - x6 * z5 - x7 * z5 + x5 * z6 - x7 * z6 + x1 * (-z0 - z3 + z5 + z6) +
            x3 * (z0 + z1 - z6 - z7) - x0 * z7 + x5 * z7 + x6 * z7);
        dz2 = sz2 * (-x5 * y0 + x7 * y0 + x0 * y1 - x5 * y1 - x6 * y1 - x0 * y3 + x6 * y3 + x7 * y3 +
            x0 * y5 - x6 * y5 - x7 * y5 + x5 * y6 - x7 * y6 + x1 * (-y0 - y3 + y5 + y6) +
            x3 * (y0 + y1 - y6 - y7) - x0 * y7 + x5 * y7 + x6 * y7);
        dx3 = sz2 * (y4 * z0 + y7 * z0 + y0 * z1 + y4 * z1 - y6 * z1 + y0 * z2 - y6 * z2 - y7 * z2 -
            y0 * z4 + y6 * z4 + y7 * z4 - y4 * z6 - y7 * z6 + y1 * (-z0 + z2 - z4 + z6) -
            y0 * z7 - y4 * z7 + y6 * z7 + y2 * (-z0 - z1 + z6 + z7));
        dy3 = -sz2 * (x4 * z0 + x7 * z0 + x0 * z1 + x4 * z1 - x6 * z1 + x0 * z2 - x6 * z2 - x7 * z2 -
            x0 * z4 + x6 * z4 + x7 * z4 - x4 * z6 - x7 * z6 + x1 * (-z0 + z2 - z4 + z6) -
            x0 * z7 - x4 * z7 + x6 * z7 + x2 * (-z0 - z1 + z6 + z7));
        dz3 = sz2 * (x4 * y0 + x7 * y0 + x0 * y1 + x4 * y1 - x6 * y1 + x0 * y2 - x6 * y2 - x7 * y2 -
            x0 * y4 + x6 * y4 + x7 * y4 - x4 * y6 - x7 * y6 + x1 * (-y0 + y2 - y4 + y6) -
            x0 * y7 - x4 * y7 + x6 * y7 + x2 * (-y0 - y1 + y6 + y7));
        dx4 = -sz2 * (-y5 * z0 + y7 * z0 + y0 * z1 - y5 * z1 - y6 * z1 - y0 * z3 + y6 * z3 + y7 * z3 +
            y0 * z5 - y6 * z5 - y7 * z5 + y5 * z6 - y7 * z6 + y1 * (-z0 - z3 + z5 + z6) +
            y3 * (z0 + z1 - z6 - z7) - y0 * z7 + y5 * z7 + y6 * z7);
        dy4 = sz2 * (-x5 * z0 + x7 * z0 + x0 * z1 - x5 * z1 - x6 * z1 - x0 * z3 + x6 * z3 + x7 * z3 +
            x0 * z5 - x6 * z5 - x7 * z5 + x5 * z6 - x7 * z6 + x1 * (-z0 - z3 + z5 + z6) +
            x3 * (z0 + z1 - z6 - z7) - x0 * z7 + x5 * z7 + x6 * z7);
        dz4 = -sz2 * (-x5 * y0 + x7 * y0 + x0 * y1 - x5 * y1 - x6 * y1 - x0 * y3 + x6 * y3 + x7 * y3 +
            x0 * y5 - x6 * y5 - x7 * y5 + x5 * y6 - x7 * y6 + x1 * (-y0 - y3 + y5 + y6) +
            x3 * (y0 + y1 - y6 - y7) - x0 * y7 + x5 * y7 + x6 * y7);
        dx5 = -sz2 * (y4 * z0 + y7 * z0 + y0 * z1 + y4 * z1 - y6 * z1 + y0 * z2 - y6 * z2 - y7 * z2 -
            y0 * z4 + y6 * z4 + y7 * z4 - y4 * z6 - y7 * z6 + y1 * (-z0 + z2 - z4 + z6) -
            y0 * z7 - y4 * z7 + y6 * z7 + y2 * (-z0 - z1 + z6 + z7));
        dy5 = sz2 * (x4 * z0 + x7 * z0 + x0 * z1 + x4 * z1 - x6 * z1 + x0 * z2 - x6 * z2 - x7 * z2 -
            x0 * z4 + x6 * z4 + x7 * z4 - x4 * z6 - x7 * z6 + x1 * (-z0 + z2 - z4 + z6) -
            x0 * z7 - x4 * z7 + x6 * z7 + x2 * (-z0 - z1 + z6 + z7));
        dz5 = -sz2 * (x4 * y0 + x7 * y0 + x0 * y1 + x4 * y1 - x6 * y1 + x0 * y2 - x6 * y2 - x7 * y2 -
            x0 * y4 + x6 * y4 + x7 * y4 - x4 * y6 - x7 * y6 + x1 * (-y0 + y2 - y4 + y6) -
            x0 * y7 - x4 * y7 + x6 * y7 + x2 * (-y0 - y1 + y6 + y7));
        dx6 = -sz2 * (y4 * z1 + y5 * z1 + y1 * z2 + y5 * z2 - y7 * z2 + y1 * z3 - y4 * z3 - y7 * z3 -
            y1 * z4 - y5 * z4 + y7 * z4 - y1 * z5 + y4 * z5 + y7 * z5 - y4 * z7 - y5 * z7 +
            y3 * (-z1 - z2 + z4 + z7) + y2 * (-z1 + z3 - z5 + z7));
        dy6 = sz2 * (x4 * z1 + x5 * z1 + x1 * z2 + x5 * z2 - x7 * z2 + x1 * z3 - x4 * z3 - x7 * z3 -
            x1 * z4 - x5 * z4 + x7 * z4 - x1 * z5 + x4 * z5 + x7 * z5 - x4 * z7 - x5 * z7 +
            x3 * (-z1 - z2 + z4 + z7) + x2 * (-z1 + z3 - z5 + z7));
        dz6 = -sz2 * (x4 * y1 + x5 * y1 + x1 * y2 + x5 * y2 - x7 * y2 + x1 * y3 - x4 * y3 - x7 * y3 -
            x1 * y4 - x5 * y4 + x7 * y4 - x1 * y5 + x4 * y5 + x7 * y5 - x4 * y7 - x5 * y7 +
            x3 * (-y1 - y2 + y4 + y7) + x2 * (-y1 + y3 - y5 + y7));
        dx7 = sz2 * (y4 * z0 + y5 * z0 + y0 * z2 - y5 * z2 - y6 * z2 + y0 * z3 + y4 * z3 - y6 * z3 -
            y0 * z4 + y5 * z4 + y6 * z4 - y0 * z5 - y4 * z5 + y6 * z5 - y4 * z6 - y5 * z6 +
            y3 * (-z0 + z2 - z4 + z6) + y2 * (-z0 - z3 + z5 + z6));
        dy7 = -sz2 * (x4 * z0 + x5 * z0 + x0 * z2 - x5 * z2 - x6 * z2 + x0 * z3 + x4 * z3 - x6 * z3 -
            x0 * z4 + x5 * z4 + x6 * z4 - x0 * z5 - x4 * z5 + x6 * z5 - x4 * z6 - x5 * z6 +
            x3 * (-z0 + z2 - z4 + z6) + x2 * (-z0 - z3 + z5 + z6));
        dz7 = sz2 * (x4 * y0 + x5 * y0 + x0 * y2 - x5 * y2 - x6 * y2 + x0 * y3 + x4 * y3 - x6 * y3 -
            x0 * y4 + x5 * y4 + x6 * y4 - x0 * y5 - x4 * y5 + x6 * y5 - x4 * y6 - x5 * y6 +
            x3 * (-y0 + y2 - y4 + y6) + x2 * (-y0 - y3 + y5 + y6));
        return;
    }
    // use Scaled Jacobian
    else if (vol > len3 * sJThres) {
        eFSJ = false;
        sz *= 4.0;
        dx0 = -sz * ((2 * vec[0][1] * (-z2 - z3 + z4 + z5) + 2 * (y2 + y3 - y4 - y5) * vec[0][2] + vec[2][1] * vec[1][2] - vec[1][1] * vec[2][2]) /
            sqrt((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) *
                (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) *
                (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])) -
            (0.5 * ((-vec[2][0] * vec[1][1] + vec[1][0] * vec[2][1]) *
                vec[0][2] + vec[0][1] *
                (-vec[1][0] * vec[2][2] + vec[2][0] * vec[1][2]) +
                vec[0][0] * (-vec[2][1] * vec[1][2] +
                    vec[1][1] * vec[2][2])) *
                (-2 * vec[2][0] * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] +
                    vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] +
                        vec[1][2] * vec[1][2]) - 2 * vec[1][0] *
                    (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) *
                    (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) -
                    2 * vec[0][0] * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] +
                        vec[1][2] * vec[1][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] +
                            vec[2][2] * vec[2][2]))) /
            pow((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) *
                (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) *
                (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]), 1.5));
        dy0 = -sz * ((2 * vec[0][0] * (z2 + z3 - z4 - z5) + vec[1][0] * vec[2][2] + 2 * (-x2 - x3 + x4 + x5) * vec[0][2] - vec[2][0] * vec[1][2]) /
            sqrt((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) *
                (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) *
                (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])) -
            (0.5 * ((-vec[2][0] * vec[1][1] + vec[1][0] * vec[2][1]) *
                vec[0][2] + vec[0][1] *
                (-vec[1][0] * vec[2][2] + vec[2][0] * vec[1][2]) +
                vec[0][0] * (-vec[2][1] * vec[1][2] +
                    vec[1][1] * vec[2][2])) *
                (-2 * vec[2][1] * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] +
                    vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] +
                        vec[1][2] * vec[1][2]) - 2 * vec[1][1] *
                    (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) *
                    (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) -
                    2 * vec[0][1] * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] +
                        vec[1][2] * vec[1][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] +
                            vec[2][2] * vec[2][2]))) /
            pow((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) *
                (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) *
                (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]), 1.5));
        dz0 = -sz * ((2 * vec[0][0] * (-y2 - y3 + y4 + y5) + 2 * (x2 + x3 - x4 - x5) * vec[0][1] + vec[2][0] * vec[1][1] - vec[1][0] * vec[2][1]) /
            sqrt((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) *
                (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) *
                (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])) -
            (0.5 * ((-vec[2][0] * vec[1][1] + vec[1][0] * vec[2][1]) *
                vec[0][2] + vec[0][1] *
                (-vec[1][0] * vec[2][2] + vec[2][0] * vec[1][2]) +
                vec[0][0] * (-vec[2][1] * vec[1][2] +
                    vec[1][1] * vec[2][2])) *
                (-2 * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) *
                    vec[2][2] * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] +
                        vec[1][2] * vec[1][2]) - 2 * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] +
                            vec[0][2] * vec[0][2]) * vec[1][2] *
                    (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) -
                    2 * vec[0][2] * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] +
                        vec[1][2] * vec[1][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] +
                            vec[2][2] * vec[2][2]))) /
            pow((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) *
                (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) *
                (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]), 1.5));
        dx1 = -sz * ((2 * vec[0][1] * (-z2 - z3 + z4 + z5) + 2 * (y2 + y3 - y4 - y5) * vec[0][2] - vec[2][1] * vec[1][2] + vec[1][1] * vec[2][2]) /
            sqrt((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) *
                (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) *
                (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])) -
            (0.5 * ((-vec[2][0] * vec[1][1] + vec[1][0] * vec[2][1]) *
                vec[0][2] + vec[0][1] *
                (-vec[1][0] * vec[2][2] + vec[2][0] * vec[1][2]) +
                vec[0][0] * (-vec[2][1] * vec[1][2] +
                    vec[1][1] * vec[2][2])) *
                (-2 * vec[2][0] * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] +
                    vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] +
                        vec[1][2] * vec[1][2]) - 2 * vec[1][0] *
                    (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) *
                    (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) +
                    2 * vec[0][0] * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] +
                        vec[1][2] * vec[1][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] +
                            vec[2][2] * vec[2][2]))) /
            pow((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) *
                (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) *
                (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]), 1.5));
        dy1 = -sz * ((2 * vec[0][0] * (z2 + z3 - z4 - z5) - vec[1][0] * vec[2][2] + 2 * (-x2 - x3 + x4 + x5) * vec[0][2] + vec[2][0] * vec[1][2]) /
            sqrt((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) *
                (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) *
                (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])) -
            (0.5 * ((-vec[2][0] * vec[1][1] + vec[1][0] * vec[2][1]) *
                vec[0][2] + vec[0][1] *
                (-vec[1][0] * vec[2][2] + vec[2][0] * vec[1][2]) +
                vec[0][0] * (-vec[2][1] * vec[1][2] +
                    vec[1][1] * vec[2][2])) *
                (-2 * vec[2][1] * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] +
                    vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] +
                        vec[1][2] * vec[1][2]) - 2 * vec[1][1] *
                    (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) *
                    (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) +
                    2 * vec[0][1] * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] +
                        vec[1][2] * vec[1][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] +
                            vec[2][2] * vec[2][2]))) /
            pow((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) *
                (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) *
                (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]), 1.5));
        dz1 = -sz * ((2 * vec[0][0] * (-y2 - y3 + y4 + y5) + 2 * (x2 + x3 - x4 - x5) * vec[0][1] - vec[2][0] * vec[1][1] + vec[1][0] * vec[2][1]) /
            sqrt((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) *
                (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) *
                (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])) -
            (0.5 * ((-vec[2][0] * vec[1][1] + vec[1][0] * vec[2][1]) *
                vec[0][2] + vec[0][1] *
                (-vec[1][0] * vec[2][2] + vec[2][0] * vec[1][2]) +
                vec[0][0] * (-vec[2][1] * vec[1][2] +
                    vec[1][1] * vec[2][2])) *
                (-2 * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) *
                    vec[2][2] * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] +
                        vec[1][2] * vec[1][2]) - 2 * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] +
                            vec[0][2] * vec[0][2]) * vec[1][2] *
                    (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) +
                    2 * vec[0][2] * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] +
                        vec[1][2] * vec[1][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] +
                            vec[2][2] * vec[2][2]))) /
            pow((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) *
                (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) *
                (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]), 1.5));
        dx2 = -sz * ((2 * vec[0][1] * (z0 + z1 - z6 - z7) + 2 * (-y0 - y1 + y6 + y7) * vec[0][2] - vec[2][1] * vec[1][2] + vec[1][1] * vec[2][2]) /
            sqrt((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) *
                (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) *
                (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])) -
            (0.5 * ((-vec[2][0] * vec[1][1] + vec[1][0] * vec[2][1]) *
                vec[0][2] + vec[0][1] *
                (-vec[1][0] * vec[2][2] + vec[2][0] * vec[1][2]) +
                vec[0][0] * (-vec[2][1] * vec[1][2] +
                    vec[1][1] * vec[2][2])) *
                (-2 * vec[2][0] * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] +
                    vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] +
                        vec[1][2] * vec[1][2]) + 2 * vec[1][0] *
                    (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) *
                    (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) +
                    2 * vec[0][0] * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] +
                        vec[1][2] * vec[1][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] +
                            vec[2][2] * vec[2][2]))) /
            pow((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) *
                (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) *
                (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]), 1.5));
        dy2 = -sz * ((-vec[1][0] * vec[2][2] + 2 * (x0 + x1 - x6 - x7) * vec[0][2] +
            vec[2][0] * vec[1][2] + 2 * vec[0][0] * (-z0 - z1 + z6 + z7)) /
            sqrt((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) *
                (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) *
                (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])) -
            (0.5 * ((-vec[2][0] * vec[1][1] + vec[1][0] * vec[2][1]) *
                vec[0][2] + vec[0][1] *
                (-vec[1][0] * vec[2][2] + vec[2][0] * vec[1][2]) +
                vec[0][0] * (-vec[2][1] * vec[1][2] +
                    vec[1][1] * vec[2][2])) *
                (-2 * vec[2][1] * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] +
                    vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] +
                        vec[1][2] * vec[1][2]) + 2 * vec[1][1] *
                    (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) *
                    (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) +
                    2 * vec[0][1] * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] +
                        vec[1][2] * vec[1][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] +
                            vec[2][2] * vec[2][2]))) /
            pow((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) *
                (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) *
                (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]), 1.5));
        dz2 = -sz * ((2 * vec[0][0] * (y0 + y1 - y6 - y7) + 2 * (-x0 - x1 + x6 + x7) * vec[0][1] - vec[2][0] * vec[1][1] + vec[1][0] * vec[2][1]) /
            sqrt((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) *
                (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) *
                (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])) -
            (0.5 * ((-vec[2][0] * vec[1][1] + vec[1][0] * vec[2][1]) *
                vec[0][2] + vec[0][1] *
                (-vec[1][0] * vec[2][2] + vec[2][0] * vec[1][2]) +
                vec[0][0] * (-vec[2][1] * vec[1][2] +
                    vec[1][1] * vec[2][2])) *
                (-2 * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) *
                    vec[2][2] * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] +
                        vec[1][2] * vec[1][2]) + 2 * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] +
                            vec[0][2] * vec[0][2]) * vec[1][2] *
                    (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) +
                    2 * vec[0][2] * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] +
                        vec[1][2] * vec[1][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] +
                            vec[2][2] * vec[2][2]))) /
            pow((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) *
                (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) *
                (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]), 1.5));
        dx3 = -sz * ((2 * vec[0][1] * (z0 + z1 - z6 - z7) + 2 * (-y0 - y1 + y6 + y7) * vec[0][2] + vec[2][1] * vec[1][2] - vec[1][1] * vec[2][2]) /
            sqrt((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) *
                (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) *
                (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])) -
            (0.5 * ((-vec[2][0] * vec[1][1] + vec[1][0] * vec[2][1]) *
                vec[0][2] + vec[0][1] *
                (-vec[1][0] * vec[2][2] + vec[2][0] * vec[1][2]) +
                vec[0][0] * (-vec[2][1] * vec[1][2] +
                    vec[1][1] * vec[2][2])) *
                (-2 * vec[2][0] * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] +
                    vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] +
                        vec[1][2] * vec[1][2]) + 2 * vec[1][0] *
                    (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) *
                    (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) -
                    2 * vec[0][0] * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] +
                        vec[1][2] * vec[1][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] +
                            vec[2][2] * vec[2][2]))) /
            pow((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) *
                (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) *
                (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]), 1.5));
        dy3 = -sz * ((vec[1][0] * vec[2][2] + 2 * (x0 + x1 - x6 - x7) * vec[0][2] -
            vec[2][0] * vec[1][2] + 2 * vec[0][0] * (-z0 - z1 + z6 + z7)) /
            sqrt((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) *
                (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) *
                (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])) -
            (0.5 * ((-vec[2][0] * vec[1][1] + vec[1][0] * vec[2][1]) *
                vec[0][2] + vec[0][1] *
                (-vec[1][0] * vec[2][2] + vec[2][0] * vec[1][2]) +
                vec[0][0] * (-vec[2][1] * vec[1][2] +
                    vec[1][1] * vec[2][2])) *
                (-2 * vec[2][1] * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] +
                    vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] +
                        vec[1][2] * vec[1][2]) + 2 * vec[1][1] *
                    (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) *
                    (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) -
                    2 * vec[0][1] * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] +
                        vec[1][2] * vec[1][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] +
                            vec[2][2] * vec[2][2]))) /
            pow((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) *
                (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) *
                (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]), 1.5));
        dz3 = -sz * ((2 * vec[0][0] * (y0 + y1 - y6 - y7) + 2 * (-x0 - x1 + x6 + x7) * vec[0][1] + vec[2][0] * vec[1][1] - vec[1][0] * vec[2][1]) /
            sqrt((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) *
                (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) *
                (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])) -
            (0.5 * ((-vec[2][0] * vec[1][1] + vec[1][0] * vec[2][1]) *
                vec[0][2] + vec[0][1] *
                (-vec[1][0] * vec[2][2] + vec[2][0] * vec[1][2]) +
                vec[0][0] * (-vec[2][1] * vec[1][2] +
                    vec[1][1] * vec[2][2])) *
                (-2 * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) *
                    vec[2][2] * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] +
                        vec[1][2] * vec[1][2]) + 2 * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] +
                            vec[0][2] * vec[0][2]) * vec[1][2] *
                    (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) -
                    2 * vec[0][2] * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] +
                        vec[1][2] * vec[1][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] +
                            vec[2][2] * vec[2][2]))) /
            pow((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) *
                (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) *
                (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]), 1.5));
        dx4 = -sz * ((2 * (y0 + y1 - y6 - y7) * vec[0][2] + vec[2][1] * vec[1][2] -
            vec[1][1] * vec[2][2] + 2 * vec[0][1] * (-z0 - z1 + z6 + z7)) /
            sqrt((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) *
                (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) *
                (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])) -
            (0.5 * ((-vec[2][0] * vec[1][1] + vec[1][0] * vec[2][1]) *
                vec[0][2] + vec[0][1] *
                (-vec[1][0] * vec[2][2] + vec[2][0] * vec[1][2]) +
                vec[0][0] * (-vec[2][1] * vec[1][2] +
                    vec[1][1] * vec[2][2])) *
                (2 * vec[2][0] * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] +
                    vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] +
                        vec[1][2] * vec[1][2]) - 2 * vec[1][0] *
                    (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) *
                    (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) -
                    2 * vec[0][0] * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] +
                        vec[1][2] * vec[1][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] +
                            vec[2][2] * vec[2][2]))) /
            pow((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) *
                (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) *
                (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]), 1.5));
        dy4 = -sz * ((2 * vec[0][0] * (z0 + z1 - z6 - z7) + vec[1][0] * vec[2][2] + 2 * (-x0 - x1 + x6 + x7) * vec[0][2] - vec[2][0] * vec[1][2]) /
            sqrt((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) *
                (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) *
                (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])) -
            (0.5 * ((-vec[2][0] * vec[1][1] + vec[1][0] * vec[2][1]) *
                vec[0][2] + vec[0][1] *
                (-vec[1][0] * vec[2][2] + vec[2][0] * vec[1][2]) +
                vec[0][0] * (-vec[2][1] * vec[1][2] +
                    vec[1][1] * vec[2][2])) *
                (2 * vec[2][1] * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] +
                    vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] +
                        vec[1][2] * vec[1][2]) - 2 * vec[1][1] *
                    (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) *
                    (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) -
                    2 * vec[0][1] * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] +
                        vec[1][2] * vec[1][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] +
                            vec[2][2] * vec[2][2]))) /
            pow((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) *
                (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) *
                (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]), 1.5));
        dz4 = -sz * ((2 * (x0 + x1 - x6 - x7) * vec[0][1] + vec[2][0] * vec[1][1] -
            vec[1][0] * vec[2][1] + 2 * vec[0][0] * (-y0 - y1 + y6 + y7)) /
            sqrt((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) *
                (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) *
                (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])) -
            (0.5 * ((-vec[2][0] * vec[1][1] + vec[1][0] * vec[2][1]) *
                vec[0][2] + vec[0][1] *
                (-vec[1][0] * vec[2][2] + vec[2][0] * vec[1][2]) +
                vec[0][0] * (-vec[2][1] * vec[1][2] +
                    vec[1][1] * vec[2][2])) *
                (2 * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) *
                    vec[2][2] * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] +
                        vec[1][2] * vec[1][2]) - 2 * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] +
                            vec[0][2] * vec[0][2]) * vec[1][2] *
                    (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) -
                    2 * vec[0][2] * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] +
                        vec[1][2] * vec[1][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] +
                            vec[2][2] * vec[2][2]))) /
            pow((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) *
                (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) *
                (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]), 1.5));
        dx5 = -sz * ((2 * (y0 + y1 - y6 - y7) * vec[0][2] - vec[2][1] * vec[1][2] +
            vec[1][1] * vec[2][2] + 2 * vec[0][1] * (-z0 - z1 + z6 + z7)) /
            sqrt((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) *
                (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) *
                (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])) -
            (0.5 * ((-vec[2][0] * vec[1][1] + vec[1][0] * vec[2][1]) *
                vec[0][2] + vec[0][1] *
                (-vec[1][0] * vec[2][2] + vec[2][0] * vec[1][2]) +
                vec[0][0] * (-vec[2][1] * vec[1][2] +
                    vec[1][1] * vec[2][2])) *
                (2 * vec[2][0] * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] +
                    vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] +
                        vec[1][2] * vec[1][2]) - 2 * vec[1][0] *
                    (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) *
                    (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) +
                    2 * vec[0][0] * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] +
                        vec[1][2] * vec[1][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] +
                            vec[2][2] * vec[2][2]))) /
            pow((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) *
                (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) *
                (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]), 1.5));
        dy5 = -sz * ((2 * vec[0][0] * (z0 + z1 - z6 - z7) - vec[1][0] * vec[2][2] + 2 * (-x0 - x1 + x6 + x7) * vec[0][2] + vec[2][0] * vec[1][2]) /
            sqrt((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) *
                (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) *
                (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])) -
            (0.5 * ((-vec[2][0] * vec[1][1] + vec[1][0] * vec[2][1]) *
                vec[0][2] + vec[0][1] *
                (-vec[1][0] * vec[2][2] + vec[2][0] * vec[1][2]) +
                vec[0][0] * (-vec[2][1] * vec[1][2] +
                    vec[1][1] * vec[2][2])) *
                (2 * vec[2][1] * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] +
                    vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] +
                        vec[1][2] * vec[1][2]) - 2 * vec[1][1] *
                    (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) *
                    (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) +
                    2 * vec[0][1] * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] +
                        vec[1][2] * vec[1][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] +
                            vec[2][2] * vec[2][2]))) /
            pow((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) *
                (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) *
                (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]), 1.5));
        dz5 = -sz * ((2 * (x0 + x1 - x6 - x7) * vec[0][1] - vec[2][0] * vec[1][1] +
            vec[1][0] * vec[2][1] + 2 * vec[0][0] * (-y0 - y1 + y6 + y7)) /
            sqrt((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) *
                (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) *
                (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])) -
            (0.5 * ((-vec[2][0] * vec[1][1] + vec[1][0] * vec[2][1]) *
                vec[0][2] + vec[0][1] *
                (-vec[1][0] * vec[2][2] + vec[2][0] * vec[1][2]) +
                vec[0][0] * (-vec[2][1] * vec[1][2] +
                    vec[1][1] * vec[2][2])) *
                (2 * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) *
                    vec[2][2] * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] +
                        vec[1][2] * vec[1][2]) - 2 * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] +
                            vec[0][2] * vec[0][2]) * vec[1][2] *
                    (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) +
                    2 * vec[0][2] * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] +
                        vec[1][2] * vec[1][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] +
                            vec[2][2] * vec[2][2]))) /
            pow((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) *
                (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) *
                (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]), 1.5));
        dx6 = -sz * ((2 * vec[0][1] * (z2 + z3 - z4 - z5) + 2 * (-y2 - y3 + y4 + y5) * vec[0][2] - vec[2][1] * vec[1][2] + vec[1][1] * vec[2][2]) /
            sqrt((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) *
                (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) *
                (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])) -
            (0.5 * ((-vec[2][0] * vec[1][1] + vec[1][0] * vec[2][1]) *
                vec[0][2] + vec[0][1] *
                (-vec[1][0] * vec[2][2] + vec[2][0] * vec[1][2]) +
                vec[0][0] * (-vec[2][1] * vec[1][2] +
                    vec[1][1] * vec[2][2])) *
                (2 * vec[2][0] * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] +
                    vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] +
                        vec[1][2] * vec[1][2]) + 2 * vec[1][0] *
                    (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) *
                    (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) +
                    2 * vec[0][0] * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] +
                        vec[1][2] * vec[1][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] +
                            vec[2][2] * vec[2][2]))) /
            pow((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) *
                (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) *
                (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]), 1.5));
        dy6 = -sz * ((2 * vec[0][0] * (-z2 - z3 + z4 + z5) - vec[1][0] * vec[2][2] + 2 * (x2 + x3 - x4 - x5) * vec[0][2] + vec[2][0] * vec[1][2]) /
            sqrt((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) *
                (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) *
                (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])) -
            (0.5 * ((-vec[2][0] * vec[1][1] + vec[1][0] * vec[2][1]) *
                vec[0][2] + vec[0][1] *
                (-vec[1][0] * vec[2][2] + vec[2][0] * vec[1][2]) +
                vec[0][0] * (-vec[2][1] * vec[1][2] +
                    vec[1][1] * vec[2][2])) *
                (2 * vec[2][1] * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] +
                    vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] +
                        vec[1][2] * vec[1][2]) + 2 * vec[1][1] *
                    (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) *
                    (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) +
                    2 * vec[0][1] * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] +
                        vec[1][2] * vec[1][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] +
                            vec[2][2] * vec[2][2]))) /
            pow((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) *
                (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) *
                (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]), 1.5));
        dz6 = -sz * ((2 * vec[0][0] * (y2 + y3 - y4 - y5) + 2 * (-x2 - x3 + x4 + x5) * vec[0][1] - vec[2][0] * vec[1][1] + vec[1][0] * vec[2][1]) /
            sqrt((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) *
                (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) *
                (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])) -
            (0.5 * ((-vec[2][0] * vec[1][1] + vec[1][0] * vec[2][1]) *
                vec[0][2] + vec[0][1] *
                (-vec[1][0] * vec[2][2] + vec[2][0] * vec[1][2]) +
                vec[0][0] * (-vec[2][1] * vec[1][2] +
                    vec[1][1] * vec[2][2])) *
                (2 * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) *
                    vec[2][2] * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] +
                        vec[1][2] * vec[1][2]) + 2 * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] +
                            vec[0][2] * vec[0][2]) * vec[1][2] *
                    (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) +
                    2 * vec[0][2] * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] +
                        vec[1][2] * vec[1][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] +
                            vec[2][2] * vec[2][2]))) /
            pow((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) *
                (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) *
                (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]), 1.5));
        dx7 = -sz * ((2 * vec[0][1] * (z2 + z3 - z4 - z5) + 2 * (-y2 - y3 + y4 + y5) * vec[0][2] + vec[2][1] * vec[1][2] - vec[1][1] * vec[2][2]) /
            sqrt((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) *
                (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) *
                (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])) -
            (0.5 * ((-vec[2][0] * vec[1][1] + vec[1][0] * vec[2][1]) *
                vec[0][2] + vec[0][1] *
                (-vec[1][0] * vec[2][2] + vec[2][0] * vec[1][2]) +
                vec[0][0] * (-vec[2][1] * vec[1][2] +
                    vec[1][1] * vec[2][2])) *
                (2 * vec[2][0] * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] +
                    vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] +
                        vec[1][2] * vec[1][2]) + 2 * vec[1][0] *
                    (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) *
                    (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) -
                    2 * vec[0][0] * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] +
                        vec[1][2] * vec[1][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] +
                            vec[2][2] * vec[2][2]))) /
            pow((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) *
                (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) *
                (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]), 1.5));
        dy7 = -sz * ((2 * vec[0][0] * (-z2 - z3 + z4 + z5) + vec[1][0] * vec[2][2] + 2 * (x2 + x3 - x4 - x5) * vec[0][2] - vec[2][0] * vec[1][2]) /
            sqrt((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) *
                (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) *
                (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])) -
            (0.5 * ((-vec[2][0] * vec[1][1] + vec[1][0] * vec[2][1]) *
                vec[0][2] + vec[0][1] *
                (-vec[1][0] * vec[2][2] + vec[2][0] * vec[1][2]) +
                vec[0][0] * (-vec[2][1] * vec[1][2] +
                    vec[1][1] * vec[2][2])) *
                (2 * vec[2][1] * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] +
                    vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] +
                        vec[1][2] * vec[1][2]) + 2 * vec[1][1] *
                    (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) *
                    (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) -
                    2 * vec[0][1] * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] +
                        vec[1][2] * vec[1][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] +
                            vec[2][2] * vec[2][2]))) /
            pow((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) *
                (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) *
                (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]), 1.5));
        dz7 = -sz * ((2 * vec[0][0] * (y2 + y3 - y4 - y5) + 2 * (-x2 - x3 + x4 + x5) * vec[0][1] + vec[2][0] * vec[1][1] - vec[1][0] * vec[2][1]) /
            sqrt((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) *
                (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) *
                (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])) -
            (0.5 * ((-vec[2][0] * vec[1][1] + vec[1][0] * vec[2][1]) *
                vec[0][2] + vec[0][1] *
                (-vec[1][0] * vec[2][2] + vec[2][0] * vec[1][2]) +
                vec[0][0] * (-vec[2][1] * vec[1][2] +
                    vec[1][1] * vec[2][2])) *
                (2 * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) *
                    vec[2][2] * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] +
                        vec[1][2] * vec[1][2]) + 2 * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] +
                            vec[0][2] * vec[0][2]) * vec[1][2] *
                    (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) -
                    2 * vec[0][2] * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] +
                        vec[1][2] * vec[1][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] +
                            vec[2][2] * vec[2][2]))) /
            pow((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) *
                (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) *
                (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]), 1.5));
        return;
    }

    // 0 1 3 4
    vec[0][0] = x1 - x0; vec[0][1] = y1 - y0; vec[0][2] = z1 - z0;
    vec[1][0] = x3 - x0; vec[1][1] = y3 - y0; vec[1][2] = z3 - z0;
    vec[2][0] = x4 - x0; vec[2][1] = y4 - y0; vec[2][2] = z4 - z0;
    vol = vec[0][0] * (vec[2][1] * vec[1][2] - vec[1][1] * vec[2][2]) +
        vec[0][1] * (vec[2][2] * vec[1][0] - vec[1][2] * vec[2][0]) +
        vec[0][2] * (vec[2][0] * vec[1][1] - vec[1][0] * vec[2][1]);
    len3 = sqrt(DOT(vec[0], vec[0]) * DOT(vec[1], vec[1]) * DOT(vec[2], vec[2]));
    if (vol >= 0) {
        eFSJ = false;
        dx0 = -sz2 * (y4 * (-z1 + z3) + y3 * (z1 - z4) + y1 * (-z3 + z4));
        dy0 = -sz2 * (x4 * (z1 - z3) + x1 * (z3 - z4) + x3 * (-z1 + z4));
        dz0 = -sz2 * (x4 * (-y1 + y3) + x3 * (y1 - y4) + x1 * (-y3 + y4));
        dx1 = -sz2 * (-y4 * vec[1][2] + y0 * (z3 - z4) + y3 * vec[2][2]);
        dy1 = -sz2 * (x4 * vec[1][2] - x3 * vec[2][2] + x0 * (-z3 + z4));
        dz1 = -sz2 * (-x4 * vec[1][1] + x0 * (y3 - y4) + x3 * vec[2][1]);
        dx3 = -sz2 * (y4 * vec[0][2] - y1 * vec[2][2] + y0 * (-z1 + z4));
        dy3 = -sz2 * (-x4 * vec[0][2] + x0 * (z1 - z4) + x1 * vec[2][2]);
        dz3 = -sz2 * (x4 * vec[0][1] - x1 * vec[2][1] + x0 * (-y1 + y4));
        dx4 = -sz2 * (-y3 * vec[0][2] + y0 * (z1 - z3) + y1 * vec[1][2]);
        dy4 = -sz2 * (x3 * vec[0][2] - x1 * vec[1][2] + x0 * (-z1 + z3));
        dz4 = -sz2 * (-x3 * vec[0][1] + x0 * (y1 - y3) + x1 * vec[1][1]);
        return;
    }
    else if (vol > len3 * sJThres) {
        eFSJ = false;
        dx0 = -sz * ((y4 * (-z1 + z3) + y3 * (z1 - z4) + y1 * (-z3 + z4)) /
            sqrt((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]))
            - ((-vec[2][0] * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) -
                vec[1][0] * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) -
                vec[0][0] * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])) *
                (-(x4 * vec[1][1] - x3 * vec[2][1] + x0 * (-y3 + y4)) * vec[0][2] - vec[0][1] * (-x4 * vec[1][2] + x0 * (z3 - z4) + x3 * vec[2][2]) - vec[0][0] * (y4 * vec[1][2] - y3 * vec[2][2] + y0 * (-z3 + z4)))) /
            pow((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]),
                1.5));
        dy0 = -sz * ((x4 * (z1 - z3) + x1 * (z3 - z4) + x3 * (-z1 + z4)) /
            sqrt((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]))
            - ((-vec[2][1] * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) -
                vec[1][1] * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) -
                vec[0][1] * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])) *
                (-(x4 * vec[1][1] - x3 * vec[2][1] + x0 * (-y3 + y4)) * vec[0][2] - vec[0][1] * (-x4 * vec[1][2] + x0 * (z3 - z4) + x3 * vec[2][2]) - vec[0][0] * (y4 * vec[1][2] - y3 * vec[2][2] + y0 * (-z3 + z4)))) /
            pow((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]),
                1.5));
        dz0 = -sz * ((x4 * (-y1 + y3) + x3 * (y1 - y4) + x1 * (-y3 + y4)) /
            sqrt((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]))
            - ((-vec[0][2] * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) -
                (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * vec[1][2] * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) -
                (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * vec[2][2]) *
                (-(x4 * vec[1][1] - x3 * vec[2][1] + x0 * (-y3 + y4)) * vec[0][2] - vec[0][1] * (-x4 * vec[1][2] + x0 * (z3 - z4) + x3 * vec[2][2]) - vec[0][0] * (y4 * vec[1][2] - y3 * vec[2][2] + y0 * (-z3 + z4)))) /
            pow((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]),
                1.5));
        dx1 = -sz * ((-y4 * vec[1][2] + y0 * (z3 - z4) + y3 * vec[2][2] - (vec[0][0] * (-(x4 * vec[1][1] - x3 * vec[2][1] + x0 * (-y3 + y4)) * vec[0][2] - vec[0][1] * (-x4 * vec[1][2] + x0 * (z3 - z4) + x3 * vec[2][2])
            - vec[0][0] * (y4 * vec[1][2] - y3 * vec[2][2] + y0 * (-z3 + z4)))) / (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2])) /
            sqrt((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])));
        dy1 = -sz * ((x4 * vec[1][2] - x3 * vec[2][2] + x0 * (-z3 + z4) - (vec[0][1] * (-(x4 * vec[1][1] - x3 * vec[2][1] + x0 * (-y3 + y4)) * vec[0][2] - vec[0][1] * (-x4 * vec[1][2] + x0 * (z3 - z4) + x3 * vec[2][2])
            - vec[0][0] * (y4 * vec[1][2] - y3 * vec[2][2] + y0 * (-z3 + z4)))) / (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2])) /
            sqrt((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])));
        dz1 = -sz * ((-x4 * vec[1][1] + x0 * (y3 - y4) + x3 * vec[2][1] - (vec[0][2] * (-(x4 * vec[1][1] - x3 * vec[2][1] + x0 * (-y3 + y4)) * vec[0][2] - vec[0][1] * (-x4 * vec[1][2] + x0 * (z3 - z4) + x3 * vec[2][2])
            - vec[0][0] * (y4 * vec[1][2] - y3 * vec[2][2] + y0 * (-z3 + z4)))) / (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2])) /
            sqrt((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])));
        dx3 = -sz * ((y4 * vec[0][2] - y1 * vec[2][2] + y0 * (-z1 + z4) - (vec[1][0] * (-(x4 * vec[1][1] - x3 * vec[2][1] + x0 * (-y3 + y4)) * vec[0][2] - vec[0][1] * (-x4 * vec[1][2] + x0 * (z3 - z4) + x3 * vec[2][2])
            - vec[0][0] * (y4 * vec[1][2] - y3 * vec[2][2] + y0 * (-z3 + z4)))) / (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2])) /
            sqrt((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])));
        dy3 = -sz * ((-x4 * vec[0][2] + x0 * (z1 - z4) + x1 * vec[2][2] - (vec[1][1] * (-(x4 * vec[1][1] - x3 * vec[2][1] + x0 * (-y3 + y4)) * vec[0][2] - vec[0][1] * (-x4 * vec[1][2] + x0 * (z3 - z4) + x3 * vec[2][2])
            - vec[0][0] * (y4 * vec[1][2] - y3 * vec[2][2] + y0 * (-z3 + z4)))) / (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2])) /
            sqrt((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])));
        dz3 = -sz * ((x4 * vec[0][1] - x1 * vec[2][1] + x0 * (-y1 + y4) - (vec[1][2] * (-(x4 * vec[1][1] - x3 * vec[2][1] + x0 * (-y3 + y4)) * vec[0][2] - vec[0][1] * (-x4 * vec[1][2] + x0 * (z3 - z4) + x3 * vec[2][2])
            - vec[0][0] * (y4 * vec[1][2] - y3 * vec[2][2] + y0 * (-z3 + z4)))) / (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2])) /
            sqrt((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])));
        dx4 = -sz * ((-y3 * vec[0][2] + y0 * (z1 - z3) + y1 * vec[1][2] - (vec[2][0] * (-(x4 * vec[1][1] - x3 * vec[2][1] + x0 * (-y3 + y4)) * vec[0][2] - vec[0][1] * (-x4 * vec[1][2] + x0 * (z3 - z4) + x3 * vec[2][2])
            - vec[0][0] * (y4 * vec[1][2] - y3 * vec[2][2] + y0 * (-z3 + z4)))) / (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])) /
            sqrt((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])));
        dy4 = -sz * ((x3 * vec[0][2] - x1 * vec[1][2] + x0 * (-z1 + z3) - (vec[2][1] * (-(x4 * vec[1][1] - x3 * vec[2][1] + x0 * (-y3 + y4)) * vec[0][2] - vec[0][1] * (-x4 * vec[1][2] + x0 * (z3 - z4) + x3 * vec[2][2])
            - vec[0][0] * (y4 * vec[1][2] - y3 * vec[2][2] + y0 * (-z3 + z4)))) / (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])) /
            sqrt((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])));
        dz4 = -sz * ((-x3 * vec[0][1] + x0 * (y1 - y3) + x1 * vec[1][1] - (vec[2][2] * (-(x4 * vec[1][1] - x3 * vec[2][1] + x0 * (-y3 + y4)) * vec[0][2] - vec[0][1] * (-x4 * vec[1][2] + x0 * (z3 - z4) + x3 * vec[2][2])
            - vec[0][0] * (y4 * vec[1][2] - y3 * vec[2][2] + y0 * (-z3 + z4)))) / (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])) /
            sqrt((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])));
        return;
    }

    // 1 2 0 5
    vec[0][0] = x2 - x1; vec[0][1] = y2 - y1; vec[0][2] = z2 - z1;
    vec[1][0] = x0 - x1; vec[1][1] = y0 - y1; vec[1][2] = z0 - z1;
    vec[2][0] = x5 - x1; vec[2][1] = y5 - y1; vec[2][2] = z5 - z1;
    vol = vec[0][0] * (vec[2][1] * vec[1][2] - vec[1][1] * vec[2][2]) +
        vec[0][1] * (vec[2][2] * vec[1][0] - vec[1][2] * vec[2][0]) +
        vec[0][2] * (vec[2][0] * vec[1][1] - vec[1][0] * vec[2][1]);
    len3 = sqrt(DOT(vec[0], vec[0]) * DOT(vec[1], vec[1]) * DOT(vec[2], vec[2]));
    if (vol >= 0) {
        eFSJ = false;
        dx0 = -sz2 * (y5 * vec[0][2] - y2 * vec[2][2] + y1 * (-z2 + z5));
        dy0 = -sz2 * (-x5 * vec[0][2] + x1 * (z2 - z5) + x2 * vec[2][2]);
        dz0 = -sz2 * (x5 * vec[0][1] - x2 * vec[2][1] + x1 * (-y2 + y5));
        dx1 = -sz2 * (y5 * (z0 - z2) + y0 * (z2 - z5) + y2 * (-z0 + z5));
        dy1 = -sz2 * (x5 * (-z0 + z2) + x2 * (z0 - z5) + x0 * (-z2 + z5));
        dz1 = -sz2 * (x5 * (y0 - y2) + x0 * (y2 - y5) + x2 * (-y0 + y5));
        dx2 = -sz2 * (-y5 * vec[1][2] + y1 * (z0 - z5) + y0 * vec[2][2]);
        dy2 = -sz2 * (x5 * vec[1][2] - x0 * vec[2][2] + x1 * (-z0 + z5));
        dz2 = -sz2 * (-x5 * vec[1][1] + x1 * (y0 - y5) + x0 * vec[2][1]);
        dx5 = -sz2 * (y2 * vec[1][2] - y0 * vec[0][2] + y1 * (-z0 + z2));
        dy5 = -sz2 * (-x2 * vec[1][2] + x1 * (z0 - z2) + x0 * vec[0][2]);
        dz5 = -sz2 * (x2 * vec[1][1] - x0 * vec[0][1] + x1 * (-y0 + y2));
        return;
    }
    else if (vol > len3 * sJThres) {
        eFSJ = false;
        dx0 = -sz * ((y5 * vec[0][2] - y2 * vec[2][2] + y1 * (-z2 + z5) - (vec[1][0] * (-(x5 * vec[1][1] - x0 * vec[2][1] + x1 * (-y0 + y5)) * vec[0][2] - vec[0][0] * (y5 * vec[1][2] - y0 * vec[2][2] + y1 * (-z0 + z5)) -
            vec[0][1] * (-x5 * vec[1][2] + x1 * (z0 - z5) + x0 * vec[2][2]))) / (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2])) /
            sqrt((vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])));
        dy0 = -sz * ((-x5 * vec[0][2] + x1 * (z2 - z5) + x2 * vec[2][2] - (vec[1][1] * (-(x5 * vec[1][1] - x0 * vec[2][1] + x1 * (-y0 + y5)) * vec[0][2] - vec[0][0] * (y5 * vec[1][2] - y0 * vec[2][2] + y1 * (-z0 + z5)) -
            vec[0][1] * (-x5 * vec[1][2] + x1 * (z0 - z5) + x0 * vec[2][2]))) / (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2])) /
            sqrt((vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])));
        dz0 = -sz * ((x5 * vec[0][1] - x2 * vec[2][1] + x1 * (-y2 + y5) - (vec[1][2] * (-(x5 * vec[1][1] - x0 * vec[2][1] + x1 * (-y0 + y5)) * vec[0][2] - vec[0][0] * (y5 * vec[1][2] - y0 * vec[2][2] + y1 * (-z0 + z5)) -
            vec[0][1] * (-x5 * vec[1][2] + x1 * (z0 - z5) + x0 * vec[2][2]))) / (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2])) /
            sqrt((vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])));
        dx1 = -sz * ((y5 * (z0 - z2) + y0 * (z2 - z5) + y2 * (-z0 + z5)) /
            sqrt((vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]))
            - ((-vec[2][0] * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) -
                vec[0][0] * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) -
                vec[1][0] * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])) *
                (-(x5 * vec[1][1] - x0 * vec[2][1] + x1 * (-y0 + y5)) * vec[0][2] - vec[0][0] * (y5 * vec[1][2] - y0 * vec[2][2] + y1 * (-z0 + z5)) - vec[0][1] * (-x5 * vec[1][2] + x1 * (z0 - z5) + x0 * vec[2][2]))) /
            pow((vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]),
                1.5));
        dy1 = -sz * ((x5 * (-z0 + z2) + x2 * (z0 - z5) + x0 * (-z2 + z5)) /
            sqrt((vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]))
            - ((-vec[2][1] * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) -
                vec[0][1] * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) -
                vec[1][1] * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])) *
                (-(x5 * vec[1][1] - x0 * vec[2][1] + x1 * (-y0 + y5)) * vec[0][2] - vec[0][0] * (y5 * vec[1][2] - y0 * vec[2][2] + y1 * (-z0 + z5)) - vec[0][1] * (-x5 * vec[1][2] + x1 * (z0 - z5) + x0 * vec[2][2]))) /
            pow((vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]),
                1.5));
        dz1 = -sz * ((x5 * (y0 - y2) + x0 * (y2 - y5) + x2 * (-y0 + y5)) /
            sqrt((vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]))
            - ((-vec[1][2] * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) -
                (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * vec[0][2] * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) -
                (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * vec[2][2]) *
                (-(x5 * vec[1][1] - x0 * vec[2][1] + x1 * (-y0 + y5)) * vec[0][2] - vec[0][0] * (y5 * vec[1][2] - y0 * vec[2][2] + y1 * (-z0 + z5)) - vec[0][1] * (-x5 * vec[1][2] + x1 * (z0 - z5) + x0 * vec[2][2]))) /
            pow((vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]),
                1.5));
        dx2 = -sz * ((-y5 * vec[1][2] + y1 * (z0 - z5) + y0 * vec[2][2] - (vec[0][0] * (-(x5 * vec[1][1] - x0 * vec[2][1] + x1 * (-y0 + y5)) * vec[0][2] - vec[0][0] * (y5 * vec[1][2] - y0 * vec[2][2] + y1 * (-z0 + z5)) -
            vec[0][1] * (-x5 * vec[1][2] + x1 * (z0 - z5) + x0 * vec[2][2]))) / (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2])) /
            sqrt((vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])));
        dy2 = -sz * ((x5 * vec[1][2] - x0 * vec[2][2] + x1 * (-z0 + z5) - (vec[0][1] * (-(x5 * vec[1][1] - x0 * vec[2][1] + x1 * (-y0 + y5)) * vec[0][2] - vec[0][0] * (y5 * vec[1][2] - y0 * vec[2][2] + y1 * (-z0 + z5)) -
            vec[0][1] * (-x5 * vec[1][2] + x1 * (z0 - z5) + x0 * vec[2][2]))) / (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2])) /
            sqrt((vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])));
        dz2 = -sz * ((-x5 * vec[1][1] + x1 * (y0 - y5) + x0 * vec[2][1] - (vec[0][2] * (-(x5 * vec[1][1] - x0 * vec[2][1] + x1 * (-y0 + y5)) * vec[0][2] - vec[0][0] * (y5 * vec[1][2] - y0 * vec[2][2] + y1 * (-z0 + z5)) -
            vec[0][1] * (-x5 * vec[1][2] + x1 * (z0 - z5) + x0 * vec[2][2]))) / (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2])) /
            sqrt((vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])));
        dx5 = -sz * ((y2 * vec[1][2] - y0 * vec[0][2] + y1 * (-z0 + z2) - (vec[2][0] * (-(x5 * vec[1][1] - x0 * vec[2][1] + x1 * (-y0 + y5)) * vec[0][2] - vec[0][0] * (y5 * vec[1][2] - y0 * vec[2][2] + y1 * (-z0 + z5)) -
            vec[0][1] * (-x5 * vec[1][2] + x1 * (z0 - z5) + x0 * vec[2][2]))) / (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])) /
            sqrt((vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])));
        dy5 = -sz * ((-x2 * vec[1][2] + x1 * (z0 - z2) + x0 * vec[0][2] - (vec[2][1] * (-(x5 * vec[1][1] - x0 * vec[2][1] + x1 * (-y0 + y5)) * vec[0][2] - vec[0][0] * (y5 * vec[1][2] - y0 * vec[2][2] + y1 * (-z0 + z5)) -
            vec[0][1] * (-x5 * vec[1][2] + x1 * (z0 - z5) + x0 * vec[2][2]))) / (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])) /
            sqrt((vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])));
        dz5 = -sz * ((x2 * vec[1][1] - x0 * vec[0][1] + x1 * (-y0 + y2) - (vec[2][2] * (-(x5 * vec[1][1] - x0 * vec[2][1] + x1 * (-y0 + y5)) * vec[0][2] - vec[0][0] * (y5 * vec[1][2] - y0 * vec[2][2] + y1 * (-z0 + z5)) -
            vec[0][1] * (-x5 * vec[1][2] + x1 * (z0 - z5) + x0 * vec[2][2]))) / (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])) /
            sqrt((vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])));
        return;
    }

    // 2 3 1 6
    vec[0][0] = x3 - x2; vec[0][1] = y3 - y2; vec[0][2] = z3 - z2;
    vec[1][0] = x1 - x2; vec[1][1] = y1 - y2; vec[1][2] = z1 - z2;
    vec[2][0] = x6 - x2; vec[2][1] = y6 - y2; vec[2][2] = z6 - z2;
    vol = vec[0][0] * (vec[2][1] * vec[1][2] - vec[1][1] * vec[2][2]) +
        vec[0][1] * (vec[2][2] * vec[1][0] - vec[1][2] * vec[2][0]) +
        vec[0][2] * (vec[2][0] * vec[1][1] - vec[1][0] * vec[2][1]);
    len3 = sqrt(DOT(vec[0], vec[0]) * DOT(vec[1], vec[1]) * DOT(vec[2], vec[2]));
    if (vol >= 0) {
        eFSJ = false;
        dx1 = -sz2 * (y6 * vec[0][2] - y3 * vec[2][2] + y2 * (-z3 + z6));
        dy1 = -sz2 * (-x6 * vec[0][2] + x2 * (z3 - z6) + x3 * vec[2][2]);
        dz1 = -sz2 * (x6 * vec[0][1] - x3 * vec[2][1] + x2 * (-y3 + y6));
        dx2 = -sz2 * (y6 * (z1 - z3) + y1 * (z3 - z6) + y3 * (-z1 + z6));
        dy2 = -sz2 * (x6 * (-z1 + z3) + x3 * (z1 - z6) + x1 * (-z3 + z6));
        dz2 = -sz2 * (x6 * (y1 - y3) + x1 * (y3 - y6) + x3 * (-y1 + y6));
        dx3 = -sz2 * (-y6 * vec[1][2] + y2 * (z1 - z6) + y1 * vec[2][2]);
        dy3 = -sz2 * (x6 * vec[1][2] - x1 * vec[2][2] + x2 * (-z1 + z6));
        dz3 = -sz2 * (-x6 * vec[1][1] + x2 * (y1 - y6) + x1 * vec[2][1]);
        dx6 = -sz2 * (y3 * vec[1][2] - y1 * vec[0][2] + y2 * (-z1 + z3));
        dy6 = -sz2 * (-x3 * vec[1][2] + x2 * (z1 - z3) + x1 * vec[0][2]);
        dz6 = -sz2 * (x3 * vec[1][1] - x1 * vec[0][1] + x2 * (-y1 + y3));
        return;
    }
    else if (vol > len3 * sJThres) {
        eFSJ = false;
        dx1 = -sz * ((y6 * vec[0][2] - y3 * vec[2][2] + y2 * (-z3 + z6) - (vec[1][0] * (-(x6 * vec[1][1] - x1 * vec[2][1] + x2 * (-y1 + y6)) * vec[0][2] - vec[0][0] * (y6 * vec[1][2] - y1 * vec[2][2] + y2 * (-z1 + z6)) -
            vec[0][1] * (-x6 * vec[1][2] + x2 * (z1 - z6) + x1 * vec[2][2]))) / (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2])) /
            sqrt((vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])));
        dy1 = -sz * ((-x6 * vec[0][2] + x2 * (z3 - z6) + x3 * vec[2][2] - (vec[1][1] * (-(x6 * vec[1][1] - x1 * vec[2][1] + x2 * (-y1 + y6)) * vec[0][2] - vec[0][0] * (y6 * vec[1][2] - y1 * vec[2][2] + y2 * (-z1 + z6)) -
            vec[0][1] * (-x6 * vec[1][2] + x2 * (z1 - z6) + x1 * vec[2][2]))) / (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2])) /
            sqrt((vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])));
        dz1 = -sz * ((x6 * vec[0][1] - x3 * vec[2][1] + x2 * (-y3 + y6) - (vec[1][2] * (-(x6 * vec[1][1] - x1 * vec[2][1] + x2 * (-y1 + y6)) * vec[0][2] - vec[0][0] * (y6 * vec[1][2] - y1 * vec[2][2] + y2 * (-z1 + z6)) -
            vec[0][1] * (-x6 * vec[1][2] + x2 * (z1 - z6) + x1 * vec[2][2]))) / (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2])) /
            sqrt((vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])));
        dx2 = -sz * ((y6 * (z1 - z3) + y1 * (z3 - z6) + y3 * (-z1 + z6)) /
            sqrt((vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]))
            - ((-vec[2][0] * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) -
                vec[0][0] * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) -
                vec[1][0] * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])) *
                (-(x6 * vec[1][1] - x1 * vec[2][1] + x2 * (-y1 + y6)) * vec[0][2] - vec[0][0] * (y6 * vec[1][2] - y1 * vec[2][2] + y2 * (-z1 + z6)) - vec[0][1] * (-x6 * vec[1][2] + x2 * (z1 - z6) + x1 * vec[2][2]))) /
            pow((vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]),
                1.5));
        dy2 = -sz * ((x6 * (-z1 + z3) + x3 * (z1 - z6) + x1 * (-z3 + z6)) /
            sqrt((vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]))
            - ((-vec[2][1] * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) -
                vec[0][1] * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) -
                vec[1][1] * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])) *
                (-(x6 * vec[1][1] - x1 * vec[2][1] + x2 * (-y1 + y6)) * vec[0][2] - vec[0][0] * (y6 * vec[1][2] - y1 * vec[2][2] + y2 * (-z1 + z6)) - vec[0][1] * (-x6 * vec[1][2] + x2 * (z1 - z6) + x1 * vec[2][2]))) /
            pow((vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]),
                1.5));
        dz2 = -sz * ((x6 * (y1 - y3) + x1 * (y3 - y6) + x3 * (-y1 + y6)) /
            sqrt((vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]))
            - ((-vec[1][2] * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) -
                (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * vec[0][2] * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) -
                (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * vec[2][2]) *
                (-(x6 * vec[1][1] - x1 * vec[2][1] + x2 * (-y1 + y6)) * vec[0][2] - vec[0][0] * (y6 * vec[1][2] - y1 * vec[2][2] + y2 * (-z1 + z6)) - vec[0][1] * (-x6 * vec[1][2] + x2 * (z1 - z6) + x1 * vec[2][2]))) /
            pow((vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]),
                1.5));
        dx3 = -sz * ((-y6 * vec[1][2] + y2 * (z1 - z6) + y1 * vec[2][2] - (vec[0][0] * (-(x6 * vec[1][1] - x1 * vec[2][1] + x2 * (-y1 + y6)) * vec[0][2] - vec[0][0] * (y6 * vec[1][2] - y1 * vec[2][2] + y2 * (-z1 + z6)) -
            vec[0][1] * (-x6 * vec[1][2] + x2 * (z1 - z6) + x1 * vec[2][2]))) / (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2])) /
            sqrt((vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])));
        dy3 = -sz * ((x6 * vec[1][2] - x1 * vec[2][2] + x2 * (-z1 + z6) - (vec[0][1] * (-(x6 * vec[1][1] - x1 * vec[2][1] + x2 * (-y1 + y6)) * vec[0][2] - vec[0][0] * (y6 * vec[1][2] - y1 * vec[2][2] + y2 * (-z1 + z6)) -
            vec[0][1] * (-x6 * vec[1][2] + x2 * (z1 - z6) + x1 * vec[2][2]))) / (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2])) /
            sqrt((vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])));
        dz3 = -sz * ((-x6 * vec[1][1] + x2 * (y1 - y6) + x1 * vec[2][1] - (vec[0][2] * (-(x6 * vec[1][1] - x1 * vec[2][1] + x2 * (-y1 + y6)) * vec[0][2] - vec[0][0] * (y6 * vec[1][2] - y1 * vec[2][2] + y2 * (-z1 + z6)) -
            vec[0][1] * (-x6 * vec[1][2] + x2 * (z1 - z6) + x1 * vec[2][2]))) / (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2])) /
            sqrt((vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])));
        dx6 = -sz * ((y3 * vec[1][2] - y1 * vec[0][2] + y2 * (-z1 + z3) - (vec[2][0] * (-(x6 * vec[1][1] - x1 * vec[2][1] + x2 * (-y1 + y6)) * vec[0][2] - vec[0][0] * (y6 * vec[1][2] - y1 * vec[2][2] + y2 * (-z1 + z6)) -
            vec[0][1] * (-x6 * vec[1][2] + x2 * (z1 - z6) + x1 * vec[2][2]))) / (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])) /
            sqrt((vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])));
        dy6 = -sz * ((-x3 * vec[1][2] + x2 * (z1 - z3) + x1 * vec[0][2] - (vec[2][1] * (-(x6 * vec[1][1] - x1 * vec[2][1] + x2 * (-y1 + y6)) * vec[0][2] - vec[0][0] * (y6 * vec[1][2] - y1 * vec[2][2] + y2 * (-z1 + z6)) -
            vec[0][1] * (-x6 * vec[1][2] + x2 * (z1 - z6) + x1 * vec[2][2]))) / (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])) /
            sqrt((vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])));
        dz6 = -sz * ((x3 * vec[1][1] - x1 * vec[0][1] + x2 * (-y1 + y3) - (vec[2][2] * (-(x6 * vec[1][1] - x1 * vec[2][1] + x2 * (-y1 + y6)) * vec[0][2] - vec[0][0] * (y6 * vec[1][2] - y1 * vec[2][2] + y2 * (-z1 + z6)) -
            vec[0][1] * (-x6 * vec[1][2] + x2 * (z1 - z6) + x1 * vec[2][2]))) / (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])) /
            sqrt((vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])));
        return;
    }

    // 3 0 2 7
    vec[0][0] = x0 - x3; vec[0][1] = y0 - y3; vec[0][2] = z0 - z3;
    vec[1][0] = x2 - x3; vec[1][1] = y2 - y3; vec[1][2] = z2 - z3;
    vec[2][0] = x7 - x3; vec[2][1] = y7 - y3; vec[2][2] = z7 - z3;
    vol = vec[0][0] * (vec[2][1] * vec[1][2] - vec[1][1] * vec[2][2]) +
        vec[0][1] * (vec[2][2] * vec[1][0] - vec[1][2] * vec[2][0]) +
        vec[0][2] * (vec[2][0] * vec[1][1] - vec[1][0] * vec[2][1]);
    len3 = sqrt(DOT(vec[0], vec[0]) * DOT(vec[1], vec[1]) * DOT(vec[2], vec[2]));
    if (vol >= 0) {
        eFSJ = false;
        dx0 = -sz2 * (-y7 * vec[1][2] + y3 * (z2 - z7) + y2 * vec[2][2]);
        dy0 = -sz2 * (x7 * vec[1][2] - x2 * vec[2][2] + x3 * (-z2 + z7));
        dz0 = -sz2 * (-x7 * vec[1][1] + x3 * (y2 - y7) + x2 * vec[2][1]);
        dx2 = -sz2 * (y7 * vec[0][2] - y0 * vec[2][2] + y3 * (-z0 + z7));
        dy2 = -sz2 * (-x7 * vec[0][2] + x3 * (z0 - z7) + x0 * vec[2][2]);
        dz2 = -sz2 * (x7 * vec[0][1] - x0 * vec[2][1] + x3 * (-y0 + y7));
        dx3 = -sz2 * (y7 * (-z0 + z2) + y2 * (z0 - z7) + y0 * (-z2 + z7));
        dy3 = -sz2 * (x7 * (z0 - z2) + x0 * (z2 - z7) + x2 * (-z0 + z7));
        dz3 = -sz2 * (x7 * (-y0 + y2) + x2 * (y0 - y7) + x0 * (-y2 + y7));
        dx7 = -sz2 * (y3 * (z0 - z2) + y0 * vec[1][2] - y2 * vec[0][2]);
        dy7 = -sz2 * (x3 * (-z0 + z2) + x2 * vec[0][2] - x0 * vec[1][2]);
        dz7 = -sz2 * (x3 * (y0 - y2) + x0 * vec[1][1] - x2 * vec[0][1]);
        return;
    }
    else if (vol > len3 * sJThres) {
        eFSJ = false;
        dx0 = -sz * ((-y7 * vec[1][2] + y3 * (z2 - z7) + y2 * vec[2][2] - (vec[0][0] * ((-x7 * vec[1][1] + x3 * (y2 - y7) + x2 * vec[2][1]) * vec[0][2] + vec[0][1] * (x7 * vec[1][2] - x2 * vec[2][2] + x3 * (-z2 + z7)) +
            vec[0][0] * (-y7 * vec[1][2] + y3 * (z2 - z7) + y2 * vec[2][2]))) / (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2])) /
            sqrt((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])));
        dy0 = -sz * ((x7 * vec[1][2] - x2 * vec[2][2] + x3 * (-z2 + z7) - (vec[0][1] * ((-x7 * vec[1][1] + x3 * (y2 - y7) + x2 * vec[2][1]) * vec[0][2] + vec[0][1] * (x7 * vec[1][2] - x2 * vec[2][2] + x3 * (-z2 + z7)) +
            vec[0][0] * (-y7 * vec[1][2] + y3 * (z2 - z7) + y2 * vec[2][2]))) / (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2])) /
            sqrt((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])));
        dz0 = -sz * ((-x7 * vec[1][1] + x3 * (y2 - y7) + x2 * vec[2][1] - (vec[0][2] * ((-x7 * vec[1][1] + x3 * (y2 - y7) + x2 * vec[2][1]) * vec[0][2] + vec[0][1] * (x7 * vec[1][2] - x2 * vec[2][2] + x3 * (-z2 + z7)) +
            vec[0][0] * (-y7 * vec[1][2] + y3 * (z2 - z7) + y2 * vec[2][2]))) / (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2])) /
            sqrt((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])));
        dx2 = -sz * ((y7 * vec[0][2] - y0 * vec[2][2] + y3 * (-z0 + z7) - (vec[1][0] * ((-x7 * vec[1][1] + x3 * (y2 - y7) + x2 * vec[2][1]) * vec[0][2] + vec[0][1] * (x7 * vec[1][2] - x2 * vec[2][2] + x3 * (-z2 + z7)) +
            vec[0][0] * (-y7 * vec[1][2] + y3 * (z2 - z7) + y2 * vec[2][2]))) / (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2])) /
            sqrt((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])));
        dy2 = -sz * ((-x7 * vec[0][2] + x3 * (z0 - z7) + x0 * vec[2][2] - (vec[1][1] * ((-x7 * vec[1][1] + x3 * (y2 - y7) + x2 * vec[2][1]) * vec[0][2] + vec[0][1] * (x7 * vec[1][2] - x2 * vec[2][2] + x3 * (-z2 + z7)) +
            vec[0][0] * (-y7 * vec[1][2] + y3 * (z2 - z7) + y2 * vec[2][2]))) / (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2])) /
            sqrt((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])));
        dz2 = -sz * ((x7 * vec[0][1] - x0 * vec[2][1] + x3 * (-y0 + y7) - (vec[1][2] * ((-x7 * vec[1][1] + x3 * (y2 - y7) + x2 * vec[2][1]) * vec[0][2] + vec[0][1] * (x7 * vec[1][2] - x2 * vec[2][2] + x3 * (-z2 + z7)) +
            vec[0][0] * (-y7 * vec[1][2] + y3 * (z2 - z7) + y2 * vec[2][2]))) / (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2])) /
            sqrt((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])));
        dx3 = -sz * ((y7 * (-z0 + z2) + y2 * (z0 - z7) + y0 * (-z2 + z7)) /
            sqrt((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]))
            - ((-vec[2][0] * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) -
                vec[1][0] * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) -
                vec[0][0] * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])) *
                ((-x7 * vec[1][1] + x3 * (y2 - y7) + x2 * vec[2][1]) * vec[0][2] + vec[0][1] * (x7 * vec[1][2] - x2 * vec[2][2] + x3 * (-z2 + z7)) + vec[0][0] * (-y7 * vec[1][2] + y3 * (z2 - z7) + y2 * vec[2][2]))) /
            pow((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]),
                1.5));
        dy3 = -sz * ((x7 * (z0 - z2) + x0 * (z2 - z7) + x2 * (-z0 + z7)) /
            sqrt((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]))
            - ((-vec[2][1] * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) -
                vec[1][1] * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) -
                vec[0][1] * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])) *
                ((-x7 * vec[1][1] + x3 * (y2 - y7) + x2 * vec[2][1]) * vec[0][2] + vec[0][1] * (x7 * vec[1][2] - x2 * vec[2][2] + x3 * (-z2 + z7)) + vec[0][0] * (-y7 * vec[1][2] + y3 * (z2 - z7) + y2 * vec[2][2]))) /
            pow((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]),
                1.5));
        dz3 = -sz * ((x7 * (-y0 + y2) + x2 * (y0 - y7) + x0 * (-y2 + y7)) /
            sqrt((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]))
            - ((-(vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * vec[0][2] * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) -
                (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * vec[1][2] * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) -
                (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * vec[2][2]) *
                ((-x7 * vec[1][1] + x3 * (y2 - y7) + x2 * vec[2][1]) * vec[0][2] + vec[0][1] * (x7 * vec[1][2] - x2 * vec[2][2] + x3 * (-z2 + z7)) + vec[0][0] * (-y7 * vec[1][2] + y3 * (z2 - z7) + y2 * vec[2][2]))) /
            pow((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]),
                1.5));
        dx7 = -sz * ((y3 * (z0 - z2) + y0 * vec[1][2] - y2 * vec[0][2] - (vec[2][0] * ((-x7 * vec[1][1] + x3 * (y2 - y7) + x2 * vec[2][1]) * vec[0][2] + vec[0][1] * (x7 * vec[1][2] - x2 * vec[2][2] + x3 * (-z2 + z7)) +
            vec[0][0] * (-y7 * vec[1][2] + y3 * (z2 - z7) + y2 * vec[2][2]))) / (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])) /
            sqrt((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])));
        dy7 = -sz * ((x3 * (-z0 + z2) + x2 * vec[0][2] - x0 * vec[1][2] - (vec[2][1] * ((-x7 * vec[1][1] + x3 * (y2 - y7) + x2 * vec[2][1]) * vec[0][2] + vec[0][1] * (x7 * vec[1][2] - x2 * vec[2][2] + x3 * (-z2 + z7)) +
            vec[0][0] * (-y7 * vec[1][2] + y3 * (z2 - z7) + y2 * vec[2][2]))) / (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])) /
            sqrt((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])));
        dz7 = -sz * ((x3 * (y0 - y2) + x0 * vec[1][1] - x2 * vec[0][1] - (vec[2][2] * ((-x7 * vec[1][1] + x3 * (y2 - y7) + x2 * vec[2][1]) * vec[0][2] + vec[0][1] * (x7 * vec[1][2] - x2 * vec[2][2] + x3 * (-z2 + z7)) +
            vec[0][0] * (-y7 * vec[1][2] + y3 * (z2 - z7) + y2 * vec[2][2]))) / (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])) /
            sqrt((vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])));
        return;
    }

    // 4 7 5 0
    vec[0][0] = x7 - x4; vec[0][1] = y7 - y4; vec[0][2] = z7 - z4;
    vec[1][0] = x5 - x4; vec[1][1] = y5 - y4; vec[1][2] = z5 - z4;
    vec[2][0] = x0 - x4; vec[2][1] = y0 - y4; vec[2][2] = z0 - z4;
    vol = vec[0][0] * (vec[2][1] * vec[1][2] - vec[1][1] * vec[2][2]) +
        vec[0][1] * (vec[2][2] * vec[1][0] - vec[1][2] * vec[2][0]) +
        vec[0][2] * (vec[2][0] * vec[1][1] - vec[1][0] * vec[2][1]);
    len3 = sqrt(DOT(vec[0], vec[0]) * DOT(vec[1], vec[1]) * DOT(vec[2], vec[2]));
    if (vol >= 0) {
        eFSJ = false;
        dx0 = -sz2 * (y7 * vec[1][2] - y5 * vec[0][2] + y4 * (-z5 + z7));
        dy0 = -sz2 * (-x7 * vec[1][2] + x4 * (z5 - z7) + x5 * vec[0][2]);
        dz0 = -sz2 * (x7 * vec[1][1] - x5 * vec[0][1] + x4 * (-y5 + y7));
        dx4 = -sz2 * (y7 * (z0 - z5) + y0 * (z5 - z7) + y5 * (-z0 + z7));
        dy4 = -sz2 * (x7 * (-z0 + z5) + x5 * (z0 - z7) + x0 * (-z5 + z7));
        dz4 = -sz2 * (x7 * (y0 - y5) + x0 * (y5 - y7) + x5 * (-y0 + y7));
        dx5 = -sz2 * (-y7 * vec[2][2] + y4 * (z0 - z7) + y0 * vec[0][2]);
        dy5 = -sz2 * (x7 * vec[2][2] - x0 * vec[0][2] + x4 * (-z0 + z7));
        dz5 = -sz2 * (-x7 * vec[2][1] + x4 * (y0 - y7) + x0 * vec[0][1]);
        dx7 = -sz2 * (y5 * vec[2][2] - y0 * vec[1][2] + y4 * (-z0 + z5));
        dy7 = -sz2 * (-x5 * vec[2][2] + x4 * (z0 - z5) + x0 * vec[1][2]);
        dz7 = -sz2 * (x5 * vec[2][1] - x0 * vec[1][1] + x4 * (-y0 + y5));
        return;
    }
    else if (vol > len3 * sJThres) {
        eFSJ = false;
        dx0 = -sz * ((y7 * vec[1][2] - (vec[2][0] * (-vec[0][1] * (x5 * vec[2][2] - x0 * vec[1][2] + x4 * (-z0 + z5)) - vec[0][0] * (-y5 * vec[2][2] + y4 * (z0 - z5) + y0 * vec[1][2]) -
            (-x5 * vec[2][1] + x4 * (y0 - y5) + x0 * vec[1][1]) * vec[0][2])) / (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) - y5 * vec[0][2] + y4 * (-z5 + z7)) /
            sqrt((vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2])));
        dy0 = -sz * ((-x7 * vec[1][2] - (vec[2][1] * (-vec[0][1] * (x5 * vec[2][2] - x0 * vec[1][2] + x4 * (-z0 + z5)) - vec[0][0] * (-y5 * vec[2][2] + y4 * (z0 - z5) + y0 * vec[1][2]) -
            (-x5 * vec[2][1] + x4 * (y0 - y5) + x0 * vec[1][1]) * vec[0][2])) / (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) + x4 * (z5 - z7) + x5 * vec[0][2]) /
            sqrt((vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2])));
        dz0 = -sz * ((x7 * vec[1][1] - x5 * vec[0][1] + x4 * (-y5 + y7) - (vec[2][2] * (-vec[0][1] * (x5 * vec[2][2] - x0 * vec[1][2] + x4 * (-z0 + z5)) - vec[0][0] * (-y5 * vec[2][2] + y4 * (z0 - z5) + y0 * vec[1][2]) -
            (-x5 * vec[2][1] + x4 * (y0 - y5) + x0 * vec[1][1]) * vec[0][2])) / (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])) /
            sqrt((vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2])));
        dx4 = -sz * (((vec[0][0] * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) +
            vec[1][0] * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) +
            vec[2][0] * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2])) *
            (-vec[0][1] * (x5 * vec[2][2] - x0 * vec[1][2] + x4 * (-z0 + z5)) - vec[0][0] * (-y5 * vec[2][2] + y4 * (z0 - z5) + y0 * vec[1][2]) - (-x5 * vec[2][1] + x4 * (y0 - y5) + x0 * vec[1][1]) * vec[0][2])) /
            pow((vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]),
                1.5) + (y7 * (z0 - z5) + y0 * (z5 - z7) + y5 * (-z0 + z7)) /
            sqrt((vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2])));
        dy4 = -sz * (((vec[0][1] * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) +
            vec[1][1] * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) +
            vec[2][1] * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2])) *
            (-vec[0][1] * (x5 * vec[2][2] - x0 * vec[1][2] + x4 * (-z0 + z5)) - vec[0][0] * (-y5 * vec[2][2] + y4 * (z0 - z5) + y0 * vec[1][2]) - (-x5 * vec[2][1] + x4 * (y0 - y5) + x0 * vec[1][1]) * vec[0][2])) /
            pow((vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]),
                1.5) + (x7 * (-z0 + z5) + x5 * (z0 - z7) + x0 * (-z5 + z7)) /
            sqrt((vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2])));
        dz4 = -sz * ((x7 * (y0 - y5) + x0 * (y5 - y7) + x5 * (-y0 + y7)) /
            sqrt((vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2])) -
            ((-vec[0][1] * (x5 * vec[2][2] - x0 * vec[1][2] + x4 * (-z0 + z5)) - vec[0][0] * (-y5 * vec[2][2] + y4 * (z0 - z5) + y0 * vec[1][2]) -
                (-x5 * vec[2][1] + x4 * (y0 - y5) + x0 * vec[1][1]) * vec[0][2]) * (-vec[2][2] * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) *
                    (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) - (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * vec[1][2] *
                    (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) - (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) *
                    vec[0][2])) / pow((vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) *
                        (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]), 1.5));
        dx5 = -sz * ((-y7 * vec[2][2] - (vec[1][0] * (-vec[0][1] * (x5 * vec[2][2] - x0 * vec[1][2] + x4 * (-z0 + z5)) - vec[0][0] * (-y5 * vec[2][2] + y4 * (z0 - z5) + y0 * vec[1][2]) -
            (-x5 * vec[2][1] + x4 * (y0 - y5) + x0 * vec[1][1]) * vec[0][2])) / (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) + y4 * (z0 - z7) + y0 * vec[0][2]) /
            sqrt((vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2])));
        dy5 = -sz * ((x7 * vec[2][2] - (vec[1][1] * (-vec[0][1] * (x5 * vec[2][2] - x0 * vec[1][2] + x4 * (-z0 + z5)) - vec[0][0] * (-y5 * vec[2][2] + y4 * (z0 - z5) + y0 * vec[1][2]) -
            (-x5 * vec[2][1] + x4 * (y0 - y5) + x0 * vec[1][1]) * vec[0][2])) / (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) - x0 * vec[0][2] + x4 * (-z0 + z7)) /
            sqrt((vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2])));
        dz5 = -sz * ((-x7 * vec[2][1] + x4 * (y0 - y7) + x0 * vec[0][1] - (vec[1][2] * (-vec[0][1] * (x5 * vec[2][2] - x0 * vec[1][2] + x4 * (-z0 + z5)) - vec[0][0] * (-y5 * vec[2][2] + y4 * (z0 - z5) + y0 * vec[1][2]) -
            (-x5 * vec[2][1] + x4 * (y0 - y5) + x0 * vec[1][1]) * vec[0][2])) / (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2])) /
            sqrt((vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2])));
        dx7 = -sz * ((y5 * vec[2][2] - y0 * vec[1][2] + y4 * (-z0 + z5) - (vec[0][0] * (-vec[0][1] * (x5 * vec[2][2] - x0 * vec[1][2] + x4 * (-z0 + z5)) - vec[0][0] * (-y5 * vec[2][2] + y4 * (z0 - z5) + y0 * vec[1][2]) -
            (-x5 * vec[2][1] + x4 * (y0 - y5) + x0 * vec[1][1]) * vec[0][2])) / (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2])) /
            sqrt((vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2])));
        dy7 = -sz * ((-x5 * vec[2][2] + x4 * (z0 - z5) + x0 * vec[1][2] - (vec[0][1] * (-vec[0][1] * (x5 * vec[2][2] - x0 * vec[1][2] + x4 * (-z0 + z5)) - vec[0][0] * (-y5 * vec[2][2] + y4 * (z0 - z5) + y0 * vec[1][2]) -
            (-x5 * vec[2][1] + x4 * (y0 - y5) + x0 * vec[1][1]) * vec[0][2])) / (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2])) /
            sqrt((vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2])));
        dz7 = -sz * ((x5 * vec[2][1] - x0 * vec[1][1] + x4 * (-y0 + y5) - ((-vec[0][1] * (x5 * vec[2][2] - x0 * vec[1][2] + x4 * (-z0 + z5)) - vec[0][0] * (-y5 * vec[2][2] + y4 * (z0 - z5) + y0 * vec[1][2]) -
            (-x5 * vec[2][1] + x4 * (y0 - y5) + x0 * vec[1][1]) * vec[0][2]) * vec[0][2]) / (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2])) /
            sqrt((vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2])));
        return;
    }

    // 5 4 6 1
    vec[0][0] = x4 - x5; vec[0][1] = y4 - y5; vec[0][2] = z4 - z5;
    vec[1][0] = x6 - x5; vec[1][1] = y6 - y5; vec[1][2] = z6 - z5;
    vec[2][0] = x1 - x5; vec[2][1] = y1 - y5; vec[2][2] = z1 - z5;
    vol = vec[0][0] * (vec[2][1] * vec[1][2] - vec[1][1] * vec[2][2]) +
        vec[0][1] * (vec[2][2] * vec[1][0] - vec[1][2] * vec[2][0]) +
        vec[0][2] * (vec[2][0] * vec[1][1] - vec[1][0] * vec[2][1]);
    len3 = sqrt(DOT(vec[0], vec[0]) * DOT(vec[1], vec[1]) * DOT(vec[2], vec[2]));
    if (vol >= 0) {
        eFSJ = false;
        dx1 = -sz2 * (-y6 * vec[0][2] + y5 * (z4 - z6) + y4 * vec[1][2]);
        dy1 = -sz2 * (x6 * vec[0][2] - x4 * vec[1][2] + x5 * (-z4 + z6));
        dz1 = -sz2 * (-x6 * vec[0][1] + x5 * (y4 - y6) + x4 * vec[1][1]);
        dx4 = -sz2 * (y6 * vec[2][2] - y1 * vec[1][2] + y5 * (-z1 + z6));
        dy4 = -sz2 * (-x6 * vec[2][2] + x5 * (z1 - z6) + x1 * vec[1][2]);
        dz4 = -sz2 * (x6 * vec[2][1] - x1 * vec[1][1] + x5 * (-y1 + y6));
        dx5 = -sz2 * (y6 * (-z1 + z4) + y4 * (z1 - z6) + y1 * (-z4 + z6));
        dy5 = -sz2 * (x6 * (z1 - z4) + x1 * (z4 - z6) + x4 * (-z1 + z6));
        dz5 = -sz2 * (x6 * (-y1 + y4) + x4 * (y1 - y6) + x1 * (-y4 + y6));
        dx6 = -sz2 * (y5 * (z1 - z4) + y1 * vec[0][2] - y4 * vec[2][2]);
        dy6 = -sz2 * (x5 * (-z1 + z4) + x4 * vec[2][2] - x1 * vec[0][2]);
        dz6 = -sz2 * (x5 * (y1 - y4) + x1 * vec[0][1] - x4 * vec[2][1]);
        return;
    }
    else if (vol > len3 * sJThres) {
        eFSJ = false;
        dx1 = -sz * ((-y6 * vec[0][2] + y5 * (z4 - z6) + y4 * vec[1][2] - (vec[2][0] * ((x6 * vec[2][1] - x1 * vec[1][1] + x5 * (-y1 + y6)) * vec[0][2] + vec[0][0] * (y6 * vec[2][2] - y1 * vec[1][2] + y5 * (-z1 + z6)) +
            vec[0][1] * (-x6 * vec[2][2] + x5 * (z1 - z6) + x1 * vec[1][2]))) / (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])) /
            sqrt((vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2])));
        dy1 = -sz * ((x6 * vec[0][2] - x4 * vec[1][2] + x5 * (-z4 + z6) - (vec[2][1] * ((x6 * vec[2][1] - x1 * vec[1][1] + x5 * (-y1 + y6)) * vec[0][2] + vec[0][0] * (y6 * vec[2][2] - y1 * vec[1][2] + y5 * (-z1 + z6)) +
            vec[0][1] * (-x6 * vec[2][2] + x5 * (z1 - z6) + x1 * vec[1][2]))) / (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])) /
            sqrt((vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2])));
        dz1 = -sz * ((-x6 * vec[0][1] + x5 * (y4 - y6) + x4 * vec[1][1] - (vec[2][2] * ((x6 * vec[2][1] - x1 * vec[1][1] + x5 * (-y1 + y6)) * vec[0][2] + vec[0][0] * (y6 * vec[2][2] - y1 * vec[1][2] + y5 * (-z1 + z6)) +
            vec[0][1] * (-x6 * vec[2][2] + x5 * (z1 - z6) + x1 * vec[1][2]))) / (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])) /
            sqrt((vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2])));
        dx4 = -sz * ((y6 * vec[2][2] - y1 * vec[1][2] + y5 * (-z1 + z6) - (vec[0][0] * ((x6 * vec[2][1] - x1 * vec[1][1] + x5 * (-y1 + y6)) * vec[0][2] + vec[0][0] * (y6 * vec[2][2] - y1 * vec[1][2] + y5 * (-z1 + z6)) +
            vec[0][1] * (-x6 * vec[2][2] + x5 * (z1 - z6) + x1 * vec[1][2]))) / (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2])) /
            sqrt((vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2])));
        dy4 = -sz * ((-x6 * vec[2][2] + x5 * (z1 - z6) + x1 * vec[1][2] - (vec[0][1] * ((x6 * vec[2][1] - x1 * vec[1][1] + x5 * (-y1 + y6)) * vec[0][2] + vec[0][0] * (y6 * vec[2][2] - y1 * vec[1][2] + y5 * (-z1 + z6)) +
            vec[0][1] * (-x6 * vec[2][2] + x5 * (z1 - z6) + x1 * vec[1][2]))) / (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2])) /
            sqrt((vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2])));
        dz4 = -sz * ((x6 * vec[2][1] - x1 * vec[1][1] + x5 * (-y1 + y6) - (vec[0][2] * ((x6 * vec[2][1] - x1 * vec[1][1] + x5 * (-y1 + y6)) * vec[0][2] + vec[0][0] * (y6 * vec[2][2] - y1 * vec[1][2] + y5 * (-z1 + z6)) +
            vec[0][1] * (-x6 * vec[2][2] + x5 * (z1 - z6) + x1 * vec[1][2]))) / (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2])) /
            sqrt((vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2])));
        dx5 = -sz * ((y6 * (-z1 + z4) + y4 * (z1 - z6) + y1 * (-z4 + z6)) /
            sqrt((vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]))
            - ((-vec[1][0] * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) -
                vec[0][0] * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) -
                vec[2][0] * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2])) *
                ((x6 * vec[2][1] - x1 * vec[1][1] + x5 * (-y1 + y6)) * vec[0][2] + vec[0][0] * (y6 * vec[2][2] - y1 * vec[1][2] + y5 * (-z1 + z6)) + vec[0][1] * (-x6 * vec[2][2] + x5 * (z1 - z6) + x1 * vec[1][2]))) /
            pow((vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]),
                1.5));
        dy5 = -sz * ((x6 * (z1 - z4) + x1 * (z4 - z6) + x4 * (-z1 + z6)) /
            sqrt((vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]))
            - ((-vec[1][1] * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) -
                vec[0][1] * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) -
                vec[2][1] * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2])) *
                ((x6 * vec[2][1] - x1 * vec[1][1] + x5 * (-y1 + y6)) * vec[0][2] + vec[0][0] * (y6 * vec[2][2] - y1 * vec[1][2] + y5 * (-z1 + z6)) + vec[0][1] * (-x6 * vec[2][2] + x5 * (z1 - z6) + x1 * vec[1][2]))) /
            pow((vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]),
                1.5));
        dz5 = -sz * ((x6 * (-y1 + y4) + x4 * (y1 - y6) + x1 * (-y4 + y6)) /
            sqrt((vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]))
            - ((-(vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * vec[2][2] * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) -
                (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * vec[0][2] * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) -
                (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * vec[1][2]) *
                ((x6 * vec[2][1] - x1 * vec[1][1] + x5 * (-y1 + y6)) * vec[0][2] + vec[0][0] * (y6 * vec[2][2] - y1 * vec[1][2] + y5 * (-z1 + z6)) + vec[0][1] * (-x6 * vec[2][2] + x5 * (z1 - z6) + x1 * vec[1][2]))) /
            pow((vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]),
                1.5));
        dx6 = -sz * ((y5 * (z1 - z4) + y1 * vec[0][2] - y4 * vec[2][2] - (vec[1][0] * ((x6 * vec[2][1] - x1 * vec[1][1] + x5 * (-y1 + y6)) * vec[0][2] + vec[0][0] * (y6 * vec[2][2] - y1 * vec[1][2] + y5 * (-z1 + z6)) +
            vec[0][1] * (-x6 * vec[2][2] + x5 * (z1 - z6) + x1 * vec[1][2]))) / (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2])) /
            sqrt((vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2])));
        dy6 = -sz * ((x5 * (-z1 + z4) + x4 * vec[2][2] - x1 * vec[0][2] - (vec[1][1] * ((x6 * vec[2][1] - x1 * vec[1][1] + x5 * (-y1 + y6)) * vec[0][2] + vec[0][0] * (y6 * vec[2][2] - y1 * vec[1][2] + y5 * (-z1 + z6)) +
            vec[0][1] * (-x6 * vec[2][2] + x5 * (z1 - z6) + x1 * vec[1][2]))) / (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2])) /
            sqrt((vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2])));
        dz6 = -sz * ((x5 * (y1 - y4) + x1 * vec[0][1] - x4 * vec[2][1] - (vec[1][2] * ((x6 * vec[2][1] - x1 * vec[1][1] + x5 * (-y1 + y6)) * vec[0][2] + vec[0][0] * (y6 * vec[2][2] - y1 * vec[1][2] + y5 * (-z1 + z6)) +
            vec[0][1] * (-x6 * vec[2][2] + x5 * (z1 - z6) + x1 * vec[1][2]))) / (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2])) /
            sqrt((vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2])));
        return;
    }

    // 6 5 7 2
    vec[0][0] = x5 - x6; vec[0][1] = y5 - y6; vec[0][2] = z5 - z6;
    vec[1][0] = x7 - x6; vec[1][1] = y7 - y6; vec[1][2] = z7 - z6;
    vec[2][0] = x2 - x6; vec[2][1] = y2 - y6; vec[2][2] = z2 - z6;
    vol = vec[0][0] * (vec[2][1] * vec[1][2] - vec[1][1] * vec[2][2]) +
        vec[0][1] * (vec[2][2] * vec[1][0] - vec[1][2] * vec[2][0]) +
        vec[0][2] * (vec[2][0] * vec[1][1] - vec[1][0] * vec[2][1]);
    len3 = sqrt(DOT(vec[0], vec[0]) * DOT(vec[1], vec[1]) * DOT(vec[2], vec[2]));
    if (vol >= 0) {
        eFSJ = false;
        dx2 = -sz2 * (-y7 * vec[0][2] + y6 * (z5 - z7) + y5 * vec[1][2]);
        dy2 = -sz2 * (x7 * vec[0][2] - x5 * vec[1][2] + x6 * (-z5 + z7));
        dz2 = -sz2 * (-x7 * vec[0][1] + x6 * (y5 - y7) + x5 * vec[1][1]);
        dx5 = -sz2 * (y7 * vec[2][2] - y2 * vec[1][2] + y6 * (-z2 + z7));
        dy5 = -sz2 * (-x7 * vec[2][2] + x6 * (z2 - z7) + x2 * vec[1][2]);
        dz5 = -sz2 * (x7 * vec[2][1] - x2 * vec[1][1] + x6 * (-y2 + y7));
        dx6 = -sz2 * (y7 * (-z2 + z5) + y5 * (z2 - z7) + y2 * (-z5 + z7));
        dy6 = -sz2 * (x7 * (z2 - z5) + x2 * (z5 - z7) + x5 * (-z2 + z7));
        dz6 = -sz2 * (x7 * (-y2 + y5) + x5 * (y2 - y7) + x2 * (-y5 + y7));
        dx7 = -sz2 * (y6 * (z2 - z5) + y2 * vec[0][2] - y5 * vec[2][2]);
        dy7 = -sz2 * (x6 * (-z2 + z5) + x5 * vec[2][2] - x2 * vec[0][2]);
        dz7 = -sz2 * (x6 * (y2 - y5) + x2 * vec[0][1] - x5 * vec[2][1]);
        return;
    }
    else if (vol > len3 * sJThres) {
        eFSJ = false;
        dx2 = -sz * ((-y7 * vec[0][2] + y6 * (z5 - z7) + y5 * vec[1][2] - (vec[2][0] * ((x7 * vec[2][1] - x2 * vec[1][1] + x6 * (-y2 + y7)) * vec[0][2] + vec[0][0] * (y7 * vec[2][2] - y2 * vec[1][2] + y6 * (-z2 + z7)) +
            vec[0][1] * (-x7 * vec[2][2] + x6 * (z2 - z7) + x2 * vec[1][2]))) / (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])) /
            sqrt((vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2])));
        dy2 = -sz * ((x7 * vec[0][2] - x5 * vec[1][2] + x6 * (-z5 + z7) - (vec[2][1] * ((x7 * vec[2][1] - x2 * vec[1][1] + x6 * (-y2 + y7)) * vec[0][2] + vec[0][0] * (y7 * vec[2][2] - y2 * vec[1][2] + y6 * (-z2 + z7)) +
            vec[0][1] * (-x7 * vec[2][2] + x6 * (z2 - z7) + x2 * vec[1][2]))) / (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])) /
            sqrt((vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2])));
        dz2 = -sz * ((-x7 * vec[0][1] + x6 * (y5 - y7) + x5 * vec[1][1] - (vec[2][2] * ((x7 * vec[2][1] - x2 * vec[1][1] + x6 * (-y2 + y7)) * vec[0][2] + vec[0][0] * (y7 * vec[2][2] - y2 * vec[1][2] + y6 * (-z2 + z7)) +
            vec[0][1] * (-x7 * vec[2][2] + x6 * (z2 - z7) + x2 * vec[1][2]))) / (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])) /
            sqrt((vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2])));
        dx5 = -sz * ((y7 * vec[2][2] - y2 * vec[1][2] + y6 * (-z2 + z7) - (vec[0][0] * ((x7 * vec[2][1] - x2 * vec[1][1] + x6 * (-y2 + y7)) * vec[0][2] + vec[0][0] * (y7 * vec[2][2] - y2 * vec[1][2] + y6 * (-z2 + z7)) +
            vec[0][1] * (-x7 * vec[2][2] + x6 * (z2 - z7) + x2 * vec[1][2]))) / (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2])) /
            sqrt((vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2])));
        dy5 = -sz * ((-x7 * vec[2][2] + x6 * (z2 - z7) + x2 * vec[1][2] - (vec[0][1] * ((x7 * vec[2][1] - x2 * vec[1][1] + x6 * (-y2 + y7)) * vec[0][2] + vec[0][0] * (y7 * vec[2][2] - y2 * vec[1][2] + y6 * (-z2 + z7)) +
            vec[0][1] * (-x7 * vec[2][2] + x6 * (z2 - z7) + x2 * vec[1][2]))) / (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2])) /
            sqrt((vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2])));
        dz5 = -sz * ((x7 * vec[2][1] - x2 * vec[1][1] + x6 * (-y2 + y7) - (vec[0][2] * ((x7 * vec[2][1] - x2 * vec[1][1] + x6 * (-y2 + y7)) * vec[0][2] + vec[0][0] * (y7 * vec[2][2] - y2 * vec[1][2] + y6 * (-z2 + z7)) +
            vec[0][1] * (-x7 * vec[2][2] + x6 * (z2 - z7) + x2 * vec[1][2]))) / (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2])) /
            sqrt((vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2])));
        dx6 = -sz * ((y7 * (-z2 + z5) + y5 * (z2 - z7) + y2 * (-z5 + z7)) /
            sqrt((vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]))
            - ((-vec[1][0] * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) -
                vec[0][0] * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) -
                vec[2][0] * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2])) *
                ((x7 * vec[2][1] - x2 * vec[1][1] + x6 * (-y2 + y7)) * vec[0][2] + vec[0][0] * (y7 * vec[2][2] - y2 * vec[1][2] + y6 * (-z2 + z7)) + vec[0][1] * (-x7 * vec[2][2] + x6 * (z2 - z7) + x2 * vec[1][2]))) /
            pow((vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]),
                1.5));
        dy6 = -sz * ((x7 * (z2 - z5) + x2 * (z5 - z7) + x5 * (-z2 + z7)) /
            sqrt((vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]))
            - ((-vec[1][1] * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) -
                vec[0][1] * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) -
                vec[2][1] * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2])) *
                ((x7 * vec[2][1] - x2 * vec[1][1] + x6 * (-y2 + y7)) * vec[0][2] + vec[0][0] * (y7 * vec[2][2] - y2 * vec[1][2] + y6 * (-z2 + z7)) + vec[0][1] * (-x7 * vec[2][2] + x6 * (z2 - z7) + x2 * vec[1][2]))) /
            pow((vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]),
                1.5));
        dz6 = -sz * ((x7 * (-y2 + y5) + x5 * (y2 - y7) + x2 * (-y5 + y7)) /
            sqrt((vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]))
            - ((-(vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * vec[2][2] * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) -
                (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * vec[0][2] * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) -
                (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * vec[1][2]) *
                ((x7 * vec[2][1] - x2 * vec[1][1] + x6 * (-y2 + y7)) * vec[0][2] + vec[0][0] * (y7 * vec[2][2] - y2 * vec[1][2] + y6 * (-z2 + z7)) + vec[0][1] * (-x7 * vec[2][2] + x6 * (z2 - z7) + x2 * vec[1][2]))) /
            pow((vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]),
                1.5));
        dx7 = -sz * ((y6 * (z2 - z5) + y2 * vec[0][2] - y5 * vec[2][2] - (vec[1][0] * ((x7 * vec[2][1] - x2 * vec[1][1] + x6 * (-y2 + y7)) * vec[0][2] + vec[0][0] * (y7 * vec[2][2] - y2 * vec[1][2] + y6 * (-z2 + z7)) +
            vec[0][1] * (-x7 * vec[2][2] + x6 * (z2 - z7) + x2 * vec[1][2]))) / (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2])) /
            sqrt((vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2])));
        dy7 = -sz * ((x6 * (-z2 + z5) + x5 * vec[2][2] - x2 * vec[0][2] - (vec[1][1] * ((x7 * vec[2][1] - x2 * vec[1][1] + x6 * (-y2 + y7)) * vec[0][2] + vec[0][0] * (y7 * vec[2][2] - y2 * vec[1][2] + y6 * (-z2 + z7)) +
            vec[0][1] * (-x7 * vec[2][2] + x6 * (z2 - z7) + x2 * vec[1][2]))) / (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2])) /
            sqrt((vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2])));
        dz7 = -sz * ((x6 * (y2 - y5) + x2 * vec[0][1] - x5 * vec[2][1] - (vec[1][2] * ((x7 * vec[2][1] - x2 * vec[1][1] + x6 * (-y2 + y7)) * vec[0][2] + vec[0][0] * (y7 * vec[2][2] - y2 * vec[1][2] + y6 * (-z2 + z7)) +
            vec[0][1] * (-x7 * vec[2][2] + x6 * (z2 - z7) + x2 * vec[1][2]))) / (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2])) /
            sqrt((vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2])));
        return;
    }

    // 7 6 4 3
    vec[0][0] = x6 - x7; vec[0][1] = y6 - y7; vec[0][2] = z6 - z7;
    vec[1][0] = x4 - x7; vec[1][1] = y4 - y7; vec[1][2] = z4 - z7;
    vec[2][0] = x3 - x7; vec[2][1] = y3 - y7; vec[2][2] = z3 - z7;
    vol = vec[0][0] * (vec[2][1] * vec[1][2] - vec[1][1] * vec[2][2]) +
        vec[0][1] * (vec[2][2] * vec[1][0] - vec[1][2] * vec[2][0]) +
        vec[0][2] * (vec[2][0] * vec[1][1] - vec[1][0] * vec[2][1]);
    len3 = sqrt(DOT(vec[0], vec[0]) * DOT(vec[1], vec[1]) * DOT(vec[2], vec[2]));
    if (vol >= 0) {
        eFSJ = false;
        dx3 = -sz2 * (y7 * (-z4 + z6) + y6 * vec[1][2] - y4 * vec[0][2]);
        dy3 = -sz2 * (x7 * (z4 - z6) + x4 * vec[0][2] - x6 * vec[1][2]);
        dz3 = -sz2 * (x7 * (-y4 + y6) + x6 * vec[1][1] - x4 * vec[0][1]);
        dx4 = -sz2 * (y7 * (z3 - z6) + y3 * vec[0][2] - y6 * vec[2][2]);
        dy4 = -sz2 * (x7 * (-z3 + z6) + x6 * vec[2][2] - x3 * vec[0][2]);
        dz4 = -sz2 * (x7 * (y3 - y6) + x3 * vec[0][1] - x6 * vec[2][1]);
        dx6 = -sz2 * (y7 * (-z3 + z4) + y4 * vec[2][2] - y3 * vec[1][2]);
        dy6 = -sz2 * (x7 * (z3 - z4) + x3 * vec[1][2] - x4 * vec[2][2]);
        dz6 = -sz2 * (x7 * (-y3 + y4) + x4 * vec[2][1] - x3 * vec[1][1]);
        dx7 = -sz2 * (y6 * (z3 - z4) + y3 * (z4 - z6) + y4 * (-z3 + z6));
        dy7 = -sz2 * (x6 * (-z3 + z4) + x4 * (z3 - z6) + x3 * (-z4 + z6));
        dz7 = -sz2 * (x6 * (y3 - y4) + x3 * (y4 - y6) + x4 * (-y3 + y6));
        return;
    }
    else if (vol > len3 * sJThres) {
        eFSJ = false;
        dx3 = -sz * ((y7 * (-z4 + z6) + y6 * vec[1][2] - y4 * vec[0][2] - (vec[2][0] * ((x7 * (-y3 + y4) + x4 * vec[2][1] - x3 * vec[1][1]) * vec[0][2] + vec[0][1] * (x7 * (z3 - z4) + x3 * vec[1][2] - x4 * vec[2][2]) +
            vec[0][0] * (y7 * (-z3 + z4) + y4 * vec[2][2] - y3 * vec[1][2]))) / (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])) /
            sqrt((vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2])));
        dy3 = -sz * ((x7 * (z4 - z6) + x4 * vec[0][2] - x6 * vec[1][2] - (vec[2][1] * ((x7 * (-y3 + y4) + x4 * vec[2][1] - x3 * vec[1][1]) * vec[0][2] + vec[0][1] * (x7 * (z3 - z4) + x3 * vec[1][2] - x4 * vec[2][2]) +
            vec[0][0] * (y7 * (-z3 + z4) + y4 * vec[2][2] - y3 * vec[1][2]))) / (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])) /
            sqrt((vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2])));
        dz3 = -sz * ((x7 * (-y4 + y6) + x6 * vec[1][1] - x4 * vec[0][1] - (vec[2][2] * ((x7 * (-y3 + y4) + x4 * vec[2][1] - x3 * vec[1][1]) * vec[0][2] + vec[0][1] * (x7 * (z3 - z4) + x3 * vec[1][2] - x4 * vec[2][2]) +
            vec[0][0] * (y7 * (-z3 + z4) + y4 * vec[2][2] - y3 * vec[1][2]))) / (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2])) /
            sqrt((vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2])));
        dx4 = -sz * ((y7 * (z3 - z6) + y3 * vec[0][2] - y6 * vec[2][2] - (vec[1][0] * ((x7 * (-y3 + y4) + x4 * vec[2][1] - x3 * vec[1][1]) * vec[0][2] + vec[0][1] * (x7 * (z3 - z4) + x3 * vec[1][2] - x4 * vec[2][2]) +
            vec[0][0] * (y7 * (-z3 + z4) + y4 * vec[2][2] - y3 * vec[1][2]))) / (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2])) /
            sqrt((vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2])));
        dy4 = -sz * ((x7 * (-z3 + z6) + x6 * vec[2][2] - x3 * vec[0][2] - (vec[1][1] * ((x7 * (-y3 + y4) + x4 * vec[2][1] - x3 * vec[1][1]) * vec[0][2] + vec[0][1] * (x7 * (z3 - z4) + x3 * vec[1][2] - x4 * vec[2][2]) +
            vec[0][0] * (y7 * (-z3 + z4) + y4 * vec[2][2] - y3 * vec[1][2]))) / (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2])) /
            sqrt((vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2])));
        dz4 = -sz * ((x7 * (y3 - y6) + x3 * vec[0][1] - x6 * vec[2][1] - (vec[1][2] * ((x7 * (-y3 + y4) + x4 * vec[2][1] - x3 * vec[1][1]) * vec[0][2] + vec[0][1] * (x7 * (z3 - z4) + x3 * vec[1][2] - x4 * vec[2][2]) +
            vec[0][0] * (y7 * (-z3 + z4) + y4 * vec[2][2] - y3 * vec[1][2]))) / (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2])) /
            sqrt((vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2])));
        dx6 = -sz * ((y7 * (-z3 + z4) + y4 * vec[2][2] - y3 * vec[1][2] - (vec[0][0] * ((x7 * (-y3 + y4) + x4 * vec[2][1] - x3 * vec[1][1]) * vec[0][2] + vec[0][1] * (x7 * (z3 - z4) + x3 * vec[1][2] - x4 * vec[2][2]) +
            vec[0][0] * (y7 * (-z3 + z4) + y4 * vec[2][2] - y3 * vec[1][2]))) / (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2])) /
            sqrt((vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2])));
        dy6 = -sz * ((x7 * (z3 - z4) + x3 * vec[1][2] - x4 * vec[2][2] - (vec[0][1] * ((x7 * (-y3 + y4) + x4 * vec[2][1] - x3 * vec[1][1]) * vec[0][2] + vec[0][1] * (x7 * (z3 - z4) + x3 * vec[1][2] - x4 * vec[2][2]) +
            vec[0][0] * (y7 * (-z3 + z4) + y4 * vec[2][2] - y3 * vec[1][2]))) / (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2])) /
            sqrt((vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2])));
        dz6 = -sz * ((x7 * (-y3 + y4) + x4 * vec[2][1] - x3 * vec[1][1] - (vec[0][2] * ((x7 * (-y3 + y4) + x4 * vec[2][1] - x3 * vec[1][1]) * vec[0][2] + vec[0][1] * (x7 * (z3 - z4) + x3 * vec[1][2] - x4 * vec[2][2]) +
            vec[0][0] * (y7 * (-z3 + z4) + y4 * vec[2][2] - y3 * vec[1][2]))) / (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2])) /
            sqrt((vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2])));
        dx7 = -sz * ((y6 * (z3 - z4) + y3 * (z4 - z6) + y4 * (-z3 + z6)) /
            sqrt((vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]))
            - ((-vec[0][0] * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) -
                vec[1][0] * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) -
                vec[2][0] * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2])) *
                ((x7 * (-y3 + y4) + x4 * vec[2][1] - x3 * vec[1][1]) * vec[0][2] + vec[0][1] * (x7 * (z3 - z4) + x3 * vec[1][2] - x4 * vec[2][2]) + vec[0][0] * (y7 * (-z3 + z4) + y4 * vec[2][2] - y3 * vec[1][2]))) /
            pow((vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]),
                1.5));
        dy7 = -sz * ((x6 * (-z3 + z4) + x4 * (z3 - z6) + x3 * (-z4 + z6)) /
            sqrt((vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]))
            - ((-vec[0][1] * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) -
                vec[1][1] * (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) -
                vec[2][1] * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2])) *
                ((x7 * (-y3 + y4) + x4 * vec[2][1] - x3 * vec[1][1]) * vec[0][2] + vec[0][1] * (x7 * (z3 - z4) + x3 * vec[1][2] - x4 * vec[2][2]) + vec[0][0] * (y7 * (-z3 + z4) + y4 * vec[2][2] - y3 * vec[1][2]))) /
            pow((vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]),
                1.5));
        dz7 = -sz * ((x6 * (y3 - y4) + x3 * (y4 - y6) + x4 * (-y3 + y6)) /
            sqrt((vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]))
            - (((vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) - vec[2][2] -
                (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]) * vec[1][2] -
                (vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * vec[0][2]) *
                ((x7 * (-y3 + y4) + x4 * vec[2][1] - x3 * vec[1][1]) * vec[0][2] + vec[0][1] * (x7 * (z3 - z4) + x3 * vec[1][2] - x4 * vec[2][2]) + vec[0][0] * (y7 * (-z3 + z4) + y4 * vec[2][2] - y3 * vec[1][2]))) /
            pow((vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]) * (vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]) * (vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]),
                1.5));
        return;
    }
}

bool gradient(const std::vector<std::vector<int>>& tri, const std::vector<std::vector<double>>& triX, const std::vector<std::vector<std::vector<double>>>& triEdgeX, const int& triENum, const std::vector<std::vector<double>>& edgeFtrX, const std::vector<std::vector<int>>& hexPtType, std::vector<double>& x, const std::vector<std::vector<int>>& hex, const int& eNum, const int& pNum3, const std::vector<int>& surf, const std::vector<int>& surfD3, const int& sPNum, const std::vector<int>& sPNumM3, std::vector<double>& surfX, std::vector<double>& totGrad, std::vector<int>& chkElem, const double& sJThres, bool& fitting) {
    // initialize vars
    int i;
    bool eFSJ = true;
    double maxGlobDist = 0; int aaaa;

    // update chkElem
#pragma omp parallel for
    for (i = 0; i < eNum; ++i)
        if (totGrad[hex[i][0]] != 0 || totGrad[hex[i][0] + 1] != 0 || totGrad[hex[i][0] + 2] != 0 ||
            totGrad[hex[i][1]] != 0 || totGrad[hex[i][1] + 1] != 0 || totGrad[hex[i][1] + 2] != 0 ||
            totGrad[hex[i][2]] != 0 || totGrad[hex[i][2] + 1] != 0 || totGrad[hex[i][2] + 2] != 0 ||
            totGrad[hex[i][3]] != 0 || totGrad[hex[i][3] + 1] != 0 || totGrad[hex[i][3] + 2] != 0 ||
            totGrad[hex[i][4]] != 0 || totGrad[hex[i][4] + 1] != 0 || totGrad[hex[i][4] + 2] != 0 ||
            totGrad[hex[i][5]] != 0 || totGrad[hex[i][5] + 1] != 0 || totGrad[hex[i][5] + 2] != 0 ||
            totGrad[hex[i][6]] != 0 || totGrad[hex[i][6] + 1] != 0 || totGrad[hex[i][6] + 2] != 0 ||
            totGrad[hex[i][7]] != 0 || totGrad[hex[i][7] + 1] != 0 || totGrad[hex[i][7] + 2] != 0)
            chkElem[i] = 1;
        else chkElem[i] = 0;
    std::fill(totGrad.begin(), totGrad.end(), 0.0);

    // calc mesh quality grad
#pragma omp parallel for
    for (i = 0; i < eNum; ++i) {
        if (chkElem[i]) {
            sJGrad(x[hex[i][0]], x[hex[i][0] + 1], x[hex[i][0] + 2],
                x[hex[i][1]], x[hex[i][1] + 1], x[hex[i][1] + 2],
                x[hex[i][2]], x[hex[i][2] + 1], x[hex[i][2] + 2],
                x[hex[i][3]], x[hex[i][3] + 1], x[hex[i][3] + 2],
                x[hex[i][4]], x[hex[i][4] + 1], x[hex[i][4] + 2],
                x[hex[i][5]], x[hex[i][5] + 1], x[hex[i][5] + 2],
                x[hex[i][6]], x[hex[i][6] + 1], x[hex[i][6] + 2],
                x[hex[i][7]], x[hex[i][7] + 1], x[hex[i][7] + 2],
                totGrad[hex[i][0]], totGrad[hex[i][0] + 1], totGrad[hex[i][0] + 2],
                totGrad[hex[i][1]], totGrad[hex[i][1] + 1], totGrad[hex[i][1] + 2],
                totGrad[hex[i][2]], totGrad[hex[i][2] + 1], totGrad[hex[i][2] + 2],
                totGrad[hex[i][3]], totGrad[hex[i][3] + 1], totGrad[hex[i][3] + 2],
                totGrad[hex[i][4]], totGrad[hex[i][4] + 1], totGrad[hex[i][4] + 2],
                totGrad[hex[i][5]], totGrad[hex[i][5] + 1], totGrad[hex[i][5] + 2],
                totGrad[hex[i][6]], totGrad[hex[i][6] + 1], totGrad[hex[i][6] + 2],
                totGrad[hex[i][7]], totGrad[hex[i][7] + 1], totGrad[hex[i][7] + 2],
                sJThres, eFSJ);
        }
    }

    // calc target pts
    if (eFSJ) {
        // calculate target pts
#pragma omp parallel for
        for (i = 0; i < sPNum; ++i)
            if (x[surf[i]] != surfX[sPNumM3[i]] || x[surf[i] + 1] != surfX[sPNumM3[i] + 1] || x[surf[i] + 2] != surfX[sPNumM3[i] + 2]) {
                if (hexPtType[surfD3[i]][0] == -2) {
                    surfX[sPNumM3[i]] = triX[hexPtType[surfD3[i]][1]][0];
                    surfX[sPNumM3[i] + 1] = triX[hexPtType[surfD3[i]][1]][1];
                    surfX[sPNumM3[i] + 2] = triX[hexPtType[surfD3[i]][1]][2];
                }
                else if (hexPtType[surfD3[i]][0] == -1) {
                    double dist, minDist = DBL_MAX, q[3];
                    for (int j = 1; j < hexPtType[surfD3[i]].size(); ++j) {
                        ptLnDist({ x[surf[i]], x[surf[i] + 1], x[surf[i] + 2] },
                            { edgeFtrX[hexPtType[surfD3[i]][j]][0], edgeFtrX[hexPtType[surfD3[i]][j]][1], edgeFtrX[hexPtType[surfD3[i]][j]][2] },
                            { edgeFtrX[hexPtType[surfD3[i]][j]][3], edgeFtrX[hexPtType[surfD3[i]][j]][4], edgeFtrX[hexPtType[surfD3[i]][j]][5] },
                            { edgeFtrX[hexPtType[surfD3[i]][j]][6], edgeFtrX[hexPtType[surfD3[i]][j]][7], edgeFtrX[hexPtType[surfD3[i]][j]][8], edgeFtrX[hexPtType[surfD3[i]][j]][9] },
                            q, dist);
                        if (dist < minDist) {
                            minDist = dist;
                            surfX[sPNumM3[i]] = q[0];
                            surfX[sPNumM3[i] + 1] = q[1];
                            surfX[sPNumM3[i] + 2] = q[2];
                        }
                    }
                }
                else {
                    double dist, minDist = DBL_MAX, q[3];
                    for (int j = 0; j < triENum; ++j) {
                        ptTriDist({ x[surf[i]], x[surf[i] + 1], x[surf[i] + 2] },
                            { triX[tri[j][0]][0], triX[tri[j][0]][1], triX[tri[j][0]][2] },
                            { triX[tri[j][1]][0], triX[tri[j][1]][1], triX[tri[j][1]][2] },
                            { triX[tri[j][2]][0], triX[tri[j][2]][1], triX[tri[j][2]][2] },
                            { triEdgeX[j][0][0], triEdgeX[j][0][1], triEdgeX[j][0][2] },
                            { triEdgeX[j][1][0], triEdgeX[j][1][1], triEdgeX[j][1][2] },
                            { triEdgeX[j][2][0], triEdgeX[j][2][1], triEdgeX[j][2][2] },
                            q, dist);
                        if (dist < minDist) {
                            minDist = dist;
                            surfX[sPNumM3[i]] = q[0];
                            surfX[sPNumM3[i] + 1] = q[1];
                            surfX[sPNumM3[i] + 2] = q[2];
                        }
                    }
                }
            }

        // calc fitting grad
        fitting = true;
#pragma omp parallel for
        for (i = 0; i < sPNum; ++i) {
            double dist = sqrt((x[surf[i]] - surfX[sPNumM3[i]]) * (x[surf[i]] - surfX[sPNumM3[i]]) + (x[surf[i] + 1] - surfX[sPNumM3[i] + 1]) * (x[surf[i] + 1] - surfX[sPNumM3[i] + 1]) + (x[surf[i] + 2] - surfX[sPNumM3[i] + 2]) * (x[surf[i] + 2] - surfX[sPNumM3[i] + 2]));
            if (maxGlobDist < dist) {
                maxGlobDist = dist;
                aaaa = surfD3[i];
            }
            if (dist > 1e-8) {
                fitting = false;
                totGrad[surf[i]] = 1.0;
                x[surf[i]] = 0.5 * x[surf[i]] + 0.5 * surfX[sPNumM3[i]];
                x[surf[i] + 1] = 0.5 * x[surf[i] + 1] + 0.5 * surfX[sPNumM3[i] + 1];
                x[surf[i] + 2] = 0.5 * x[surf[i] + 2] + 0.5 * surfX[sPNumM3[i] + 2];
            }
        }
        std::cout << maxGlobDist << " " << aaaa << " " << fitting << std::endl;
        return true;
    }
    else {
#pragma omp parallel for
        for (i = 0; i < pNum3; i+=3) {
            double norm = totGrad[i] * totGrad[i] + totGrad[i + 1] * totGrad[i + 1] + totGrad[i + 2] * totGrad[i + 2];
            if (norm > 0) {
                //norm = 1.0 / sqrt(norm);
                norm = std::min(1.0, 1.0 / sqrt(norm));
                x[i] -= totGrad[i] * norm * 1e-3;
                x[i + 1] -= totGrad[i + 1] * norm * 1e-3;
                x[i + 2] -= totGrad[i + 2] * norm * 1e-3;
            }
        }
        return false;
    }
}