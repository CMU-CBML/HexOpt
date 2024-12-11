#pragma once

#include <iostream>
#include <omp.h>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include "constants.h"
#include "geometry.h"

// Laplacian smoothing
void lapSmth(std::vector<double>& x, const std::vector<std::vector<int>>& adjPts, const int& pNum, const std::vector<int>& pNumM3);

// calculate (scaled) Jacobian gradient
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
    const double& sJThres, bool& eFSJ);

// get the gradient of all the points
bool gradient(const std::vector<std::vector<int>>& tri, const std::vector<std::vector<double>>& triX, const std::vector<std::vector<std::vector<double>>>& triEdgeX, const int& triENum, const std::vector<std::vector<double>>& edgeFtrX, const std::vector<std::vector<int>>& hexPtType, std::vector<double>& x, const std::vector<std::vector<int>>& hex, const int& eNum, const int& pNum3, const std::vector<int>& surf, const std::vector<int>& surfD3, const int& sPNum, const std::vector<int>& sPNumM3, std::vector<double>& surfX, std::vector<double>& totGrad, std::vector<int>& chkElem, const double& sJThres, bool& fitting);