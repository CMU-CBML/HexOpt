#pragma once

#include <string>
#include "constants.h"

//Get the closest point on a triangle to a point and the distance
void ptTriDist(const double(&p)[3], const double(&a)[3], const double(&b)[3], const double(&c)[3],
	const double(&ab)[3], const double(&ac)[3], const double(&bc)[3],
	double(&q)[3], double& minDist);

//Get the closest point on a line to a point and the distance
void ptLnDist(const double(&p)[3], const double(&a)[3], const double(&b)[3],
	const double(&ab)[4],
	double(&q)[3], double& minDist);

//Save face given the indexes of four points
size_t hashFace(const int& idx0, const int& idx1, const int& idx2, const int& idx3);

//Save edge given the indexes of two points
size_t hashEdge(const int& idx0, const int& idx1);