#include "geometry.h"

void ptTriDist(const double(&p)[3], const double(&a)[3], const double(&b)[3], const double(&c)[3],
	const double(&ab)[3], const double(&ac)[3], const double(&bc)[3],
	double(&q)[3], double& minDist) {
	const double ap[3] = { p[0] - a[0], p[1] - a[1], p[2] - a[2] };

	// corner1
	const double d1 = DOT(ab, ap), d2 = DOT(ac, ap);
	if (d1 <= 0 && d2 <= 0) {
		q[0] = a[0];
		q[1] = a[1];
		q[2] = a[2];
		minDist = (p[0] - q[0]) * (p[0] - q[0]) + (p[1] - q[1]) * (p[1] - q[1]) + (p[2] - q[2]) * (p[2] - q[2]);
		return;
	}

	// corner2
	const double bp[3] = { p[0] - b[0], p[1] - b[1], p[2] - b[2] },
		d3 = DOT(ab, bp), d4 = DOT(ac, bp);
	if (d3 >= 0 && d4 <= d3) {
		q[0] = b[0];
		q[1] = b[1];
		q[2] = b[2];
		minDist = (p[0] - q[0]) * (p[0] - q[0]) + (p[1] - q[1]) * (p[1] - q[1]) + (p[2] - q[2]) * (p[2] - q[2]);
		return;
	}

	// corner3
	const double cp[3] = { p[0] - c[0], p[1] - c[1], p[2] - c[2] },
		d5 = DOT(ab, cp), d6 = DOT(ac, cp);
	if (d6 >= 0 && d5 <= d6) {
		q[0] = c[0];
		q[1] = c[1];
		q[2] = c[2];
		minDist = (p[0] - q[0]) * (p[0] - q[0]) + (p[1] - q[1]) * (p[1] - q[1]) + (p[2] - q[2]) * (p[2] - q[2]);
		return;
	}

	// edge1
	const double vc = d1 * d4 - d3 * d2;
	if (vc <= 0 && d1 >= 0 && d3 <= 0) {
		const double v = d1 / (d1 - d3);
		q[0] = a[0] + v * ab[0];
		q[1] = a[1] + v * ab[1];
		q[2] = a[2] + v * ab[2];
		minDist = (p[0] - q[0]) * (p[0] - q[0]) + (p[1] - q[1]) * (p[1] - q[1]) + (p[2] - q[2]) * (p[2] - q[2]);
		return;
	}

	// edge2
	const double vb = d5 * d2 - d1 * d6;
	if (vb <= 0 && d2 >= 0 && d6 <= 0) {
		const double v = d2 / (d2 - d6);
		q[0] = a[0] + v * ac[0];
		q[1] = a[1] + v * ac[1];
		q[2] = a[2] + v * ac[2];
		minDist = (p[0] - q[0]) * (p[0] - q[0]) + (p[1] - q[1]) * (p[1] - q[1]) + (p[2] - q[2]) * (p[2] - q[2]);
		return;
	}

	// edge3
	const double va = d3 * d6 - d5 * d4;
	if (va <= 0 && d4 >= d3 && d5 >= d6) {
		const double v = (d4 - d3) / ((d4 - d3) + (d5 - d6));
		q[0] = b[0] + v * bc[0];
		q[1] = b[1] + v * bc[1];
		q[2] = b[2] + v * bc[2];
		minDist = (p[0] - q[0]) * (p[0] - q[0]) + (p[1] - q[1]) * (p[1] - q[1]) + (p[2] - q[2]) * (p[2] - q[2]);
		return;
	}

	// inside
	const double denom = 1.0 / (va + vb + vc);
	const double v = vb * denom;
	const double w = vc * denom;
	q[0] = a[0] + v * ab[0] + w * ac[0];
	q[1] = a[1] + v * ab[1] + w * ac[1];
	q[2] = a[2] + v * ab[2] + w * ac[2];
	minDist = (p[0] - q[0]) * (p[0] - q[0]) + (p[1] - q[1]) * (p[1] - q[1]) + (p[2] - q[2]) * (p[2] - q[2]);
}

void ptLnDist(const double(&p)[3], const double(&a)[3], const double(&b)[3],
	const double(&ab)[4],
	double(&q)[3], double& minDist) {
	const double ap[3] = { p[0] - a[0], p[1] - a[1], p[2] - a[2] };

	//left
	const double d = DOT(ab, ap);
	if (d <= 0) {
		q[0] = a[0];
		q[1] = a[1];
		q[2] = a[2];
		minDist = (p[0] - q[0]) * (p[0] - q[0]) + (p[1] - q[1]) * (p[1] - q[1]) + (p[2] - q[2]) * (p[2] - q[2]);
		return;
	}

	// right
	const double bp[3] = { p[0] - b[0], p[1] - b[1], p[2] - b[2] };
	if (DOT(ab, bp) >= 0) {
		q[0] = b[0];
		q[1] = b[1];
		q[2] = b[2];
		minDist = (p[0] - q[0]) * (p[0] - q[0]) + (p[1] - q[1]) * (p[1] - q[1]) + (p[2] - q[2]) * (p[2] - q[2]);
		return;
	}

	// inside
	const double v = d / ab[3];
	q[0] = a[0] + v * ab[0];
	q[1] = a[1] + v * ab[1];
	q[2] = a[2] + v * ab[2];
	minDist = (p[0] - q[0]) * (p[0] - q[0]) + (p[1] - q[1]) * (p[1] - q[1]) + (p[2] - q[2]) * (p[2] - q[2]);
}

size_t hashFace(const int& idx0, const int& idx1, const int& idx2, const int& idx3) {
	int pts[4] = {idx0, idx1, idx2, idx3};
	if (pts[0] > pts[1])
		std::swap(pts[0], pts[1]);
	if (pts[2] > pts[3])
		std::swap(pts[2], pts[3]);
	if (pts[0] > pts[2])
		std::swap(pts[0], pts[2]);
	if (pts[1] > pts[3])
		std::swap(pts[1], pts[3]);
	if (pts[1] > pts[2])
		std::swap(pts[1], pts[2]);
	return std::hash<std::string>{}(std::to_string(pts[0]) + " " + std::to_string(pts[1]) + " " + std::to_string(pts[2]) + " " + std::to_string(pts[3]));
}

size_t hashEdge(const int& idx0, const int& idx1) {
	if (idx0 > idx1) return std::hash<std::string>{}(std::to_string(idx1) + " " + std::to_string(idx0));
	else return std::hash<std::string>{}(std::to_string(idx0) + " " + std::to_string(idx1));
}