#ifndef MESH_OPTIMIZATION_H
#define MESH_OPTIMIZATION_H

#include <iostream>
#include <random>
#include <unordered_set>
#include "Mesh_Shortest_Path.h"
#include "Mesh_Input_Output.h"

namespace MeshOptimization {
	size_t MeshQualityBelowThreshold(
		const MeshTypes::Vec3& p0,
		const MeshTypes::Vec3& p1,
		const MeshTypes::Vec3& p2,
		const MeshTypes::Vec3& p3,
		const MeshTypes::Vec3& p4,
		const MeshTypes::Vec3& p5,
		const MeshTypes::Vec3& p6,
		const MeshTypes::Vec3& p7) {
		// J at the center
		double x0 = p1.x + p2.x + p5.x + p6.x - p0.x - p3.x - p4.x - p7.x;
		double y0 = p1.y + p2.y + p5.y + p6.y - p0.y - p3.y - p4.y - p7.y;
		double z0 = p1.z + p2.z + p5.z + p6.z - p0.z - p3.z - p4.z - p7.z;

		double x1 = p2.x + p3.x + p6.x + p7.x - p0.x - p1.x - p4.x - p5.x;
		double y1 = p2.y + p3.y + p6.y + p7.y - p0.y - p1.y - p4.y - p5.y;
		double z1 = p2.z + p3.z + p6.z + p7.z - p0.z - p1.z - p4.z - p5.z;

		double x2 = p4.x + p5.x + p6.x + p7.x - p0.x - p1.x - p2.x - p3.x;
		double y2 = p4.y + p5.y + p6.y + p7.y - p0.y - p1.y - p2.y - p3.y;
		double z2 = p4.z + p5.z + p6.z + p7.z - p0.z - p1.z - p2.z - p3.z;

		double volume = x0 * (y1 * z2 - y2 * z1) + y0 * (z1 * x2 - z2 * x1) + z0 * (x1 * y2 - x2 * y1);

		if (volume <= 0) {
			return 0;
		}

		// J(0,0,0)
		x0 = p1.x - p0.x;
		y0 = p1.y - p0.y;
		z0 = p1.z - p0.z;

		x1 = p3.x - p0.x;
		y1 = p3.y - p0.y;
		z1 = p3.z - p0.z;

		x2 = p4.x - p0.x;
		y2 = p4.y - p0.y;
		z2 = p4.z - p0.z;

		volume = x0 * (y1 * z2 - y2 * z1) + y0 * (z1 * x2 - z2 * x1) + z0 * (x1 * y2 - x2 * y1);

		if (volume <= 0) {
			return 1;
		}

		// J(1,0,0)
		x0 = p2.x - p1.x;
		y0 = p2.y - p1.y;
		z0 = p2.z - p1.z;

		x1 = p0.x - p1.x;
		y1 = p0.y - p1.y;
		z1 = p0.z - p1.z;

		x2 = p5.x - p1.x;
		y2 = p5.y - p1.y;
		z2 = p5.z - p1.z;

		volume = x0 * (y1 * z2 - y2 * z1) + y0 * (z1 * x2 - z2 * x1) + z0 * (x1 * y2 - x2 * y1);

		if (volume <= 0) {
			return 2;
		}

		// J(1,1,0)
		x0 = p3.x - p2.x;
		y0 = p3.y - p2.y;
		z0 = p3.z - p2.z;

		x1 = p1.x - p2.x;
		y1 = p1.y - p2.y;
		z1 = p1.z - p2.z;

		x2 = p6.x - p2.x;
		y2 = p6.y - p2.y;
		z2 = p6.z - p2.z;

		volume = x0 * (y1 * z2 - y2 * z1) + y0 * (z1 * x2 - z2 * x1) + z0 * (x1 * y2 - x2 * y1);

		if (volume <= 0) {
			return 3;
		}

		// J(0,1,0)
		x0 = p0.x - p3.x;
		y0 = p0.y - p3.y;
		z0 = p0.z - p3.z;

		x1 = p2.x - p3.x;
		y1 = p2.y - p3.y;
		z1 = p2.z - p3.z;

		x2 = p7.x - p3.x;
		y2 = p7.y - p3.y;
		z2 = p7.z - p3.z;

		volume = x0 * (y1 * z2 - y2 * z1) + y0 * (z1 * x2 - z2 * x1) + z0 * (x1 * y2 - x2 * y1);

		if (volume <= 0) {
			return 4;
		}

		// J(0,0,1)
		x0 = p7.x - p4.x;
		y0 = p7.y - p4.y;
		z0 = p7.z - p4.z;

		x1 = p5.x - p4.x;
		y1 = p5.y - p4.y;
		z1 = p5.z - p4.z;

		x2 = p0.x - p4.x;
		y2 = p0.y - p4.y;
		z2 = p0.z - p4.z;

		volume = x0 * (y1 * z2 - y2 * z1) + y0 * (z1 * x2 - z2 * x1) + z0 * (x1 * y2 - x2 * y1);

		if (volume <= 0) {
			return 5;
		}

		// J(1,0,1)
		x0 = p4.x - p5.x;
		y0 = p4.y - p5.y;
		z0 = p4.z - p5.z;

		x1 = p6.x - p5.x;
		y1 = p6.y - p5.y;
		z1 = p6.z - p5.z;

		x2 = p1.x - p5.x;
		y2 = p1.y - p5.y;
		z2 = p1.z - p5.z;

		volume = x0 * (y1 * z2 - y2 * z1) + y0 * (z1 * x2 - z2 * x1) + z0 * (x1 * y2 - x2 * y1);

		if (volume <= 0) {
			return 6;
		}

		// J(1,1,1)
		x0 = p5.x - p6.x;
		y0 = p5.y - p6.y;
		z0 = p5.z - p6.z;

		x1 = p7.x - p6.x;
		y1 = p7.y - p6.y;
		z1 = p7.z - p6.z;

		x2 = p2.x - p6.x;
		y2 = p2.y - p6.y;
		z2 = p2.z - p6.z;

		volume = x0 * (y1 * z2 - y2 * z1) + y0 * (z1 * x2 - z2 * x1) + z0 * (x1 * y2 - x2 * y1);

		if (volume <= 0) {
			return 7;
		}

		// J(0,1,1)
		x0 = p6.x - p7.x;
		y0 = p6.y - p7.y;
		z0 = p6.z - p7.z;

		x1 = p4.x - p7.x;
		y1 = p4.y - p7.y;
		z1 = p4.z - p7.z;

		x2 = p3.x - p7.x;
		y2 = p3.y - p7.y;
		z2 = p3.z - p7.z;

		volume = x0 * (y1 * z2 - y2 * z1) + y0 * (z1 * x2 - z2 * x1) + z0 * (x1 * y2 - x2 * y1);

		if (volume <= 0) {
			return 8;
		}

		return 9;
	}

	std::pair<size_t, double> MinMeshQuality(
		const MeshTypes::Vec3& p0,
		const MeshTypes::Vec3& p1,
		const MeshTypes::Vec3& p2,
		const MeshTypes::Vec3& p3,
		const MeshTypes::Vec3& p4,
		const MeshTypes::Vec3& p5,
		const MeshTypes::Vec3& p6,
		const MeshTypes::Vec3& p7) {
		// J at the center
		double x0 = p1.x + p2.x + p5.x + p6.x - p0.x - p3.x - p4.x - p7.x;
		double y0 = p1.y + p2.y + p5.y + p6.y - p0.y - p3.y - p4.y - p7.y;
		double z0 = p1.z + p2.z + p5.z + p6.z - p0.z - p3.z - p4.z - p7.z;

		double x1 = p2.x + p3.x + p6.x + p7.x - p0.x - p1.x - p4.x - p5.x;
		double y1 = p2.y + p3.y + p6.y + p7.y - p0.y - p1.y - p4.y - p5.y;
		double z1 = p2.z + p3.z + p6.z + p7.z - p0.z - p1.z - p4.z - p5.z;

		double x2 = p4.x + p5.x + p6.x + p7.x - p0.x - p1.x - p2.x - p3.x;
		double y2 = p4.y + p5.y + p6.y + p7.y - p0.y - p1.y - p2.y - p3.y;
		double z2 = p4.z + p5.z + p6.z + p7.z - p0.z - p1.z - p2.z - p3.z;

		double min_volume = (x0 * (y1 * z2 - y2 * z1) + y0 * (z1 * x2 - z2 * x1) + z0 * (x1 * y2 - x2 * y1)) / 64.0;

		size_t min_volume_index = 0;

		// J(0,0,0)
		x0 = p1.x - p0.x;
		y0 = p1.y - p0.y;
		z0 = p1.z - p0.z;

		x1 = p3.x - p0.x;
		y1 = p3.y - p0.y;
		z1 = p3.z - p0.z;

		x2 = p4.x - p0.x;
		y2 = p4.y - p0.y;
		z2 = p4.z - p0.z;

		double volume = x0 * (y1 * z2 - y2 * z1) + y0 * (z1 * x2 - z2 * x1) + z0 * (x1 * y2 - x2 * y1);

		if (volume < min_volume) {
			min_volume = volume;
			min_volume_index = 1;
		}

		// J(1,0,0)
		x0 = p2.x - p1.x;
		y0 = p2.y - p1.y;
		z0 = p2.z - p1.z;

		x1 = p0.x - p1.x;
		y1 = p0.y - p1.y;
		z1 = p0.z - p1.z;

		x2 = p5.x - p1.x;
		y2 = p5.y - p1.y;
		z2 = p5.z - p1.z;

		volume = x0 * (y1 * z2 - y2 * z1) + y0 * (z1 * x2 - z2 * x1) + z0 * (x1 * y2 - x2 * y1);

		if (volume < min_volume) {
			min_volume = volume;
			min_volume_index = 2;
		}

		// J(1,1,0)
		x0 = p3.x - p2.x;
		y0 = p3.y - p2.y;
		z0 = p3.z - p2.z;

		x1 = p1.x - p2.x;
		y1 = p1.y - p2.y;
		z1 = p1.z - p2.z;

		x2 = p6.x - p2.x;
		y2 = p6.y - p2.y;
		z2 = p6.z - p2.z;

		volume = x0 * (y1 * z2 - y2 * z1) + y0 * (z1 * x2 - z2 * x1) + z0 * (x1 * y2 - x2 * y1);

		if (volume < min_volume) {
			min_volume = volume;
			min_volume_index = 3;
		}

		// J(0,1,0)
		x0 = p0.x - p3.x;
		y0 = p0.y - p3.y;
		z0 = p0.z - p3.z;

		x1 = p2.x - p3.x;
		y1 = p2.y - p3.y;
		z1 = p2.z - p3.z;

		x2 = p7.x - p3.x;
		y2 = p7.y - p3.y;
		z2 = p7.z - p3.z;

		volume = x0 * (y1 * z2 - y2 * z1) + y0 * (z1 * x2 - z2 * x1) + z0 * (x1 * y2 - x2 * y1);

		if (volume < min_volume) {
			min_volume = volume;
			min_volume_index = 4;
		}

		// J(0,0,1)
		x0 = p7.x - p4.x;
		y0 = p7.y - p4.y;
		z0 = p7.z - p4.z;

		x1 = p5.x - p4.x;
		y1 = p5.y - p4.y;
		z1 = p5.z - p4.z;

		x2 = p0.x - p4.x;
		y2 = p0.y - p4.y;
		z2 = p0.z - p4.z;

		volume = x0 * (y1 * z2 - y2 * z1) + y0 * (z1 * x2 - z2 * x1) + z0 * (x1 * y2 - x2 * y1);

		if (volume < min_volume) {
			min_volume = volume;
			min_volume_index = 5;
		}

		// J(1,0,1)
		x0 = p4.x - p5.x;
		y0 = p4.y - p5.y;
		z0 = p4.z - p5.z;

		x1 = p6.x - p5.x;
		y1 = p6.y - p5.y;
		z1 = p6.z - p5.z;

		x2 = p1.x - p5.x;
		y2 = p1.y - p5.y;
		z2 = p1.z - p5.z;

		volume = x0 * (y1 * z2 - y2 * z1) + y0 * (z1 * x2 - z2 * x1) + z0 * (x1 * y2 - x2 * y1);

		if (volume < min_volume) {
			min_volume = volume;
			min_volume_index = 6;
		}

		// J(1,1,1)
		x0 = p5.x - p6.x;
		y0 = p5.y - p6.y;
		z0 = p5.z - p6.z;

		x1 = p7.x - p6.x;
		y1 = p7.y - p6.y;
		z1 = p7.z - p6.z;

		x2 = p2.x - p6.x;
		y2 = p2.y - p6.y;
		z2 = p2.z - p6.z;

		volume = x0 * (y1 * z2 - y2 * z1) + y0 * (z1 * x2 - z2 * x1) + z0 * (x1 * y2 - x2 * y1);

		if (volume < min_volume) {
			min_volume = volume;
			min_volume_index = 7;
		}

		// J(0,1,1)
		x0 = p6.x - p7.x;
		y0 = p6.y - p7.y;
		z0 = p6.z - p7.z;

		x1 = p4.x - p7.x;
		y1 = p4.y - p7.y;
		z1 = p4.z - p7.z;

		x2 = p3.x - p7.x;
		y2 = p3.y - p7.y;
		z2 = p3.z - p7.z;

		volume = x0 * (y1 * z2 - y2 * z1) + y0 * (z1 * x2 - z2 * x1) + z0 * (x1 * y2 - x2 * y1);

		if (volume < min_volume) {
			min_volume = volume;
			min_volume_index = 8;
		}

		return std::make_pair(min_volume_index, min_volume);
	}

	size_t RayTriangleIntersect(
		const MeshTypes::Vec3& origin,
		const MeshTypes::Vec3& ray,
		const MeshTypes::Vec3& v0,
		const MeshTypes::Vec3& v1,
		const MeshTypes::Vec3& v2) {
		const MeshTypes::Vec3 edge1 = v1 - v0;
		const MeshTypes::Vec3 edge2 = v2 - v0;
		const MeshTypes::Vec3 pvec = ray.cross(edge2);
		const double det = edge1.dot(pvec);

		if (det == 0) {
			return 0;
		}
		const double inv_det = 1.0 / det;

		const MeshTypes::Vec3 tvec = origin - v0;
		const double u = tvec.dot(pvec) * inv_det;
		const MeshTypes::Vec3 qvec = tvec.cross(edge1);
		const double v = ray.dot(qvec) * inv_det;
		const double t = (edge2.dot(qvec)) * inv_det;
		if (t >= 0 && u >= 0 && v >= 0 && u + v <= 1) {
			if (u > 0 && v > 0 && u + v < 1) {
				return 1;
			}
			return 2;
		}
		return 0;
	}

	MeshTypes::Vec3 VertexTriangleDistance(
		const MeshTypes::Vec3& p,
		const MeshTypes::Vec3& a,
		const MeshTypes::Vec3& b,
		const MeshTypes::Vec3& c) {
		const MeshTypes::Vec3 ap = p - a;
		const MeshTypes::Vec3 ab = b - a;
		const MeshTypes::Vec3 ac = c - a;
		const double d1 = ab.dot(ap);
		const double d2 = ac.dot(ap);
		MeshTypes::Vec3 pq;
		if (d1 <= 0 && d2 <= 0) {
			// q = a
			pq = a - p;
			return pq;
		}

		const MeshTypes::Vec3 bp = p - b;
		const double d3 = ab.dot(bp);
		const double d4 = ac.dot(bp);
		if (d3 >= 0 && d4 <= d3) {
			// q = b
			pq = b - p;
			return pq;
		}

		const MeshTypes::Vec3 cp = p - c;
		const double d5 = ab.dot(cp);
		const double d6 = ac.dot(cp);
		if (d6 >= 0 && d5 <= d6) {
			// q = c
			pq = c - p;
			return pq;
		}

		const double vc = d1 * d4 - d3 * d2;
		if (vc <= 0 && d1 >= 0 && d3 <= 0) {
			const double v = d1 / (d1 - d3);
			// q = a + ab * v
			pq = a + ab * v - p;
			return pq;
		}

		const double vb = d5 * d2 - d1 * d6;
		if (vb <= 0 && d2 >= 0 && d6 <= 0) {
			const double w = d2 / (d2 - d6);
			// q = a + ac * w
			pq = a + ac * w - p;
			return pq;
		}

		const double va = d3 * d6 - d5 * d4;
		if (va <= 0 && d4 >= d3 && d5 >= d6) {
			const double v = (d4 - d3) / ((d4 - d3) + (d5 - d6));
			// q = b + (c - b) * v
			pq = b + (c - b) * v - p;
			return pq;
		}

		const double denom = 1.0 / (va + vb + vc);
		const double v = vb * denom;
		const double w = vc * denom;
		// q = a + ab * v + ac * w
		pq = a + ab * v + ac * w - p;
		return pq;
	}

	class MeshOptimizer {
	public:
		MeshOptimizer(const MeshTypes::Mesh3D& triangle_mesh,
			const MeshTypes::SharpFeature& sharp_feature,
			const MeshTypes::Mesh3D& hexahedra_mesh)
			: triangle_mesh_(triangle_mesh),
			sharp_feature_(sharp_feature),
			hexahedra_mesh_(hexahedra_mesh) {
		}

		MeshTypes::Mesh3D Optimize(
			double learning_rate,
			double lambda_projection,
			const double max_distance_threshold
		) {
			const double internal_learning_rate = learning_rate;
			const double learning_rate_projection = learning_rate * lambda_projection;
			const size_t original_element_number = hexahedra_mesh_.faces.size();

			std::unordered_map<size_t, std::vector<size_t>> faces_on_surface_set;
			faces_on_surface_set.reserve(original_element_number);
			const std::vector<std::vector<size_t>> hex_face_vertices =
			{ { 0, 1, 2, 3, 4, 5, 6, 7 },
			  { 4, 5, 1, 0, 7, 6, 2, 3 },
			  { 4, 0, 3, 7, 5, 1, 2, 6 },
			  { 6, 2, 1, 5, 7, 3, 0, 4 },
			  { 3, 2, 6, 7, 0, 1, 5, 4 },
			  { 7, 6, 5, 4, 3, 2, 1, 0 } };
			for (size_t i = 0; i < original_element_number; ++i) {
				for (size_t j = 0; j < 6; ++j) {
					std::vector<size_t> face_array = {
						hexahedra_mesh_.faces[i][hex_face_vertices[j][0]],
						hexahedra_mesh_.faces[i][hex_face_vertices[j][1]],
						hexahedra_mesh_.faces[i][hex_face_vertices[j][2]],
						hexahedra_mesh_.faces[i][hex_face_vertices[j][3]] };
					const size_t face_array_hash = HASH_FACE(face_array);
					face_array.push_back(i);
					face_array.push_back(j);
					if (!faces_on_surface_set.count(face_array_hash)) {
						faces_on_surface_set[face_array_hash] = face_array;
					}
					else {
						faces_on_surface_set[face_array_hash].clear();
					}
				}
			}
			size_t original_vertex_number = hexahedra_mesh_.vertices.size();

			// handle sharp feature
			MeshTypes::Mesh3D surface_quad_mesh_;
			surface_quad_mesh_.vertices = hexahedra_mesh_.vertices;
			std::vector<bool> is_vertex_active_(original_vertex_number, false);
			std::vector<bool> is_vertex_on_surface_(original_vertex_number, false);
			for (auto i = faces_on_surface_set.begin(); i != faces_on_surface_set.end();) {
				const std::vector<size_t>& face_array = i->second;
				if (!face_array.empty()) {
					for (size_t j = 0; j < 4; ++j) {
						is_vertex_active_[face_array[j]] = true;
						is_vertex_on_surface_[face_array[j]] = true;
					}
					surface_quad_mesh_.faces.push_back(MeshTypes::VecNI({ face_array[0], face_array[1], face_array[2], face_array[3] }));
					++i;
				}
				else {
					i = faces_on_surface_set.erase(i);
				}
			}
			std::unordered_map<size_t, size_t> sharp_vertex_indices;
			std::vector<long long> where_should_end_vertex_go(original_vertex_number, -2);
			std::vector<std::pair<size_t, size_t>> where_should_middle_vertex_go(original_vertex_number);
			std::vector<size_t> vertex_candidate_distance(original_vertex_number, std::numeric_limits<double>::max());
			auto FindClosestVertex = [&](size_t triangle_vert) {
				double min_distance = std::numeric_limits<double>::max();
				size_t closest_vert = 0;
				for (size_t i = 0; i < original_vertex_number; ++i) {
					if (!is_vertex_on_surface_[i]) {
						continue;
					}
					const auto diff = surface_quad_mesh_.vertices[i] - triangle_mesh_.vertices[triangle_vert];
					const auto d2 = diff.dot(diff);
					if (d2 < min_distance) {
						min_distance = d2;
						closest_vert = i;
					}
				}
				is_vertex_active_[closest_vert] = false;
				sharp_vertex_indices[triangle_vert] = closest_vert;
				if (min_distance < vertex_candidate_distance[closest_vert]) {
					vertex_candidate_distance[closest_vert] = min_distance;
					where_should_end_vertex_go[closest_vert] = triangle_vert;
				}
				};
			for (const auto& [start, end] : sharp_feature_.edges) {
				if (!sharp_vertex_indices.count(start)) {
					FindClosestVertex(start);
				}
				if (!sharp_vertex_indices.count(end)) {
					FindClosestVertex(end);
				}
			}
			for (const auto& vert : sharp_feature_.vertices) {
				if (!sharp_vertex_indices.count(vert)) {
					FindClosestVertex(vert);
				}
			}

			MeshShortestPath::DijkstraSolver dijkstra(surface_quad_mesh_, is_vertex_active_);
			for (const auto& [start, end] : sharp_feature_.edges) {
				const auto& path_indices = dijkstra.ComputePath(sharp_vertex_indices[start],
					sharp_vertex_indices[end],
					true);
				if (path_indices.empty()) {
					throw std::runtime_error("Cannot compute valid sharp edges.");
				}
				for (size_t i = 1; i < path_indices.size() - 1; ++i) {
					where_should_end_vertex_go[path_indices[i]] = -1;
					where_should_middle_vertex_go[path_indices[i]] = std::make_pair(path_indices[i - 1], path_indices[i + 1]);
				}
			}
			std::vector<size_t> need_refine(original_element_number, 6);
			for (auto i = faces_on_surface_set.begin(); i != faces_on_surface_set.end();) {
				const auto& face_array = i->second;
				bool goto_next = false;
				for (size_t j = 0; j < 4; ++j) {
					if (!is_vertex_active_[face_array[j]] &&
						!is_vertex_active_[face_array[(j + 1) % 4]] &&
						!is_vertex_active_[face_array[(j + 2) % 4]]) {
						need_refine[face_array[4]] = face_array[5];
						i = faces_on_surface_set.erase(i);
						goto_next = true;
						break;
					}
				}
				if (!goto_next) {
					++i;
				}
			}
			
			// refine elements that need to be refined
			// operate directly on optimized_hexahedra_mesh_
			// variables to be updated:
			// optimized_hexahedra_mesh
			// original_vertex_number
			// faces_on_surface_set
			MeshTypes::Mesh3D optimized_hexahedra_mesh_ = hexahedra_mesh_;
			for (size_t i = 0; i < original_element_number; ++i) {
				if (need_refine[i] < 6) {
					MeshTypes::Vec3 face_center =
						0.25 * hexahedra_mesh_.vertices[hexahedra_mesh_.faces[i][hex_face_vertices[need_refine[i]][0]]] +
						0.25 * hexahedra_mesh_.vertices[hexahedra_mesh_.faces[i][hex_face_vertices[need_refine[i]][1]]] +
						0.25 * hexahedra_mesh_.vertices[hexahedra_mesh_.faces[i][hex_face_vertices[need_refine[i]][2]]] +
						0.25 * hexahedra_mesh_.vertices[hexahedra_mesh_.faces[i][hex_face_vertices[need_refine[i]][3]]];
					optimized_hexahedra_mesh_.vertices.push_back(
						0.5 * hexahedra_mesh_.vertices[hexahedra_mesh_.faces[i][hex_face_vertices[need_refine[i]][0]]] +
						0.5 * face_center);
					optimized_hexahedra_mesh_.vertices.push_back(
						0.5 * hexahedra_mesh_.vertices[hexahedra_mesh_.faces[i][hex_face_vertices[need_refine[i]][1]]] +
						0.5 * face_center);
					optimized_hexahedra_mesh_.vertices.push_back(
						0.5 * hexahedra_mesh_.vertices[hexahedra_mesh_.faces[i][hex_face_vertices[need_refine[i]][2]]] +
						0.5 * face_center);
					optimized_hexahedra_mesh_.vertices.push_back(
						0.5 * hexahedra_mesh_.vertices[hexahedra_mesh_.faces[i][hex_face_vertices[need_refine[i]][3]]] +
						0.5 * face_center);
					optimized_hexahedra_mesh_.vertices.push_back(
						0.5 * hexahedra_mesh_.vertices[hexahedra_mesh_.faces[i][hex_face_vertices[need_refine[i]][4]]] +
						0.5 * face_center);
					optimized_hexahedra_mesh_.vertices.push_back(
						0.5 * hexahedra_mesh_.vertices[hexahedra_mesh_.faces[i][hex_face_vertices[need_refine[i]][5]]] +
						0.5 * face_center);
					optimized_hexahedra_mesh_.vertices.push_back(
						0.5 * hexahedra_mesh_.vertices[hexahedra_mesh_.faces[i][hex_face_vertices[need_refine[i]][6]]] +
						0.5 * face_center);
					optimized_hexahedra_mesh_.vertices.push_back(
						0.5 * hexahedra_mesh_.vertices[hexahedra_mesh_.faces[i][hex_face_vertices[need_refine[i]][7]]] +
						0.5 * face_center);
					optimized_hexahedra_mesh_.faces[i] = {
						original_vertex_number,
						original_vertex_number + 1,
						original_vertex_number + 2,
						original_vertex_number + 3,
						original_vertex_number + 4,
						original_vertex_number + 5,
						original_vertex_number + 6,
						original_vertex_number + 7 };
					optimized_hexahedra_mesh_.faces.push_back({
						original_vertex_number + 4,
						original_vertex_number + 5,
						original_vertex_number + 6,
						original_vertex_number + 7,
						hexahedra_mesh_.faces[i][hex_face_vertices[need_refine[i]][4]],
						hexahedra_mesh_.faces[i][hex_face_vertices[need_refine[i]][5]],
						hexahedra_mesh_.faces[i][hex_face_vertices[need_refine[i]][6]],
						hexahedra_mesh_.faces[i][hex_face_vertices[need_refine[i]][7]] });
					optimized_hexahedra_mesh_.faces.push_back({
						hexahedra_mesh_.faces[i][hex_face_vertices[need_refine[i]][4]],
						hexahedra_mesh_.faces[i][hex_face_vertices[need_refine[i]][5]],
						hexahedra_mesh_.faces[i][hex_face_vertices[need_refine[i]][1]],
						hexahedra_mesh_.faces[i][hex_face_vertices[need_refine[i]][0]],
						original_vertex_number + 4,
						original_vertex_number + 5,
						original_vertex_number + 1,
						original_vertex_number });
					optimized_hexahedra_mesh_.faces.push_back({
						hexahedra_mesh_.faces[i][hex_face_vertices[need_refine[i]][7]],
						hexahedra_mesh_.faces[i][hex_face_vertices[need_refine[i]][4]],
						hexahedra_mesh_.faces[i][hex_face_vertices[need_refine[i]][0]],
						hexahedra_mesh_.faces[i][hex_face_vertices[need_refine[i]][3]],
						original_vertex_number + 7,
						original_vertex_number + 4,
						original_vertex_number,
						original_vertex_number + 3 });
					optimized_hexahedra_mesh_.faces.push_back({
						hexahedra_mesh_.faces[i][hex_face_vertices[need_refine[i]][3]],
						hexahedra_mesh_.faces[i][hex_face_vertices[need_refine[i]][2]],
						hexahedra_mesh_.faces[i][hex_face_vertices[need_refine[i]][6]],
						hexahedra_mesh_.faces[i][hex_face_vertices[need_refine[i]][7]],
						original_vertex_number + 3,
						original_vertex_number + 2,
						original_vertex_number + 6,
						original_vertex_number + 7 });
					optimized_hexahedra_mesh_.faces.push_back({
						hexahedra_mesh_.faces[i][hex_face_vertices[need_refine[i]][5]],
						hexahedra_mesh_.faces[i][hex_face_vertices[need_refine[i]][6]],
						hexahedra_mesh_.faces[i][hex_face_vertices[need_refine[i]][2]],
						hexahedra_mesh_.faces[i][hex_face_vertices[need_refine[i]][1]],
						original_vertex_number + 5,
						original_vertex_number + 6,
						original_vertex_number + 2,
						original_vertex_number + 1 });

					std::vector<size_t> face_array = {
						hexahedra_mesh_.faces[i][hex_face_vertices[need_refine[i]][0]],
						hexahedra_mesh_.faces[i][hex_face_vertices[need_refine[i]][1]],
						original_vertex_number + 1,
						original_vertex_number
					};
					faces_on_surface_set[HASH_FACE(face_array)] = face_array;
					face_array = {
						hexahedra_mesh_.faces[i][hex_face_vertices[need_refine[i]][1]],
						hexahedra_mesh_.faces[i][hex_face_vertices[need_refine[i]][2]],
						original_vertex_number + 2,
						original_vertex_number + 1
					};
					faces_on_surface_set[HASH_FACE(face_array)] = face_array;
					face_array = {
						hexahedra_mesh_.faces[i][hex_face_vertices[need_refine[i]][2]],
						hexahedra_mesh_.faces[i][hex_face_vertices[need_refine[i]][3]],
						original_vertex_number + 3,
						original_vertex_number + 2
					};
					faces_on_surface_set[HASH_FACE(face_array)] = face_array;
					face_array = {
						hexahedra_mesh_.faces[i][hex_face_vertices[need_refine[i]][3]],
						hexahedra_mesh_.faces[i][hex_face_vertices[need_refine[i]][0]],
						original_vertex_number,
						original_vertex_number + 3
					};
					faces_on_surface_set[HASH_FACE(face_array)] = face_array;
					face_array = {
						original_vertex_number,
						original_vertex_number + 1,
						original_vertex_number + 2,
						original_vertex_number + 3
					};
					faces_on_surface_set[HASH_FACE(face_array)] = face_array;
					original_vertex_number += 8;
				}
			}
			// finish sharp feature handling

			std::unordered_map<size_t, size_t> vertices_on_surface_map;
			size_t vertex_number = original_vertex_number;

			for (const auto& i : faces_on_surface_set) {
				const std::vector<size_t>& face_array = { i.second[0], i.second[1], i.second[2], i.second[3] };
				for (const auto& face_vertices : face_array) {
					if (!vertices_on_surface_map.count(face_vertices)) {
						vertices_on_surface_map[face_vertices] = vertex_number;
						optimized_hexahedra_mesh_.vertices.push_back(optimized_hexahedra_mesh_.vertices[face_vertices]);
						++vertex_number;
					}
				}

				optimized_hexahedra_mesh_.faces.push_back(
					{ vertices_on_surface_map[face_array[0]], vertices_on_surface_map[face_array[1]],
					vertices_on_surface_map[face_array[2]], vertices_on_surface_map[face_array[3]],
					face_array[0], face_array[1], face_array[2], face_array[3] });
			}
			size_t element_number = optimized_hexahedra_mesh_.faces.size();
			// update  where_should_end_vertex_go and where_should_middle_vertex_go
			where_should_end_vertex_go.resize(vertex_number, -2);
			where_should_middle_vertex_go.resize(vertex_number);
			for (size_t i = 0; i < original_vertex_number; ++i) {
				if (where_should_end_vertex_go[i] > -2) {
					where_should_end_vertex_go[vertices_on_surface_map[i]] = where_should_end_vertex_go[i];
					if (where_should_end_vertex_go[i] == -1) {
						where_should_middle_vertex_go[i].first = vertices_on_surface_map[where_should_middle_vertex_go[i].first];
						where_should_middle_vertex_go[i].second = vertices_on_surface_map[where_should_middle_vertex_go[i].second];
						where_should_middle_vertex_go[vertices_on_surface_map[i]] = where_should_middle_vertex_go[i];
					}
				}
			}

			// right hand normal indices around all the surface vertices
			std::vector<std::vector<std::vector<size_t>>> surface_vertex_normals(vertex_number);
			for (const auto& face : optimized_hexahedra_mesh_.faces) {
				for (size_t face_vertex = 0; face_vertex < 4; ++face_vertex) {
					surface_vertex_normals[face[face_vertex]].push_back(
						{ face[(face_vertex + 3) % 4],
						face[(face_vertex + 1) % 4],
						face[(face_vertex + 2) % 4] });
				}
			}

			// record vertex neighbors
			const std::vector<std::vector<size_t>> hex_vertex_neighbors =
			{ { 1, 3, 4 },
			  { 0, 2, 5 },
			  { 1, 3, 6 },
			  { 0, 2, 7 },
			  { 0, 5, 7 },
			  { 1, 4, 6 },
			  { 2, 5, 7 },
			  { 3, 4, 6 } };
			std::vector<std::unordered_set<size_t>> vertex_neighbors(vertex_number);
			for (const auto& i : optimized_hexahedra_mesh_.faces) {
				for (size_t j = 0; j < 8; ++j) {
					for (size_t k = 0; k < 3; ++k) {
						if (i[j] < original_vertex_number || i[hex_vertex_neighbors[j][k]] >= original_vertex_number) {
							vertex_neighbors[i[j]].insert(i[hex_vertex_neighbors[j][k]]);
						}
					}
				}
			}

			std::vector<std::vector<size_t>> elements_around_vertex(vertex_number);
			for (long long hex_index = 0; hex_index < element_number; ++hex_index) {
				for (size_t i = 0; i < 8; ++i) {
					elements_around_vertex[optimized_hexahedra_mesh_.faces[hex_index][i]].push_back(hex_index);
				}
			}

			std::vector<MeshTypes::Vec3> triangle_normal(triangle_mesh_.faces.size());
#pragma omp parallel for
			for (long long triangle_index = 0; triangle_index < triangle_mesh_.faces.size(); ++triangle_index) {
				MeshTypes::Vec3 v0 = triangle_mesh_.vertices[triangle_mesh_.faces[triangle_index][0]];
				MeshTypes::Vec3 v1 = triangle_mesh_.vertices[triangle_mesh_.faces[triangle_index][1]];
				MeshTypes::Vec3 v2 = triangle_mesh_.vertices[triangle_mesh_.faces[triangle_index][2]];
				triangle_normal[triangle_index] = (v1 - v0).cross(v2 - v0);
			}

			// get a vertex that is outside the triangle mesh
			MeshTypes::Vec3 triangle_low_bound = triangle_mesh_.vertices[triangle_mesh_.faces[0][0]];
			MeshTypes::Vec3 triangle_up_bound = triangle_mesh_.vertices[triangle_mesh_.faces[0][0]];
			for (const auto& face : triangle_mesh_.faces) {
				for (const auto& vertex : face) {
					triangle_low_bound.x = std::fmin(triangle_low_bound.x, triangle_mesh_.vertices[vertex].x);
					triangle_low_bound.y = std::fmin(triangle_low_bound.y, triangle_mesh_.vertices[vertex].y);
					triangle_low_bound.z = std::fmin(triangle_low_bound.z, triangle_mesh_.vertices[vertex].z);
					triangle_up_bound.x = std::fmax(triangle_up_bound.x, triangle_mesh_.vertices[vertex].x);
					triangle_up_bound.y = std::fmax(triangle_up_bound.y, triangle_mesh_.vertices[vertex].y);
					triangle_up_bound.z = std::fmax(triangle_up_bound.z, triangle_mesh_.vertices[vertex].z);
				}
			}
			const MeshTypes::Vec3 triangle_middle = 0.5 * (triangle_low_bound + triangle_up_bound);
			triangle_low_bound = 1.1 * triangle_low_bound - triangle_middle;
			triangle_up_bound = 1.1 * triangle_up_bound - triangle_middle;

			// variable
			std::vector<MeshTypes::Vec3> gradient(vertex_number);
			size_t below_threshold_count;
			std::vector<double> vertex_minimum_volume(vertex_number);
			std::vector<MeshTypes::Vec3> closest_projection_vertex(vertex_number, triangle_low_bound);
			double next_smooth_threshold = std::numeric_limits<double>::max();
			size_t smooth_count = 0;

			// optimization
			while (true) {
				// zero all gradients
				std::fill(gradient.begin(), gradient.end(), MeshTypes::Vec3());

				// traverse all elements, update mesh quality gradients
				below_threshold_count = 0;
#pragma omp parallel for reduction(+: below_threshold_count)
				for (long long element_index = 0; element_index < element_number; ++element_index) {
					const MeshTypes::Vec3& p0 = optimized_hexahedra_mesh_.vertices[optimized_hexahedra_mesh_.faces[element_index][0]];
					const MeshTypes::Vec3& p1 = optimized_hexahedra_mesh_.vertices[optimized_hexahedra_mesh_.faces[element_index][1]];
					const MeshTypes::Vec3& p2 = optimized_hexahedra_mesh_.vertices[optimized_hexahedra_mesh_.faces[element_index][2]];
					const MeshTypes::Vec3& p3 = optimized_hexahedra_mesh_.vertices[optimized_hexahedra_mesh_.faces[element_index][3]];
					const MeshTypes::Vec3& p4 = optimized_hexahedra_mesh_.vertices[optimized_hexahedra_mesh_.faces[element_index][4]];
					const MeshTypes::Vec3& p5 = optimized_hexahedra_mesh_.vertices[optimized_hexahedra_mesh_.faces[element_index][5]];
					const MeshTypes::Vec3& p6 = optimized_hexahedra_mesh_.vertices[optimized_hexahedra_mesh_.faces[element_index][6]];
					const MeshTypes::Vec3& p7 = optimized_hexahedra_mesh_.vertices[optimized_hexahedra_mesh_.faces[element_index][7]];

					const size_t below_threshold_vertex = MeshQualityBelowThreshold(p0, p1, p2, p3, p4, p5, p6, p7);

					// mesh quality is above threshold
					if (below_threshold_vertex < 9) {
						++below_threshold_count;
						const double scaling_factor =
							((p0 - p1).norm() + (p1 - p2).norm() + (p2 - p3).norm() + (p3 - p0).norm() +
								(p0 - p4).norm() + (p1 - p5).norm() + (p2 - p6).norm() + (p3 - p7).norm() +
								(p4 - p5).norm() + (p5 - p6).norm() + (p6 - p7).norm() + (p7 - p4).norm()) / 12.0;

						// 0
						if (below_threshold_vertex == 0) {
							gradient[optimized_hexahedra_mesh_.faces[element_index][0]].x += -(p4.y * p1.z + p5.y * p1.z + p1.y * p2.z + p5.y * p2.z - p7.y * p2.z + p1.y * p3.z - p4.y * p3.z - p7.y * p3.z -
								p1.y * p4.z - p5.y * p4.z + p7.y * p4.z - p1.y * p5.z + p4.y * p5.z + p7.y * p5.z - p4.y * p7.z - p5.y * p7.z +
								p3.y * (-p1.z - p2.z + p4.z + p7.z) + p2.y * (-p1.z + p3.z - p5.z + p7.z)) / 16.0 / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][0]].y += (p4.x * p1.z + p5.x * p1.z + p1.x * p2.z + p5.x * p2.z - p7.x * p2.z + p1.x * p3.z - p4.x * p3.z - p7.x * p3.z -
								p1.x * p4.z - p5.x * p4.z + p7.x * p4.z - p1.x * p5.z + p4.x * p5.z + p7.x * p5.z - p4.x * p7.z - p5.x * p7.z +
								p3.x * (-p1.z - p2.z + p4.z + p7.z) + p2.x * (-p1.z + p3.z - p5.z + p7.z)) / 16.0 / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][0]].z += -(p4.x * p1.y + p5.x * p1.y + p1.x * p2.y + p5.x * p2.y - p7.x * p2.y + p1.x * p3.y - p4.x * p3.y - p7.x * p3.y -
								p1.x * p4.y - p5.x * p4.y + p7.x * p4.y - p1.x * p5.y + p4.x * p5.y + p7.x * p5.y - p4.x * p7.y - p5.x * p7.y +
								p3.x * (-p1.y - p2.y + p4.y + p7.y) + p2.x * (-p1.y + p3.y - p5.y + p7.y)) / 16.0 / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][1]].x += (p4.y * p0.z + p5.y * p0.z + p0.y * p2.z - p5.y * p2.z - p6.y * p2.z + p0.y * p3.z + p4.y * p3.z - p6.y * p3.z -
								p0.y * p4.z + p5.y * p4.z + p6.y * p4.z - p0.y * p5.z - p4.y * p5.z + p6.y * p5.z - p4.y * p6.z - p5.y * p6.z +
								p3.y * (-p0.z + p2.z - p4.z + p6.z) + p2.y * (-p0.z - p3.z + p5.z + p6.z)) / 16.0 / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][1]].y += -(p4.x * p0.z + p5.x * p0.z + p0.x * p2.z - p5.x * p2.z - p6.x * p2.z + p0.x * p3.z + p4.x * p3.z - p6.x * p3.z -
								p0.x * p4.z + p5.x * p4.z + p6.x * p4.z - p0.x * p5.z - p4.x * p5.z + p6.x * p5.z - p4.x * p6.z - p5.x * p6.z +
								p3.x * (-p0.z + p2.z - p4.z + p6.z) + p2.x * (-p0.z - p3.z + p5.z + p6.z)) / 16.0 / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][1]].z += (p4.x * p0.y + p5.x * p0.y + p0.x * p2.y - p5.x * p2.y - p6.x * p2.y + p0.x * p3.y + p4.x * p3.y - p6.x * p3.y -
								p0.x * p4.y + p5.x * p4.y + p6.x * p4.y - p0.x * p5.y - p4.x * p5.y + p6.x * p5.y - p4.x * p6.y - p5.x * p6.y +
								p3.x * (-p0.y + p2.y - p4.y + p6.y) + p2.x * (-p0.y - p3.y + p5.y + p6.y)) / 16.0 / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][2]].x += -(-p5.y * p0.z + p7.y * p0.z + p0.y * p1.z - p5.y * p1.z - p6.y * p1.z - p0.y * p3.z + p6.y * p3.z + p7.y * p3.z +
								p0.y * p5.z - p6.y * p5.z - p7.y * p5.z + p5.y * p6.z - p7.y * p6.z + p1.y * (-p0.z - p3.z + p5.z + p6.z) +
								p3.y * (p0.z + p1.z - p6.z - p7.z) - p0.y * p7.z + p5.y * p7.z + p6.y * p7.z) / 16.0 / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][2]].y += (-p5.x * p0.z + p7.x * p0.z + p0.x * p1.z - p5.x * p1.z - p6.x * p1.z - p0.x * p3.z + p6.x * p3.z + p7.x * p3.z +
								p0.x * p5.z - p6.x * p5.z - p7.x * p5.z + p5.x * p6.z - p7.x * p6.z + p1.x * (-p0.z - p3.z + p5.z + p6.z) +
								p3.x * (p0.z + p1.z - p6.z - p7.z) - p0.x * p7.z + p5.x * p7.z + p6.x * p7.z) / 16.0 / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][2]].z += -(-p5.x * p0.y + p7.x * p0.y + p0.x * p1.y - p5.x * p1.y - p6.x * p1.y - p0.x * p3.y + p6.x * p3.y + p7.x * p3.y +
								p0.x * p5.y - p6.x * p5.y - p7.x * p5.y + p5.x * p6.y - p7.x * p6.y + p1.x * (-p0.y - p3.y + p5.y + p6.y) +
								p3.x * (p0.y + p1.y - p6.y - p7.y) - p0.x * p7.y + p5.x * p7.y + p6.x * p7.y) / 16.0 / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][3]].x += -(p4.y * p0.z + p7.y * p0.z + p0.y * p1.z + p4.y * p1.z - p6.y * p1.z + p0.y * p2.z - p6.y * p2.z - p7.y * p2.z -
								p0.y * p4.z + p6.y * p4.z + p7.y * p4.z - p4.y * p6.z - p7.y * p6.z + p1.y * (-p0.z + p2.z - p4.z + p6.z) -
								p0.y * p7.z - p4.y * p7.z + p6.y * p7.z + p2.y * (-p0.z - p1.z + p6.z + p7.z)) / 16.0 / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][3]].y += (p4.x * p0.z + p7.x * p0.z + p0.x * p1.z + p4.x * p1.z - p6.x * p1.z + p0.x * p2.z - p6.x * p2.z - p7.x * p2.z -
								p0.x * p4.z + p6.x * p4.z + p7.x * p4.z - p4.x * p6.z - p7.x * p6.z + p1.x * (-p0.z + p2.z - p4.z + p6.z) -
								p0.x * p7.z - p4.x * p7.z + p6.x * p7.z + p2.x * (-p0.z - p1.z + p6.z + p7.z)) / 16.0 / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][3]].z += -(p4.x * p0.y + p7.x * p0.y + p0.x * p1.y + p4.x * p1.y - p6.x * p1.y + p0.x * p2.y - p6.x * p2.y - p7.x * p2.y -
								p0.x * p4.y + p6.x * p4.y + p7.x * p4.y - p4.x * p6.y - p7.x * p6.y + p1.x * (-p0.y + p2.y - p4.y + p6.y) -
								p0.x * p7.y - p4.x * p7.y + p6.x * p7.y + p2.x * (-p0.y - p1.y + p6.y + p7.y)) / 16.0 / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][4]].x += (-p5.y * p0.z + p7.y * p0.z + p0.y * p1.z - p5.y * p1.z - p6.y * p1.z - p0.y * p3.z + p6.y * p3.z + p7.y * p3.z +
								p0.y * p5.z - p6.y * p5.z - p7.y * p5.z + p5.y * p6.z - p7.y * p6.z + p1.y * (-p0.z - p3.z + p5.z + p6.z) +
								p3.y * (p0.z + p1.z - p6.z - p7.z) - p0.y * p7.z + p5.y * p7.z + p6.y * p7.z) / 16.0 / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][4]].y += -(-p5.x * p0.z + p7.x * p0.z + p0.x * p1.z - p5.x * p1.z - p6.x * p1.z - p0.x * p3.z + p6.x * p3.z + p7.x * p3.z +
								p0.x * p5.z - p6.x * p5.z - p7.x * p5.z + p5.x * p6.z - p7.x * p6.z + p1.x * (-p0.z - p3.z + p5.z + p6.z) +
								p3.x * (p0.z + p1.z - p6.z - p7.z) - p0.x * p7.z + p5.x * p7.z + p6.x * p7.z) / 16.0 / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][4]].z += (-p5.x * p0.y + p7.x * p0.y + p0.x * p1.y - p5.x * p1.y - p6.x * p1.y - p0.x * p3.y + p6.x * p3.y + p7.x * p3.y +
								p0.x * p5.y - p6.x * p5.y - p7.x * p5.y + p5.x * p6.y - p7.x * p6.y + p1.x * (-p0.y - p3.y + p5.y + p6.y) +
								p3.x * (p0.y + p1.y - p6.y - p7.y) - p0.x * p7.y + p5.x * p7.y + p6.x * p7.y) / 16.0 / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][5]].x += (p4.y * p0.z + p7.y * p0.z + p0.y * p1.z + p4.y * p1.z - p6.y * p1.z + p0.y * p2.z - p6.y * p2.z - p7.y * p2.z -
								p0.y * p4.z + p6.y * p4.z + p7.y * p4.z - p4.y * p6.z - p7.y * p6.z + p1.y * (-p0.z + p2.z - p4.z + p6.z) -
								p0.y * p7.z - p4.y * p7.z + p6.y * p7.z + p2.y * (-p0.z - p1.z + p6.z + p7.z)) / 16.0 / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][5]].y += -(p4.x * p0.z + p7.x * p0.z + p0.x * p1.z + p4.x * p1.z - p6.x * p1.z + p0.x * p2.z - p6.x * p2.z - p7.x * p2.z -
								p0.x * p4.z + p6.x * p4.z + p7.x * p4.z - p4.x * p6.z - p7.x * p6.z + p1.x * (-p0.z + p2.z - p4.z + p6.z) -
								p0.x * p7.z - p4.x * p7.z + p6.x * p7.z + p2.x * (-p0.z - p1.z + p6.z + p7.z)) / 16.0 / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][5]].z += (p4.x * p0.y + p7.x * p0.y + p0.x * p1.y + p4.x * p1.y - p6.x * p1.y + p0.x * p2.y - p6.x * p2.y - p7.x * p2.y -
								p0.x * p4.y + p6.x * p4.y + p7.x * p4.y - p4.x * p6.y - p7.x * p6.y + p1.x * (-p0.y + p2.y - p4.y + p6.y) -
								p0.x * p7.y - p4.x * p7.y + p6.x * p7.y + p2.x * (-p0.y - p1.y + p6.y + p7.y)) / 16.0 / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][6]].x += (p4.y * p1.z + p5.y * p1.z + p1.y * p2.z + p5.y * p2.z - p7.y * p2.z + p1.y * p3.z - p4.y * p3.z - p7.y * p3.z -
								p1.y * p4.z - p5.y * p4.z + p7.y * p4.z - p1.y * p5.z + p4.y * p5.z + p7.y * p5.z - p4.y * p7.z - p5.y * p7.z +
								p3.y * (-p1.z - p2.z + p4.z + p7.z) + p2.y * (-p1.z + p3.z - p5.z + p7.z)) / 16.0 / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][6]].y += -(p4.x * p1.z + p5.x * p1.z + p1.x * p2.z + p5.x * p2.z - p7.x * p2.z + p1.x * p3.z - p4.x * p3.z - p7.x * p3.z -
								p1.x * p4.z - p5.x * p4.z + p7.x * p4.z - p1.x * p5.z + p4.x * p5.z + p7.x * p5.z - p4.x * p7.z - p5.x * p7.z +
								p3.x * (-p1.z - p2.z + p4.z + p7.z) + p2.x * (-p1.z + p3.z - p5.z + p7.z)) / 16.0 / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][6]].z += (p4.x * p1.y + p5.x * p1.y + p1.x * p2.y + p5.x * p2.y - p7.x * p2.y + p1.x * p3.y - p4.x * p3.y - p7.x * p3.y -
								p1.x * p4.y - p5.x * p4.y + p7.x * p4.y - p1.x * p5.y + p4.x * p5.y + p7.x * p5.y - p4.x * p7.y - p5.x * p7.y +
								p3.x * (-p1.y - p2.y + p4.y + p7.y) + p2.x * (-p1.y + p3.y - p5.y + p7.y)) / 16.0 / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][7]].x += -(p4.y * p0.z + p5.y * p0.z + p0.y * p2.z - p5.y * p2.z - p6.y * p2.z + p0.y * p3.z + p4.y * p3.z - p6.y * p3.z -
								p0.y * p4.z + p5.y * p4.z + p6.y * p4.z - p0.y * p5.z - p4.y * p5.z + p6.y * p5.z - p4.y * p6.z - p5.y * p6.z +
								p3.y * (-p0.z + p2.z - p4.z + p6.z) + p2.y * (-p0.z - p3.z + p5.z + p6.z)) / 16.0 / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][7]].y += (p4.x * p0.z + p5.x * p0.z + p0.x * p2.z - p5.x * p2.z - p6.x * p2.z + p0.x * p3.z + p4.x * p3.z - p6.x * p3.z -
								p0.x * p4.z + p5.x * p4.z + p6.x * p4.z - p0.x * p5.z - p4.x * p5.z + p6.x * p5.z - p4.x * p6.z - p5.x * p6.z +
								p3.x * (-p0.z + p2.z - p4.z + p6.z) + p2.x * (-p0.z - p3.z + p5.z + p6.z)) / 16.0 / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][7]].z += -(p4.x * p0.y + p5.x * p0.y + p0.x * p2.y - p5.x * p2.y - p6.x * p2.y + p0.x * p3.y + p4.x * p3.y - p6.x * p3.y -
								p0.x * p4.y + p5.x * p4.y + p6.x * p4.y - p0.x * p5.y - p4.x * p5.y + p6.x * p5.y - p4.x * p6.y - p5.x * p6.y +
								p3.x * (-p0.y + p2.y - p4.y + p6.y) + p2.x * (-p0.y - p3.y + p5.y + p6.y)) / 16.0 / scaling_factor;
						}
						// 1
						else if (below_threshold_vertex == 1) {
							gradient[optimized_hexahedra_mesh_.faces[element_index][0]].x += (p4.y * (-p1.z + p3.z) + p3.y * (p1.z - p4.z) + p1.y * (-p3.z + p4.z)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][0]].y += (p4.x * (p1.z - p3.z) + p1.x * (p3.z - p4.z) + p3.x * (-p1.z + p4.z)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][0]].z += (p4.x * (-p1.y + p3.y) + p3.x * (p1.y - p4.y) + p1.x * (-p3.y + p4.y)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][1]].x += (p4.y * (p0.z - p3.z) + p0.y * (p3.z - p4.z) + p3.y * (-p0.z + p4.z)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][1]].y += (p4.x * (-p0.z + p3.z) + p3.x * (p0.z - p4.z) + p0.x * (-p3.z + p4.z)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][1]].z += (p4.x * (p0.y - p3.y) + p0.x * (p3.y - p4.y) + p3.x * (-p0.y + p4.y)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][3]].x += (p4.y * (-p0.z + p1.z) + p1.y * (p0.z - p4.z) + p0.y * (-p1.z + p4.z)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][3]].y += (p4.x * (p0.z - p1.z) + p0.x * (p1.z - p4.z) + p1.x * (-p0.z + p4.z)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][3]].z += (p4.x * (-p0.y + p1.y) + p1.x * (p0.y - p4.y) + p0.x * (-p1.y + p4.y)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][4]].x += (p3.y * (p0.z - p1.z) + p0.y * (p1.z - p3.z) + p1.y * (-p0.z + p3.z)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][4]].y += (p3.x * (-p0.z + p1.z) + p1.x * (p0.z - p3.z) + p0.x * (-p1.z + p3.z)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][4]].z += (p3.x * (p0.y - p1.y) + p0.x * (p1.y - p3.y) + p1.x * (-p0.y + p3.y)) / scaling_factor;
						}
						// 2
						else if (below_threshold_vertex == 2) {
							gradient[optimized_hexahedra_mesh_.faces[element_index][0]].x += (p5.y * (-p1.z + p2.z) + p2.y * (p1.z - p5.z) + p1.y * (-p2.z + p5.z)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][0]].y += (p5.x * (p1.z - p2.z) + p1.x * (p2.z - p5.z) + p2.x * (-p1.z + p5.z)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][0]].z += (p5.x * (-p1.y + p2.y) + p2.x * (p1.y - p5.y) + p1.x * (-p2.y + p5.y)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][1]].x += (p5.y * (p0.z - p2.z) + p0.y * (p2.z - p5.z) + p2.y * (-p0.z + p5.z)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][1]].y += (p5.x * (-p0.z + p2.z) + p2.x * (p0.z - p5.z) + p0.x * (-p2.z + p5.z)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][1]].z += (p5.x * (p0.y - p2.y) + p0.x * (p2.y - p5.y) + p2.x * (-p0.y + p5.y)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][2]].x += (p5.y * (-p0.z + p1.z) + p1.y * (p0.z - p5.z) + p0.y * (-p1.z + p5.z)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][2]].y += (p5.x * (p0.z - p1.z) + p0.x * (p1.z - p5.z) + p1.x * (-p0.z + p5.z)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][2]].z += (p5.x * (-p0.y + p1.y) + p1.x * (p0.y - p5.y) + p0.x * (-p1.y + p5.y)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][5]].x += (p2.y * (p0.z - p1.z) + p0.y * (p1.z - p2.z) + p1.y * (-p0.z + p2.z)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][5]].y += (p2.x * (-p0.z + p1.z) + p1.x * (p0.z - p2.z) + p0.x * (-p1.z + p2.z)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][5]].z += (p2.x * (p0.y - p1.y) + p0.x * (p1.y - p2.y) + p1.x * (-p0.y + p2.y)) / scaling_factor;
						}
						// 3
						else if (below_threshold_vertex == 3) {
							gradient[optimized_hexahedra_mesh_.faces[element_index][1]].x += (p6.y * (-p2.z + p3.z) + p3.y * (p2.z - p6.z) + p2.y * (-p3.z + p6.z)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][1]].y += (p6.x * (p2.z - p3.z) + p2.x * (p3.z - p6.z) + p3.x * (-p2.z + p6.z)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][1]].z += (p6.x * (-p2.y + p3.y) + p3.x * (p2.y - p6.y) + p2.x * (-p3.y + p6.y)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][2]].x += (p6.y * (p1.z - p3.z) + p1.y * (p3.z - p6.z) + p3.y * (-p1.z + p6.z)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][2]].y += (p6.x * (-p1.z + p3.z) + p3.x * (p1.z - p6.z) + p1.x * (-p3.z + p6.z)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][2]].z += (p6.x * (p1.y - p3.y) + p1.x * (p3.y - p6.y) + p3.x * (-p1.y + p6.y)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][3]].x += (p6.y * (-p1.z + p2.z) + p2.y * (p1.z - p6.z) + p1.y * (-p2.z + p6.z)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][3]].y += (p6.x * (p1.z - p2.z) + p1.x * (p2.z - p6.z) + p2.x * (-p1.z + p6.z)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][3]].z += (p6.x * (-p1.y + p2.y) + p2.x * (p1.y - p6.y) + p1.x * (-p2.y + p6.y)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][6]].x += (p3.y * (p1.z - p2.z) + p1.y * (p2.z - p3.z) + p2.y * (-p1.z + p3.z)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][6]].y += (p3.x * (-p1.z + p2.z) + p2.x * (p1.z - p3.z) + p1.x * (-p2.z + p3.z)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][6]].z += (p3.x * (p1.y - p2.y) + p1.x * (p2.y - p3.y) + p2.x * (-p1.y + p3.y)) / scaling_factor;
						}
						// 4
						else if (below_threshold_vertex == 4) {
							gradient[optimized_hexahedra_mesh_.faces[element_index][0]].x += (p7.y * (-p2.z + p3.z) + p3.y * (p2.z - p7.z) + p2.y * (-p3.z + p7.z)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][0]].y += (p7.x * (p2.z - p3.z) + p2.x * (p3.z - p7.z) + p3.x * (-p2.z + p7.z)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][0]].z += (p7.x * (-p2.y + p3.y) + p3.x * (p2.y - p7.y) + p2.x * (-p3.y + p7.y)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][2]].x += (p7.y * (p0.z - p3.z) + p0.y * (p3.z - p7.z) + p3.y * (-p0.z + p7.z)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][2]].y += (p7.x * (-p0.z + p3.z) + p3.x * (p0.z - p7.z) + p0.x * (-p3.z + p7.z)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][2]].z += (p7.x * (p0.y - p3.y) + p0.x * (p3.y - p7.y) + p3.x * (-p0.y + p7.y)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][3]].x += (p7.y * (-p0.z + p2.z) + p2.y * (p0.z - p7.z) + p0.y * (-p2.z + p7.z)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][3]].y += (p7.x * (p0.z - p2.z) + p0.x * (p2.z - p7.z) + p2.x * (-p0.z + p7.z)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][3]].z += (p7.x * (-p0.y + p2.y) + p2.x * (p0.y - p7.y) + p0.x * (-p2.y + p7.y)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][7]].x += (p3.y * (p0.z - p2.z) + p0.y * (p2.z - p3.z) + p2.y * (-p0.z + p3.z)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][7]].y += (p3.x * (-p0.z + p2.z) + p2.x * (p0.z - p3.z) + p0.x * (-p2.z + p3.z)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][7]].z += (p3.x * (p0.y - p2.y) + p0.x * (p2.y - p3.y) + p2.x * (-p0.y + p3.y)) / scaling_factor;
						}
						// 5
						else if (below_threshold_vertex == 5) {
							gradient[optimized_hexahedra_mesh_.faces[element_index][0]].x += (p7.y * (-p4.z + p5.z) + p5.y * (p4.z - p7.z) + p4.y * (-p5.z + p7.z)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][0]].y += (p7.x * (p4.z - p5.z) + p4.x * (p5.z - p7.z) + p5.x * (-p4.z + p7.z)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][0]].z += (p7.x * (-p4.y + p5.y) + p5.x * (p4.y - p7.y) + p4.x * (-p5.y + p7.y)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][4]].x += (p7.y * (p0.z - p5.z) + p0.y * (p5.z - p7.z) + p5.y * (-p0.z + p7.z)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][4]].y += (p7.x * (-p0.z + p5.z) + p5.x * (p0.z - p7.z) + p0.x * (-p5.z + p7.z)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][4]].z += (p7.x * (p0.y - p5.y) + p0.x * (p5.y - p7.y) + p5.x * (-p0.y + p7.y)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][5]].x += (p7.y * (-p0.z + p4.z) + p4.y * (p0.z - p7.z) + p0.y * (-p4.z + p7.z)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][5]].y += (p7.x * (p0.z - p4.z) + p0.x * (p4.z - p7.z) + p4.x * (-p0.z + p7.z)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][5]].z += (p7.x * (-p0.y + p4.y) + p4.x * (p0.y - p7.y) + p0.x * (-p4.y + p7.y)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][7]].x += (p5.y * (p0.z - p4.z) + p0.y * (p4.z - p5.z) + p4.y * (-p0.z + p5.z)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][7]].y += (p5.x * (-p0.z + p4.z) + p4.x * (p0.z - p5.z) + p0.x * (-p4.z + p5.z)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][7]].z += (p5.x * (p0.y - p4.y) + p0.x * (p4.y - p5.y) + p4.x * (-p0.y + p5.y)) / scaling_factor;
						}
						// 6
						else if (below_threshold_vertex == 6) {
							gradient[optimized_hexahedra_mesh_.faces[element_index][1]].x += (p6.y * (-p4.z + p5.z) + p5.y * (p4.z - p6.z) + p4.y * (-p5.z + p6.z)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][1]].y += (p6.x * (p4.z - p5.z) + p4.x * (p5.z - p6.z) + p5.x * (-p4.z + p6.z)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][1]].z += (p6.x * (-p4.y + p5.y) + p5.x * (p4.y - p6.y) + p4.x * (-p5.y + p6.y)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][4]].x += (p6.y * (p1.z - p5.z) + p1.y * (p5.z - p6.z) + p5.y * (-p1.z + p6.z)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][4]].y += (p6.x * (-p1.z + p5.z) + p5.x * (p1.z - p6.z) + p1.x * (-p5.z + p6.z)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][4]].z += (p6.x * (p1.y - p5.y) + p1.x * (p5.y - p6.y) + p5.x * (-p1.y + p6.y)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][5]].x += (p6.y * (-p1.z + p4.z) + p4.y * (p1.z - p6.z) + p1.y * (-p4.z + p6.z)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][5]].y += (p6.x * (p1.z - p4.z) + p1.x * (p4.z - p6.z) + p4.x * (-p1.z + p6.z)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][5]].z += (p6.x * (-p1.y + p4.y) + p4.x * (p1.y - p6.y) + p1.x * (-p4.y + p6.y)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][6]].x += (p5.y * (p1.z - p4.z) + p1.y * (p4.z - p5.z) + p4.y * (-p1.z + p5.z)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][6]].y += (p5.x * (-p1.z + p4.z) + p4.x * (p1.z - p5.z) + p1.x * (-p4.z + p5.z)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][6]].z += (p5.x * (p1.y - p4.y) + p1.x * (p4.y - p5.y) + p4.x * (-p1.y + p5.y)) / scaling_factor;
						}
						// 7
						else if (below_threshold_vertex == 7) {
							gradient[optimized_hexahedra_mesh_.faces[element_index][2]].x += (p7.y * (-p5.z + p6.z) + p6.y * (p5.z - p7.z) + p5.y * (-p6.z + p7.z)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][2]].y += (p7.x * (p5.z - p6.z) + p5.x * (p6.z - p7.z) + p6.x * (-p5.z + p7.z)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][2]].z += (p7.x * (-p5.y + p6.y) + p6.x * (p5.y - p7.y) + p5.x * (-p6.y + p7.y)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][5]].x += (p7.y * (p2.z - p6.z) + p2.y * (p6.z - p7.z) + p6.y * (-p2.z + p7.z)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][5]].y += (p7.x * (-p2.z + p6.z) + p6.x * (p2.z - p7.z) + p2.x * (-p6.z + p7.z)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][5]].z += (p7.x * (p2.y - p6.y) + p2.x * (p6.y - p7.y) + p6.x * (-p2.y + p7.y)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][6]].x += (p7.y * (-p2.z + p5.z) + p5.y * (p2.z - p7.z) + p2.y * (-p5.z + p7.z)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][6]].y += (p7.x * (p2.z - p5.z) + p2.x * (p5.z - p7.z) + p5.x * (-p2.z + p7.z)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][6]].z += (p7.x * (-p2.y + p5.y) + p5.x * (p2.y - p7.y) + p2.x * (-p5.y + p7.y)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][7]].x += (p6.y * (p2.z - p5.z) + p2.y * (p5.z - p6.z) + p5.y * (-p2.z + p6.z)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][7]].y += (p6.x * (-p2.z + p5.z) + p5.x * (p2.z - p6.z) + p2.x * (-p5.z + p6.z)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][7]].z += (p6.x * (p2.y - p5.y) + p2.x * (p5.y - p6.y) + p5.x * (-p2.y + p6.y)) / scaling_factor;
						}
						// 8
						else {
							gradient[optimized_hexahedra_mesh_.faces[element_index][3]].x += (p7.y * (-p4.z + p6.z) + p6.y * (p4.z - p7.z) + p4.y * (-p6.z + p7.z)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][3]].y += (p7.x * (p4.z - p6.z) + p4.x * (p6.z - p7.z) + p6.x * (-p4.z + p7.z)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][3]].z += (p7.x * (-p4.y + p6.y) + p6.x * (p4.y - p7.y) + p4.x * (-p6.y + p7.y)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][4]].x += (p7.y * (p3.z - p6.z) + p3.y * (p6.z - p7.z) + p6.y * (-p3.z + p7.z)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][4]].y += (p7.x * (-p3.z + p6.z) + p6.x * (p3.z - p7.z) + p3.x * (-p6.z + p7.z)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][4]].z += (p7.x * (p3.y - p6.y) + p3.x * (p6.y - p7.y) + p6.x * (-p3.y + p7.y)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][6]].x += (p7.y * (-p3.z + p4.z) + p4.y * (p3.z - p7.z) + p3.y * (-p4.z + p7.z)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][6]].y += (p7.x * (p3.z - p4.z) + p3.x * (p4.z - p7.z) + p4.x * (-p3.z + p7.z)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][6]].z += (p7.x * (-p3.y + p4.y) + p4.x * (p3.y - p7.y) + p3.x * (-p4.y + p7.y)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][7]].x += (p6.y * (p3.z - p4.z) + p3.y * (p4.z - p6.z) + p4.y * (-p3.z + p6.z)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][7]].y += (p6.x * (-p3.z + p4.z) + p4.x * (p3.z - p6.z) + p3.x * (-p4.z + p6.z)) / scaling_factor;
							gradient[optimized_hexahedra_mesh_.faces[element_index][7]].z += (p6.x * (p3.y - p4.y) + p3.x * (p4.y - p6.y) + p4.x * (-p3.y + p6.y)) / scaling_factor;
						}
					}
				}

				// if all Jacobian are greater than 0
				if (below_threshold_count == 0) {
					std::fill(vertex_minimum_volume.begin(),
						vertex_minimum_volume.end(), std::numeric_limits<double>::max());

					// prevent too close vertices
#pragma omp parallel for
					for (long long element_index = 0; element_index < element_number; ++element_index) {
						const size_t& x0 = optimized_hexahedra_mesh_.faces[element_index][0];
						const size_t& x1 = optimized_hexahedra_mesh_.faces[element_index][1];
						const size_t& x2 = optimized_hexahedra_mesh_.faces[element_index][2];
						const size_t& x3 = optimized_hexahedra_mesh_.faces[element_index][3];
						const size_t& x4 = optimized_hexahedra_mesh_.faces[element_index][4];
						const size_t& x5 = optimized_hexahedra_mesh_.faces[element_index][5];
						const size_t& x6 = optimized_hexahedra_mesh_.faces[element_index][6];
						const size_t& x7 = optimized_hexahedra_mesh_.faces[element_index][7];
						const MeshTypes::Vec3& p0 = optimized_hexahedra_mesh_.vertices[x0];
						const MeshTypes::Vec3& p1 = optimized_hexahedra_mesh_.vertices[x1];
						const MeshTypes::Vec3& p2 = optimized_hexahedra_mesh_.vertices[x2];
						const MeshTypes::Vec3& p3 = optimized_hexahedra_mesh_.vertices[x3];
						const MeshTypes::Vec3& p4 = optimized_hexahedra_mesh_.vertices[x4];
						const MeshTypes::Vec3& p5 = optimized_hexahedra_mesh_.vertices[x5];
						const MeshTypes::Vec3& p6 = optimized_hexahedra_mesh_.vertices[x6];
						const MeshTypes::Vec3& p7 = optimized_hexahedra_mesh_.vertices[x7];

						const std::pair<size_t, double> return_value =
							MinMeshQuality(p0, p1, p2, p3, p4, p5, p6, p7);

						const double scaling_factor =
							((p0 - p1).norm() + (p1 - p2).norm() + (p2 - p3).norm() + (p3 - p0).norm() +
								(p0 - p4).norm() + (p1 - p5).norm() + (p2 - p6).norm() + (p3 - p7).norm() +
								(p4 - p5).norm() + (p5 - p6).norm() + (p6 - p7).norm() + (p7 - p4).norm()) / 12.0;

						// 0
						if (return_value.first == 0) {
							if (return_value.second < vertex_minimum_volume[x0]) {
								vertex_minimum_volume[x0] = return_value.second;
								gradient[x0].x = -(p4.y * p1.z + p5.y * p1.z + p1.y * p2.z + p5.y * p2.z - p7.y * p2.z + p1.y * p3.z - p4.y * p3.z - p7.y * p3.z -
									p1.y * p4.z - p5.y * p4.z + p7.y * p4.z - p1.y * p5.z + p4.y * p5.z + p7.y * p5.z - p4.y * p7.z - p5.y * p7.z +
									p3.y * (-p1.z - p2.z + p4.z + p7.z) + p2.y * (-p1.z + p3.z - p5.z + p7.z)) / 16.0 / scaling_factor;
								gradient[x0].y = (p4.x * p1.z + p5.x * p1.z + p1.x * p2.z + p5.x * p2.z - p7.x * p2.z + p1.x * p3.z - p4.x * p3.z - p7.x * p3.z -
									p1.x * p4.z - p5.x * p4.z + p7.x * p4.z - p1.x * p5.z + p4.x * p5.z + p7.x * p5.z - p4.x * p7.z - p5.x * p7.z +
									p3.x * (-p1.z - p2.z + p4.z + p7.z) + p2.x * (-p1.z + p3.z - p5.z + p7.z)) / 16.0 / scaling_factor;
								gradient[x0].z = -(p4.x * p1.y + p5.x * p1.y + p1.x * p2.y + p5.x * p2.y - p7.x * p2.y + p1.x * p3.y - p4.x * p3.y - p7.x * p3.y -
									p1.x * p4.y - p5.x * p4.y + p7.x * p4.y - p1.x * p5.y + p4.x * p5.y + p7.x * p5.y - p4.x * p7.y - p5.x * p7.y +
									p3.x * (-p1.y - p2.y + p4.y + p7.y) + p2.x * (-p1.y + p3.y - p5.y + p7.y)) / 16.0 / scaling_factor;
							}
							if (return_value.second < vertex_minimum_volume[x1]) {
								vertex_minimum_volume[x1] = return_value.second;
								gradient[x1].x = (p4.y * p0.z + p5.y * p0.z + p0.y * p2.z - p5.y * p2.z - p6.y * p2.z + p0.y * p3.z + p4.y * p3.z - p6.y * p3.z -
									p0.y * p4.z + p5.y * p4.z + p6.y * p4.z - p0.y * p5.z - p4.y * p5.z + p6.y * p5.z - p4.y * p6.z - p5.y * p6.z +
									p3.y * (-p0.z + p2.z - p4.z + p6.z) + p2.y * (-p0.z - p3.z + p5.z + p6.z)) / 16.0 / scaling_factor;
								gradient[x1].y = -(p4.x * p0.z + p5.x * p0.z + p0.x * p2.z - p5.x * p2.z - p6.x * p2.z + p0.x * p3.z + p4.x * p3.z - p6.x * p3.z -
									p0.x * p4.z + p5.x * p4.z + p6.x * p4.z - p0.x * p5.z - p4.x * p5.z + p6.x * p5.z - p4.x * p6.z - p5.x * p6.z +
									p3.x * (-p0.z + p2.z - p4.z + p6.z) + p2.x * (-p0.z - p3.z + p5.z + p6.z)) / 16.0 / scaling_factor;
								gradient[x1].z = (p4.x * p0.y + p5.x * p0.y + p0.x * p2.y - p5.x * p2.y - p6.x * p2.y + p0.x * p3.y + p4.x * p3.y - p6.x * p3.y -
									p0.x * p4.y + p5.x * p4.y + p6.x * p4.y - p0.x * p5.y - p4.x * p5.y + p6.x * p5.y - p4.x * p6.y - p5.x * p6.y +
									p3.x * (-p0.y + p2.y - p4.y + p6.y) + p2.x * (-p0.y - p3.y + p5.y + p6.y)) / 16.0 / scaling_factor;
							}
							if (return_value.second < vertex_minimum_volume[x2]) {
								vertex_minimum_volume[x2] = return_value.second;
								gradient[x2].x = -(-p5.y * p0.z + p7.y * p0.z + p0.y * p1.z - p5.y * p1.z - p6.y * p1.z - p0.y * p3.z + p6.y * p3.z + p7.y * p3.z +
									p0.y * p5.z - p6.y * p5.z - p7.y * p5.z + p5.y * p6.z - p7.y * p6.z + p1.y * (-p0.z - p3.z + p5.z + p6.z) +
									p3.y * (p0.z + p1.z - p6.z - p7.z) - p0.y * p7.z + p5.y * p7.z + p6.y * p7.z) / 16.0 / scaling_factor;
								gradient[x2].y = (-p5.x * p0.z + p7.x * p0.z + p0.x * p1.z - p5.x * p1.z - p6.x * p1.z - p0.x * p3.z + p6.x * p3.z + p7.x * p3.z +
									p0.x * p5.z - p6.x * p5.z - p7.x * p5.z + p5.x * p6.z - p7.x * p6.z + p1.x * (-p0.z - p3.z + p5.z + p6.z) +
									p3.x * (p0.z + p1.z - p6.z - p7.z) - p0.x * p7.z + p5.x * p7.z + p6.x * p7.z) / 16.0 / scaling_factor;
								gradient[x2].z = -(-p5.x * p0.y + p7.x * p0.y + p0.x * p1.y - p5.x * p1.y - p6.x * p1.y - p0.x * p3.y + p6.x * p3.y + p7.x * p3.y +
									p0.x * p5.y - p6.x * p5.y - p7.x * p5.y + p5.x * p6.y - p7.x * p6.y + p1.x * (-p0.y - p3.y + p5.y + p6.y) +
									p3.x * (p0.y + p1.y - p6.y - p7.y) - p0.x * p7.y + p5.x * p7.y + p6.x * p7.y) / 16.0 / scaling_factor;
							}
							if (return_value.second < vertex_minimum_volume[x3]) {
								vertex_minimum_volume[x3] = return_value.second;
								gradient[x3].x = -(p4.y * p0.z + p7.y * p0.z + p0.y * p1.z + p4.y * p1.z - p6.y * p1.z + p0.y * p2.z - p6.y * p2.z - p7.y * p2.z -
									p0.y * p4.z + p6.y * p4.z + p7.y * p4.z - p4.y * p6.z - p7.y * p6.z + p1.y * (-p0.z + p2.z - p4.z + p6.z) -
									p0.y * p7.z - p4.y * p7.z + p6.y * p7.z + p2.y * (-p0.z - p1.z + p6.z + p7.z)) / 16.0 / scaling_factor;
								gradient[x3].y = (p4.x * p0.z + p7.x * p0.z + p0.x * p1.z + p4.x * p1.z - p6.x * p1.z + p0.x * p2.z - p6.x * p2.z - p7.x * p2.z -
									p0.x * p4.z + p6.x * p4.z + p7.x * p4.z - p4.x * p6.z - p7.x * p6.z + p1.x * (-p0.z + p2.z - p4.z + p6.z) -
									p0.x * p7.z - p4.x * p7.z + p6.x * p7.z + p2.x * (-p0.z - p1.z + p6.z + p7.z)) / 16.0 / scaling_factor;
								gradient[x3].z = -(p4.x * p0.y + p7.x * p0.y + p0.x * p1.y + p4.x * p1.y - p6.x * p1.y + p0.x * p2.y - p6.x * p2.y - p7.x * p2.y -
									p0.x * p4.y + p6.x * p4.y + p7.x * p4.y - p4.x * p6.y - p7.x * p6.y + p1.x * (-p0.y + p2.y - p4.y + p6.y) -
									p0.x * p7.y - p4.x * p7.y + p6.x * p7.y + p2.x * (-p0.y - p1.y + p6.y + p7.y)) / 16.0 / scaling_factor;
							}
							if (return_value.second < vertex_minimum_volume[x4]) {
								vertex_minimum_volume[x4] = return_value.second;
								gradient[x4].x = (-p5.y * p0.z + p7.y * p0.z + p0.y * p1.z - p5.y * p1.z - p6.y * p1.z - p0.y * p3.z + p6.y * p3.z + p7.y * p3.z +
									p0.y * p5.z - p6.y * p5.z - p7.y * p5.z + p5.y * p6.z - p7.y * p6.z + p1.y * (-p0.z - p3.z + p5.z + p6.z) +
									p3.y * (p0.z + p1.z - p6.z - p7.z) - p0.y * p7.z + p5.y * p7.z + p6.y * p7.z) / 16.0 / scaling_factor;
								gradient[x4].y = -(-p5.x * p0.z + p7.x * p0.z + p0.x * p1.z - p5.x * p1.z - p6.x * p1.z - p0.x * p3.z + p6.x * p3.z + p7.x * p3.z +
									p0.x * p5.z - p6.x * p5.z - p7.x * p5.z + p5.x * p6.z - p7.x * p6.z + p1.x * (-p0.z - p3.z + p5.z + p6.z) +
									p3.x * (p0.z + p1.z - p6.z - p7.z) - p0.x * p7.z + p5.x * p7.z + p6.x * p7.z) / 16.0 / scaling_factor;
								gradient[x4].z = (-p5.x * p0.y + p7.x * p0.y + p0.x * p1.y - p5.x * p1.y - p6.x * p1.y - p0.x * p3.y + p6.x * p3.y + p7.x * p3.y +
									p0.x * p5.y - p6.x * p5.y - p7.x * p5.y + p5.x * p6.y - p7.x * p6.y + p1.x * (-p0.y - p3.y + p5.y + p6.y) +
									p3.x * (p0.y + p1.y - p6.y - p7.y) - p0.x * p7.y + p5.x * p7.y + p6.x * p7.y) / 16.0 / scaling_factor;
							}
							if (return_value.second < vertex_minimum_volume[x5]) {
								vertex_minimum_volume[x5] = return_value.second;
								gradient[x5].x = (p4.y * p0.z + p7.y * p0.z + p0.y * p1.z + p4.y * p1.z - p6.y * p1.z + p0.y * p2.z - p6.y * p2.z - p7.y * p2.z -
									p0.y * p4.z + p6.y * p4.z + p7.y * p4.z - p4.y * p6.z - p7.y * p6.z + p1.y * (-p0.z + p2.z - p4.z + p6.z) -
									p0.y * p7.z - p4.y * p7.z + p6.y * p7.z + p2.y * (-p0.z - p1.z + p6.z + p7.z)) / 16.0 / scaling_factor;
								gradient[x5].y = -(p4.x * p0.z + p7.x * p0.z + p0.x * p1.z + p4.x * p1.z - p6.x * p1.z + p0.x * p2.z - p6.x * p2.z - p7.x * p2.z -
									p0.x * p4.z + p6.x * p4.z + p7.x * p4.z - p4.x * p6.z - p7.x * p6.z + p1.x * (-p0.z + p2.z - p4.z + p6.z) -
									p0.x * p7.z - p4.x * p7.z + p6.x * p7.z + p2.x * (-p0.z - p1.z + p6.z + p7.z)) / 16.0 / scaling_factor;
								gradient[x5].z = (p4.x * p0.y + p7.x * p0.y + p0.x * p1.y + p4.x * p1.y - p6.x * p1.y + p0.x * p2.y - p6.x * p2.y - p7.x * p2.y -
									p0.x * p4.y + p6.x * p4.y + p7.x * p4.y - p4.x * p6.y - p7.x * p6.y + p1.x * (-p0.y + p2.y - p4.y + p6.y) -
									p0.x * p7.y - p4.x * p7.y + p6.x * p7.y + p2.x * (-p0.y - p1.y + p6.y + p7.y)) / 16.0 / scaling_factor;
							}
							if (return_value.second < vertex_minimum_volume[x6]) {
								vertex_minimum_volume[x6] = return_value.second;
								gradient[x6].x = (p4.y * p1.z + p5.y * p1.z + p1.y * p2.z + p5.y * p2.z - p7.y * p2.z + p1.y * p3.z - p4.y * p3.z - p7.y * p3.z -
									p1.y * p4.z - p5.y * p4.z + p7.y * p4.z - p1.y * p5.z + p4.y * p5.z + p7.y * p5.z - p4.y * p7.z - p5.y * p7.z +
									p3.y * (-p1.z - p2.z + p4.z + p7.z) + p2.y * (-p1.z + p3.z - p5.z + p7.z)) / 16.0 / scaling_factor;
								gradient[x6].y = -(p4.x * p1.z + p5.x * p1.z + p1.x * p2.z + p5.x * p2.z - p7.x * p2.z + p1.x * p3.z - p4.x * p3.z - p7.x * p3.z -
									p1.x * p4.z - p5.x * p4.z + p7.x * p4.z - p1.x * p5.z + p4.x * p5.z + p7.x * p5.z - p4.x * p7.z - p5.x * p7.z +
									p3.x * (-p1.z - p2.z + p4.z + p7.z) + p2.x * (-p1.z + p3.z - p5.z + p7.z)) / 16.0 / scaling_factor;
								gradient[x6].z = (p4.x * p1.y + p5.x * p1.y + p1.x * p2.y + p5.x * p2.y - p7.x * p2.y + p1.x * p3.y - p4.x * p3.y - p7.x * p3.y -
									p1.x * p4.y - p5.x * p4.y + p7.x * p4.y - p1.x * p5.y + p4.x * p5.y + p7.x * p5.y - p4.x * p7.y - p5.x * p7.y +
									p3.x * (-p1.y - p2.y + p4.y + p7.y) + p2.x * (-p1.y + p3.y - p5.y + p7.y)) / 16.0 / scaling_factor;
							}
							if (return_value.second < vertex_minimum_volume[x7]) {
								vertex_minimum_volume[x7] = return_value.second;
								gradient[x7].x = -(p4.y * p0.z + p5.y * p0.z + p0.y * p2.z - p5.y * p2.z - p6.y * p2.z + p0.y * p3.z + p4.y * p3.z - p6.y * p3.z -
									p0.y * p4.z + p5.y * p4.z + p6.y * p4.z - p0.y * p5.z - p4.y * p5.z + p6.y * p5.z - p4.y * p6.z - p5.y * p6.z +
									p3.y * (-p0.z + p2.z - p4.z + p6.z) + p2.y * (-p0.z - p3.z + p5.z + p6.z)) / 16.0 / scaling_factor;
								gradient[x7].y = (p4.x * p0.z + p5.x * p0.z + p0.x * p2.z - p5.x * p2.z - p6.x * p2.z + p0.x * p3.z + p4.x * p3.z - p6.x * p3.z -
									p0.x * p4.z + p5.x * p4.z + p6.x * p4.z - p0.x * p5.z - p4.x * p5.z + p6.x * p5.z - p4.x * p6.z - p5.x * p6.z +
									p3.x * (-p0.z + p2.z - p4.z + p6.z) + p2.x * (-p0.z - p3.z + p5.z + p6.z)) / 16.0 / scaling_factor;
								gradient[x7].z = -(p4.x * p0.y + p5.x * p0.y + p0.x * p2.y - p5.x * p2.y - p6.x * p2.y + p0.x * p3.y + p4.x * p3.y - p6.x * p3.y -
									p0.x * p4.y + p5.x * p4.y + p6.x * p4.y - p0.x * p5.y - p4.x * p5.y + p6.x * p5.y - p4.x * p6.y - p5.x * p6.y +
									p3.x * (-p0.y + p2.y - p4.y + p6.y) + p2.x * (-p0.y - p3.y + p5.y + p6.y)) / 16.0 / scaling_factor;
							}
						}
						// 1
						else if (return_value.first == 1) {
							if (return_value.second < vertex_minimum_volume[x0]) {
								vertex_minimum_volume[x0] = return_value.second;
								gradient[x0].x = (p4.y * (-p1.z + p3.z) + p3.y * (p1.z - p4.z) + p1.y * (-p3.z + p4.z)) / scaling_factor;
								gradient[x0].y = (p4.x * (p1.z - p3.z) + p1.x * (p3.z - p4.z) + p3.x * (-p1.z + p4.z)) / scaling_factor;
								gradient[x0].z = (p4.x * (-p1.y + p3.y) + p3.x * (p1.y - p4.y) + p1.x * (-p3.y + p4.y)) / scaling_factor;
							}
							if (return_value.second < vertex_minimum_volume[x1]) {
								vertex_minimum_volume[x1] = return_value.second;
								gradient[x1].x = (p4.y * (p0.z - p3.z) + p0.y * (p3.z - p4.z) + p3.y * (-p0.z + p4.z)) / scaling_factor;
								gradient[x1].y = (p4.x * (-p0.z + p3.z) + p3.x * (p0.z - p4.z) + p0.x * (-p3.z + p4.z)) / scaling_factor;
								gradient[x1].z = (p4.x * (p0.y - p3.y) + p0.x * (p3.y - p4.y) + p3.x * (-p0.y + p4.y)) / scaling_factor;
							}
							if (return_value.second < vertex_minimum_volume[x3]) {
								vertex_minimum_volume[x3] = return_value.second;
								gradient[x3].x = (p4.y * (-p0.z + p1.z) + p1.y * (p0.z - p4.z) + p0.y * (-p1.z + p4.z)) / scaling_factor;
								gradient[x3].y = (p4.x * (p0.z - p1.z) + p0.x * (p1.z - p4.z) + p1.x * (-p0.z + p4.z)) / scaling_factor;
								gradient[x3].z = (p4.x * (-p0.y + p1.y) + p1.x * (p0.y - p4.y) + p0.x * (-p1.y + p4.y)) / scaling_factor;
							}
							if (return_value.second < vertex_minimum_volume[x4]) {
								vertex_minimum_volume[x4] = return_value.second;
								gradient[x4].x = (p3.y * (p0.z - p1.z) + p0.y * (p1.z - p3.z) + p1.y * (-p0.z + p3.z)) / scaling_factor;
								gradient[x4].y = (p3.x * (-p0.z + p1.z) + p1.x * (p0.z - p3.z) + p0.x * (-p1.z + p3.z)) / scaling_factor;
								gradient[x4].z = (p3.x * (p0.y - p1.y) + p0.x * (p1.y - p3.y) + p1.x * (-p0.y + p3.y)) / scaling_factor;
							}
						}
						// 2
						else if (return_value.first == 2) {
							if (return_value.second < vertex_minimum_volume[x0]) {
								vertex_minimum_volume[x0] = return_value.second;
								gradient[x0].x = (p5.y * (-p1.z + p2.z) + p2.y * (p1.z - p5.z) + p1.y * (-p2.z + p5.z)) / scaling_factor;
								gradient[x0].y = (p5.x * (p1.z - p2.z) + p1.x * (p2.z - p5.z) + p2.x * (-p1.z + p5.z)) / scaling_factor;
								gradient[x0].z = (p5.x * (-p1.y + p2.y) + p2.x * (p1.y - p5.y) + p1.x * (-p2.y + p5.y)) / scaling_factor;
							}
							if (return_value.second < vertex_minimum_volume[x1]) {
								vertex_minimum_volume[x1] = return_value.second;
								gradient[x1].x = (p5.y * (p0.z - p2.z) + p0.y * (p2.z - p5.z) + p2.y * (-p0.z + p5.z)) / scaling_factor;
								gradient[x1].y = (p5.x * (-p0.z + p2.z) + p2.x * (p0.z - p5.z) + p0.x * (-p2.z + p5.z)) / scaling_factor;
								gradient[x1].z = (p5.x * (p0.y - p2.y) + p0.x * (p2.y - p5.y) + p2.x * (-p0.y + p5.y)) / scaling_factor;
							}
							if (return_value.second < vertex_minimum_volume[x2]) {
								vertex_minimum_volume[x2] = return_value.second;
								gradient[x2].x = (p5.y * (-p0.z + p1.z) + p1.y * (p0.z - p5.z) + p0.y * (-p1.z + p5.z)) / scaling_factor;
								gradient[x2].y = (p5.x * (p0.z - p1.z) + p0.x * (p1.z - p5.z) + p1.x * (-p0.z + p5.z)) / scaling_factor;
								gradient[x2].z = (p5.x * (-p0.y + p1.y) + p1.x * (p0.y - p5.y) + p0.x * (-p1.y + p5.y)) / scaling_factor;
							}
							if (return_value.second < vertex_minimum_volume[x5]) {
								vertex_minimum_volume[x5] = return_value.second;
								gradient[x5].x = (p2.y * (p0.z - p1.z) + p0.y * (p1.z - p2.z) + p1.y * (-p0.z + p2.z)) / scaling_factor;
								gradient[x5].y = (p2.x * (-p0.z + p1.z) + p1.x * (p0.z - p2.z) + p0.x * (-p1.z + p2.z)) / scaling_factor;
								gradient[x5].z = (p2.x * (p0.y - p1.y) + p0.x * (p1.y - p2.y) + p1.x * (-p0.y + p2.y)) / scaling_factor;
							}
						}
						// 3
						else if (return_value.first == 3) {
							if (return_value.second < vertex_minimum_volume[x1]) {
								vertex_minimum_volume[x1] = return_value.second;
								gradient[x1].x = (p6.y * (-p2.z + p3.z) + p3.y * (p2.z - p6.z) + p2.y * (-p3.z + p6.z)) / scaling_factor;
								gradient[x1].y = (p6.x * (p2.z - p3.z) + p2.x * (p3.z - p6.z) + p3.x * (-p2.z + p6.z)) / scaling_factor;
								gradient[x1].z = (p6.x * (-p2.y + p3.y) + p3.x * (p2.y - p6.y) + p2.x * (-p3.y + p6.y)) / scaling_factor;
							}
							if (return_value.second < vertex_minimum_volume[x2]) {
								vertex_minimum_volume[x2] = return_value.second;
								gradient[x2].x = (p6.y * (p1.z - p3.z) + p1.y * (p3.z - p6.z) + p3.y * (-p1.z + p6.z)) / scaling_factor;
								gradient[x2].y = (p6.x * (-p1.z + p3.z) + p3.x * (p1.z - p6.z) + p1.x * (-p3.z + p6.z)) / scaling_factor;
								gradient[x2].z = (p6.x * (p1.y - p3.y) + p1.x * (p3.y - p6.y) + p3.x * (-p1.y + p6.y)) / scaling_factor;
							}
							if (return_value.second < vertex_minimum_volume[x3]) {
								vertex_minimum_volume[x3] = return_value.second;
								gradient[x3].x = (p6.y * (-p1.z + p2.z) + p2.y * (p1.z - p6.z) + p1.y * (-p2.z + p6.z)) / scaling_factor;
								gradient[x3].y = (p6.x * (p1.z - p2.z) + p1.x * (p2.z - p6.z) + p2.x * (-p1.z + p6.z)) / scaling_factor;
								gradient[x3].z = (p6.x * (-p1.y + p2.y) + p2.x * (p1.y - p6.y) + p1.x * (-p2.y + p6.y)) / scaling_factor;
							}
							if (return_value.second < vertex_minimum_volume[x6]) {
								vertex_minimum_volume[x6] = return_value.second;
								gradient[x6].x = (p3.y * (p1.z - p2.z) + p1.y * (p2.z - p3.z) + p2.y * (-p1.z + p3.z)) / scaling_factor;
								gradient[x6].y = (p3.x * (-p1.z + p2.z) + p2.x * (p1.z - p3.z) + p1.x * (-p2.z + p3.z)) / scaling_factor;
								gradient[x6].z = (p3.x * (p1.y - p2.y) + p1.x * (p2.y - p3.y) + p2.x * (-p1.y + p3.y)) / scaling_factor;
							}
						}
						// 4
						else if (return_value.first == 4) {
							if (return_value.second < vertex_minimum_volume[x0]) {
								vertex_minimum_volume[x0] = return_value.second;
								gradient[x0].x = (p7.y * (-p2.z + p3.z) + p3.y * (p2.z - p7.z) + p2.y * (-p3.z + p7.z)) / scaling_factor;
								gradient[x0].y = (p7.x * (p2.z - p3.z) + p2.x * (p3.z - p7.z) + p3.x * (-p2.z + p7.z)) / scaling_factor;
								gradient[x0].z = (p7.x * (-p2.y + p3.y) + p3.x * (p2.y - p7.y) + p2.x * (-p3.y + p7.y)) / scaling_factor;
							}
							if (return_value.second < vertex_minimum_volume[x2]) {
								vertex_minimum_volume[x2] = return_value.second;
								gradient[x2].x = (p7.y * (p0.z - p3.z) + p0.y * (p3.z - p7.z) + p3.y * (-p0.z + p7.z)) / scaling_factor;
								gradient[x2].y = (p7.x * (-p0.z + p3.z) + p3.x * (p0.z - p7.z) + p0.x * (-p3.z + p7.z)) / scaling_factor;
								gradient[x2].z = (p7.x * (p0.y - p3.y) + p0.x * (p3.y - p7.y) + p3.x * (-p0.y + p7.y)) / scaling_factor;
							}
							if (return_value.second < vertex_minimum_volume[x3]) {
								vertex_minimum_volume[x3] = return_value.second;
								gradient[x3].x = (p7.y * (-p0.z + p2.z) + p2.y * (p0.z - p7.z) + p0.y * (-p2.z + p7.z)) / scaling_factor;
								gradient[x3].y = (p7.x * (p0.z - p2.z) + p0.x * (p2.z - p7.z) + p2.x * (-p0.z + p7.z)) / scaling_factor;
								gradient[x3].z = (p7.x * (-p0.y + p2.y) + p2.x * (p0.y - p7.y) + p0.x * (-p2.y + p7.y)) / scaling_factor;
							}
							if (return_value.second < vertex_minimum_volume[x7]) {
								vertex_minimum_volume[x7] = return_value.second;
								gradient[x7].x = (p3.y * (p0.z - p2.z) + p0.y * (p2.z - p3.z) + p2.y * (-p0.z + p3.z)) / scaling_factor;
								gradient[x7].y = (p3.x * (-p0.z + p2.z) + p2.x * (p0.z - p3.z) + p0.x * (-p2.z + p3.z)) / scaling_factor;
								gradient[x7].z = (p3.x * (p0.y - p2.y) + p0.x * (p2.y - p3.y) + p2.x * (-p0.y + p3.y)) / scaling_factor;
							}
						}
						// 5
						else if (return_value.first == 5) {
							if (return_value.second < vertex_minimum_volume[x0]) {
								vertex_minimum_volume[x0] = return_value.second;
								gradient[x0].x = (p7.y * (-p4.z + p5.z) + p5.y * (p4.z - p7.z) + p4.y * (-p5.z + p7.z)) / scaling_factor;
								gradient[x0].y = (p7.x * (p4.z - p5.z) + p4.x * (p5.z - p7.z) + p5.x * (-p4.z + p7.z)) / scaling_factor;
								gradient[x0].z = (p7.x * (-p4.y + p5.y) + p5.x * (p4.y - p7.y) + p4.x * (-p5.y + p7.y)) / scaling_factor;
							}
							if (return_value.second < vertex_minimum_volume[x4]) {
								vertex_minimum_volume[x4] = return_value.second;
								gradient[x4].x = (p7.y * (p0.z - p5.z) + p0.y * (p5.z - p7.z) + p5.y * (-p0.z + p7.z)) / scaling_factor;
								gradient[x4].y = (p7.x * (-p0.z + p5.z) + p5.x * (p0.z - p7.z) + p0.x * (-p5.z + p7.z)) / scaling_factor;
								gradient[x4].z = (p7.x * (p0.y - p5.y) + p0.x * (p5.y - p7.y) + p5.x * (-p0.y + p7.y)) / scaling_factor;
							}
							if (return_value.second < vertex_minimum_volume[x5]) {
								vertex_minimum_volume[x5] = return_value.second;
								gradient[x5].x = (p7.y * (-p0.z + p4.z) + p4.y * (p0.z - p7.z) + p0.y * (-p4.z + p7.z)) / scaling_factor;
								gradient[x5].y = (p7.x * (p0.z - p4.z) + p0.x * (p4.z - p7.z) + p4.x * (-p0.z + p7.z)) / scaling_factor;
								gradient[x5].z = (p7.x * (-p0.y + p4.y) + p4.x * (p0.y - p7.y) + p0.x * (-p4.y + p7.y)) / scaling_factor;
							}
							if (return_value.second < vertex_minimum_volume[x7]) {
								vertex_minimum_volume[x7] = return_value.second;
								gradient[x7].x = (p5.y * (p0.z - p4.z) + p0.y * (p4.z - p5.z) + p4.y * (-p0.z + p5.z)) / scaling_factor;
								gradient[x7].y = (p5.x * (-p0.z + p4.z) + p4.x * (p0.z - p5.z) + p0.x * (-p4.z + p5.z)) / scaling_factor;
								gradient[x7].z = (p5.x * (p0.y - p4.y) + p0.x * (p4.y - p5.y) + p4.x * (-p0.y + p5.y)) / scaling_factor;
							}
						}
						// 6
						else if (return_value.first == 6) {
							if (return_value.second < vertex_minimum_volume[x1]) {
								vertex_minimum_volume[x1] = return_value.second;
								gradient[x1].x = (p6.y * (-p4.z + p5.z) + p5.y * (p4.z - p6.z) + p4.y * (-p5.z + p6.z)) / scaling_factor;
								gradient[x1].y = (p6.x * (p4.z - p5.z) + p4.x * (p5.z - p6.z) + p5.x * (-p4.z + p6.z)) / scaling_factor;
								gradient[x1].z = (p6.x * (-p4.y + p5.y) + p5.x * (p4.y - p6.y) + p4.x * (-p5.y + p6.y)) / scaling_factor;
							}
							if (return_value.second < vertex_minimum_volume[x4]) {
								vertex_minimum_volume[x4] = return_value.second;
								gradient[x4].x = (p6.y * (p1.z - p5.z) + p1.y * (p5.z - p6.z) + p5.y * (-p1.z + p6.z)) / scaling_factor;
								gradient[x4].y = (p6.x * (-p1.z + p5.z) + p5.x * (p1.z - p6.z) + p1.x * (-p5.z + p6.z)) / scaling_factor;
								gradient[x4].z = (p6.x * (p1.y - p5.y) + p1.x * (p5.y - p6.y) + p5.x * (-p1.y + p6.y)) / scaling_factor;
							}
							if (return_value.second < vertex_minimum_volume[x5]) {
								vertex_minimum_volume[x5] = return_value.second;
								gradient[x5].x = (p6.y * (-p1.z + p4.z) + p4.y * (p1.z - p6.z) + p1.y * (-p4.z + p6.z)) / scaling_factor;
								gradient[x5].y = (p6.x * (p1.z - p4.z) + p1.x * (p4.z - p6.z) + p4.x * (-p1.z + p6.z)) / scaling_factor;
								gradient[x5].z = (p6.x * (-p1.y + p4.y) + p4.x * (p1.y - p6.y) + p1.x * (-p4.y + p6.y)) / scaling_factor;
							}
							if (return_value.second < vertex_minimum_volume[x6]) {
								vertex_minimum_volume[x6] = return_value.second;
								gradient[x6].x = (p5.y * (p1.z - p4.z) + p1.y * (p4.z - p5.z) + p4.y * (-p1.z + p5.z)) / scaling_factor;
								gradient[x6].y = (p5.x * (-p1.z + p4.z) + p4.x * (p1.z - p5.z) + p1.x * (-p4.z + p5.z)) / scaling_factor;
								gradient[x6].z = (p5.x * (p1.y - p4.y) + p1.x * (p4.y - p5.y) + p4.x * (-p1.y + p5.y)) / scaling_factor;
							}
						}
						// 7
						else if (return_value.first == 7) {
							if (return_value.second < vertex_minimum_volume[x2]) {
								vertex_minimum_volume[x2] = return_value.second;
								gradient[x2].x = (p7.y * (-p5.z + p6.z) + p6.y * (p5.z - p7.z) + p5.y * (-p6.z + p7.z)) / scaling_factor;
								gradient[x2].y = (p7.x * (p5.z - p6.z) + p5.x * (p6.z - p7.z) + p6.x * (-p5.z + p7.z)) / scaling_factor;
								gradient[x2].z = (p7.x * (-p5.y + p6.y) + p6.x * (p5.y - p7.y) + p5.x * (-p6.y + p7.y)) / scaling_factor;
							}
							if (return_value.second < vertex_minimum_volume[x5]) {
								vertex_minimum_volume[x5] = return_value.second;
								gradient[x5].x = (p7.y * (p2.z - p6.z) + p2.y * (p6.z - p7.z) + p6.y * (-p2.z + p7.z)) / scaling_factor;
								gradient[x5].y = (p7.x * (-p2.z + p6.z) + p6.x * (p2.z - p7.z) + p2.x * (-p6.z + p7.z)) / scaling_factor;
								gradient[x5].z = (p7.x * (p2.y - p6.y) + p2.x * (p6.y - p7.y) + p6.x * (-p2.y + p7.y)) / scaling_factor;
							}
							if (return_value.second < vertex_minimum_volume[x6]) {
								vertex_minimum_volume[x6] = return_value.second;
								gradient[x6].x = (p7.y * (-p2.z + p5.z) + p5.y * (p2.z - p7.z) + p2.y * (-p5.z + p7.z)) / scaling_factor;
								gradient[x6].y = (p7.x * (p2.z - p5.z) + p2.x * (p5.z - p7.z) + p5.x * (-p2.z + p7.z)) / scaling_factor;
								gradient[x6].z = (p7.x * (-p2.y + p5.y) + p5.x * (p2.y - p7.y) + p2.x * (-p5.y + p7.y)) / scaling_factor;
							}
							if (return_value.second < vertex_minimum_volume[x7]) {
								vertex_minimum_volume[x7] = return_value.second;
								gradient[x7].x = (p6.y * (p2.z - p5.z) + p2.y * (p5.z - p6.z) + p5.y * (-p2.z + p6.z)) / scaling_factor;
								gradient[x7].y = (p6.x * (-p2.z + p5.z) + p5.x * (p2.z - p6.z) + p2.x * (-p5.z + p6.z)) / scaling_factor;
								gradient[x7].z = (p6.x * (p2.y - p5.y) + p2.x * (p5.y - p6.y) + p5.x * (-p2.y + p6.y)) / scaling_factor;
							}
						}
						// 8
						else {
							if (return_value.second < vertex_minimum_volume[x3]) {
								vertex_minimum_volume[x3] = return_value.second;
								gradient[x3].x = (p7.y * (-p4.z + p6.z) + p6.y * (p4.z - p7.z) + p4.y * (-p6.z + p7.z)) / scaling_factor;
								gradient[x3].y = (p7.x * (p4.z - p6.z) + p4.x * (p6.z - p7.z) + p6.x * (-p4.z + p7.z)) / scaling_factor;
								gradient[x3].z = (p7.x * (-p4.y + p6.y) + p6.x * (p4.y - p7.y) + p4.x * (-p6.y + p7.y)) / scaling_factor;
							}
							if (return_value.second < vertex_minimum_volume[x4]) {
								vertex_minimum_volume[x4] = return_value.second;
								gradient[x4].x = (p7.y * (p3.z - p6.z) + p3.y * (p6.z - p7.z) + p6.y * (-p3.z + p7.z)) / scaling_factor;
								gradient[x4].y = (p7.x * (-p3.z + p6.z) + p6.x * (p3.z - p7.z) + p3.x * (-p6.z + p7.z)) / scaling_factor;
								gradient[x4].z = (p7.x * (p3.y - p6.y) + p3.x * (p6.y - p7.y) + p6.x * (-p3.y + p7.y)) / scaling_factor;
							}
							if (return_value.second < vertex_minimum_volume[x6]) {
								vertex_minimum_volume[x6] = return_value.second;
								gradient[x6].x = (p7.y * (-p3.z + p4.z) + p4.y * (p3.z - p7.z) + p3.y * (-p4.z + p7.z)) / scaling_factor;
								gradient[x6].y = (p7.x * (p3.z - p4.z) + p3.x * (p4.z - p7.z) + p4.x * (-p3.z + p7.z)) / scaling_factor;
								gradient[x6].z = (p7.x * (-p3.y + p4.y) + p4.x * (p3.y - p7.y) + p3.x * (-p4.y + p7.y)) / scaling_factor;
							}
							if (return_value.second < vertex_minimum_volume[x7]) {
								vertex_minimum_volume[x7] = return_value.second;
								gradient[x7].x = (p6.y * (p3.z - p4.z) + p3.y * (p4.z - p6.z) + p4.y * (-p3.z + p6.z)) / scaling_factor;
								gradient[x7].y = (p6.x * (-p3.z + p4.z) + p4.x * (p3.z - p6.z) + p3.x * (-p4.z + p6.z)) / scaling_factor;
								gradient[x7].z = (p6.x * (p3.y - p4.y) + p3.x * (p4.y - p6.y) + p4.x * (-p3.y + p6.y)) / scaling_factor;
							}
						}
					}
					for (long long vertex_index = 0; vertex_index < original_vertex_number; ++vertex_index) {
						MeshTypes::Vec3 original_coordinate = optimized_hexahedra_mesh_.vertices[vertex_index];
						optimized_hexahedra_mesh_.vertices[vertex_index] += 250 * learning_rate * gradient[vertex_index];
						for (const auto& adj_elem : elements_around_vertex[vertex_index]) {
							const std::pair<size_t, double> return_value = MinMeshQuality(
								optimized_hexahedra_mesh_.vertices[optimized_hexahedra_mesh_.faces[adj_elem][0]],
								optimized_hexahedra_mesh_.vertices[optimized_hexahedra_mesh_.faces[adj_elem][1]],
								optimized_hexahedra_mesh_.vertices[optimized_hexahedra_mesh_.faces[adj_elem][2]],
								optimized_hexahedra_mesh_.vertices[optimized_hexahedra_mesh_.faces[adj_elem][3]],
								optimized_hexahedra_mesh_.vertices[optimized_hexahedra_mesh_.faces[adj_elem][4]],
								optimized_hexahedra_mesh_.vertices[optimized_hexahedra_mesh_.faces[adj_elem][5]],
								optimized_hexahedra_mesh_.vertices[optimized_hexahedra_mesh_.faces[adj_elem][6]],
								optimized_hexahedra_mesh_.vertices[optimized_hexahedra_mesh_.faces[adj_elem][7]]);
							if (return_value.second <= 0) {
								optimized_hexahedra_mesh_.vertices[vertex_index] = original_coordinate;
								break;
							}
						}
					}
					std::fill(gradient.begin(), gradient.end(), MeshTypes::Vec3());

#pragma omp parallel for
					for (long long vertex_index = 0; vertex_index < original_vertex_number; ++vertex_index) {
						double min_edge_length = std::numeric_limits<double>::max();
						size_t min_edge_vertex_index;
						for (const auto& neighbors : vertex_neighbors[vertex_index]) {
							MeshTypes::Vec3 edge = optimized_hexahedra_mesh_.vertices[vertex_index] -
								optimized_hexahedra_mesh_.vertices[neighbors];
							double length = edge.norm();
							if (length < min_edge_length) {
								min_edge_vertex_index = neighbors;
								min_edge_length = length;
							}
						}
						MeshTypes::Vec3 direction = optimized_hexahedra_mesh_.vertices[vertex_index] -
							optimized_hexahedra_mesh_.vertices[min_edge_vertex_index];

						MeshTypes::Vec3 weighted_direction;
						for (const auto& neighbors : vertex_neighbors[vertex_index]) {
							weighted_direction += optimized_hexahedra_mesh_.vertices[neighbors];
						}
						weighted_direction /= vertex_neighbors[vertex_index].size();
						weighted_direction -= optimized_hexahedra_mesh_.vertices[vertex_index];
						weighted_direction -= direction *
							weighted_direction.dot(direction) / direction.dot(direction);

						MeshTypes::Vec3 original_coordinate = optimized_hexahedra_mesh_.vertices[vertex_index];
						optimized_hexahedra_mesh_.vertices[vertex_index] += 500 * learning_rate * weighted_direction;
						for (const auto& adj_elem : elements_around_vertex[vertex_index]) {
							const std::pair<size_t, double> return_value = MinMeshQuality(
								optimized_hexahedra_mesh_.vertices[optimized_hexahedra_mesh_.faces[adj_elem][0]],
								optimized_hexahedra_mesh_.vertices[optimized_hexahedra_mesh_.faces[adj_elem][1]],
								optimized_hexahedra_mesh_.vertices[optimized_hexahedra_mesh_.faces[adj_elem][2]],
								optimized_hexahedra_mesh_.vertices[optimized_hexahedra_mesh_.faces[adj_elem][3]],
								optimized_hexahedra_mesh_.vertices[optimized_hexahedra_mesh_.faces[adj_elem][4]],
								optimized_hexahedra_mesh_.vertices[optimized_hexahedra_mesh_.faces[adj_elem][5]],
								optimized_hexahedra_mesh_.vertices[optimized_hexahedra_mesh_.faces[adj_elem][6]],
								optimized_hexahedra_mesh_.vertices[optimized_hexahedra_mesh_.faces[adj_elem][7]]);
							if (return_value.second <= 0) {
								optimized_hexahedra_mesh_.vertices[vertex_index] = original_coordinate;
								break;
							}
						}
					}

					// project surface vertices onto triangle surface
#pragma omp parallel for
					for (long long vertex_index = original_vertex_number; vertex_index < vertex_number; ++vertex_index) {
						// end vertex
						if (where_should_end_vertex_go[vertex_index] > -1) {
							MeshTypes::Vec3 vertex_gradient =
								triangle_mesh_.vertices[where_should_end_vertex_go[vertex_index]] -
								optimized_hexahedra_mesh_.vertices[vertex_index];
							if (vertex_gradient.norm() >= max_distance_threshold) {
								gradient[vertex_index] = lambda_projection * vertex_gradient;
							}
							continue;
						}
						// middle vertex
						else if (where_should_end_vertex_go[vertex_index] == -1) {
							MeshTypes::Vec3 vertex_gradient =
								0.5 * optimized_hexahedra_mesh_.vertices[where_should_middle_vertex_go[vertex_index].first] +
								0.5 * optimized_hexahedra_mesh_.vertices[where_should_middle_vertex_go[vertex_index].second] -
								optimized_hexahedra_mesh_.vertices[vertex_index];
							if (vertex_gradient.norm() >= max_distance_threshold) {
								gradient[vertex_index] = lambda_projection * vertex_gradient;
							}
							continue;
						}

						if ((optimized_hexahedra_mesh_.vertices[vertex_index] -
							closest_projection_vertex[vertex_index]).norm() < max_distance_threshold) {
							continue;
						}
						std::uniform_real_distribution<double> y_distribution(triangle_low_bound.y, triangle_up_bound.y);
						std::uniform_real_distribution<double> z_distribution(triangle_low_bound.z, triangle_up_bound.z);
						std::default_random_engine random_generator(rand());
						int intersect_count = -2;
						while (intersect_count == -2) {
							intersect_count = 0;
							MeshTypes::Vec3 ray = MeshTypes::Vec3(triangle_low_bound.x,
								y_distribution(random_generator),
								z_distribution(random_generator)) - optimized_hexahedra_mesh_.vertices[vertex_index];
							for (const auto& face : triangle_mesh_.faces) {
								size_t intersect = RayTriangleIntersect(
									optimized_hexahedra_mesh_.vertices[vertex_index],
									ray,
									triangle_mesh_.vertices[face[0]],
									triangle_mesh_.vertices[face[1]],
									triangle_mesh_.vertices[face[2]]);
								if (intersect == 1) {
									++intersect_count;
								}
								else if (intersect == 2) {
									intersect_count = -2;
									break;
								}
							}
						}
						bool inside = (intersect_count % 2 == 1);

						double min_metric = std::numeric_limits<double>::max();

						for (size_t triangle_index = 0; triangle_index < triangle_mesh_.faces.size(); ++triangle_index) {
							MeshTypes::Vec3 vertex_gradient = VertexTriangleDistance(
								optimized_hexahedra_mesh_.vertices[vertex_index],
								triangle_mesh_.vertices[triangle_mesh_.faces[triangle_index][0]],
								triangle_mesh_.vertices[triangle_mesh_.faces[triangle_index][1]],
								triangle_mesh_.vertices[triangle_mesh_.faces[triangle_index][2]]);
							double distance = vertex_gradient.norm();
							const auto& normals = surface_vertex_normals[vertex_index];
							double metric = 1.1 * normals.size() * distance;
							MeshTypes::Vec3 current_normal;
							for (const auto& normal : normals) {
								current_normal = ((optimized_hexahedra_mesh_.vertices[normal[0]] -
									optimized_hexahedra_mesh_.vertices[vertex_index]).cross(
										optimized_hexahedra_mesh_.vertices[normal[1]] -
										optimized_hexahedra_mesh_.vertices[vertex_index])).normalized();
								if (current_normal.norm() < 0.5) {
									current_normal = ((optimized_hexahedra_mesh_.vertices[normal[2]] -
										optimized_hexahedra_mesh_.vertices[vertex_index]).cross(
											optimized_hexahedra_mesh_.vertices[normal[1]] -
											optimized_hexahedra_mesh_.vertices[vertex_index])).normalized();
								}
								if (inside) {
									metric -= current_normal.dot(vertex_gradient);
								}
								else {
									metric += current_normal.dot(vertex_gradient);
								}
							}
							if (metric < min_metric) {
								min_metric = metric;
								gradient[vertex_index] = lambda_projection * vertex_gradient;
								closest_projection_vertex[vertex_index] = vertex_gradient +
									optimized_hexahedra_mesh_.vertices[vertex_index];
							}
						}
					}

					double max_distance = 0;
					size_t max_distance_vertex = 0;
					for (long long vertex_index = original_vertex_number; vertex_index < vertex_number; ++vertex_index) {
						double distance = gradient[vertex_index].norm() / lambda_projection;
						if (max_distance < distance) {
							max_distance = distance;
							max_distance_vertex = vertex_index;
						}
					}

					if (learning_rate > 0.1 * learning_rate_projection * max_distance) {
						learning_rate = 0.1 * learning_rate_projection * max_distance;
						lambda_projection = 10 / max_distance;
					}
					MeshOutput::WriteHexahedraMesh("current_mesh.vtk", optimized_hexahedra_mesh_);

					if (next_smooth_threshold > max_distance && smooth_count < 50) {
						std::cout << std::endl << "Smooth";
						next_smooth_threshold = max_distance;
						++smooth_count;

						// surface smoothing
#pragma omp parallel for
						for (long long vertex_index = 0; vertex_index < vertex_number; ++vertex_index) {
							MeshTypes::Vec3 weighted_gradient;
							for (const auto& neighbors : vertex_neighbors[vertex_index]) {
								weighted_gradient += optimized_hexahedra_mesh_.vertices[neighbors];
							}
							weighted_gradient /= vertex_neighbors[vertex_index].size();
							weighted_gradient -= optimized_hexahedra_mesh_.vertices[vertex_index];

							/*const auto& normals = surface_vertex_normals[vertex_index];
							MeshTypes::Vec3 average_normal;
							for (const auto& normal : normals) {
								average_normal +=
									(optimized_hexahedra_mesh_.vertices[normal[0]] -
										optimized_hexahedra_mesh_.vertices[vertex_index]).cross(
											optimized_hexahedra_mesh_.vertices[normal[1]] -
											optimized_hexahedra_mesh_.vertices[vertex_index]);
							}
							weighted_gradient -= weighted_gradient.dot(average_normal) * average_normal / average_normal.dot(average_normal);*/
							optimized_hexahedra_mesh_.vertices[vertex_index] += 0.1 * weighted_gradient;
						}
					}

					if (max_distance < max_distance_threshold) {
						std::cout << "Converged" << std::endl;
						break;
					}
					std::cout << std::endl << "Max distance: " << max_distance
						<< " at vertex: " << max_distance_vertex << std::endl;
				}
				else {
					std::cout << below_threshold_count << " ";
				}

				// gradient descent
#pragma omp parallel for
				for (long long vertex_index = 0; vertex_index < original_vertex_number; ++vertex_index) {
					optimized_hexahedra_mesh_.vertices[vertex_index] += internal_learning_rate * gradient[vertex_index];
				}
#pragma omp parallel for
				for (long long vertex_index = original_vertex_number; vertex_index < vertex_number; ++vertex_index) {
					optimized_hexahedra_mesh_.vertices[vertex_index] += learning_rate * gradient[vertex_index];
				}
			}

			return optimized_hexahedra_mesh_;
		}
	private:
		const MeshTypes::Mesh3D& triangle_mesh_;
		const MeshTypes::SharpFeature& sharp_feature_;
		const MeshTypes::Mesh3D& hexahedra_mesh_;
	};
}

#endif