#ifndef MESH_PROJECTION_H
#define MESH_PROJECTION_H

#include "Mesh_Embedding.h"
#include "Mesh_Shortest_Path.h"
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <numeric>
#include <iostream>
#include <fstream>

namespace MeshProjection {
    class MeshProjector {
    public:
        MeshProjector(const MeshTypes::Mesh3D& triangle_mesh,
            const MeshTypes::Mesh3D& polygon_mesh,
            const MeshTypes::SharpFeature& sharp_feature)
            : triangle_mesh_(triangle_mesh), polygon_mesh_(polygon_mesh),
            sharp_edges_(sharp_feature.edges), sharp_vertices_(sharp_feature.vertices) {
            triangles_.resize(triangle_mesh_.faces.size());
#pragma omp parallel for
            for (long long face_index = 0; face_index < triangle_mesh_.faces.size(); ++face_index) {
                const auto& face = triangle_mesh_.faces[face_index];
                const auto& vertexA = triangle_mesh_.vertices[face[0]];
                const auto& vertexB = triangle_mesh_.vertices[face[1]];
                const auto& vertexC = triangle_mesh_.vertices[face[2]];
                triangles_[face_index] = { vertexA, vertexB - vertexA, vertexC - vertexA, vertexC - vertexB };
            }

            const size_t vertex_number = polygon_mesh_.vertices.size();
            is_vertex_active_.resize(vertex_number, true);
            auto FindClosestPolygonVertex = [&](size_t triangle_vert) -> size_t {
                double min_distance_squared = INF;
                size_t closest_polygon_vert = 0;
                for (size_t i = 0; i < vertex_number; ++i) {
                    const auto diff = polygon_mesh_.vertices[i] - triangle_mesh_.vertices[triangle_vert];
                    const auto d2 = diff.dot(diff);
                    if (d2 < min_distance_squared) {
                        min_distance_squared = d2;
                        closest_polygon_vert = i;
                    }
                }
                return closest_polygon_vert;
                };
            for (const auto& [start, end] : sharp_edges_) {
                if (!polygon_sharp_vertex_indices.count(start)) {
                    const auto closest_polygon_vert = FindClosestPolygonVertex(start);
                    polygon_sharp_vertex_indices[start] = closest_polygon_vert;
                    is_vertex_active_[closest_polygon_vert] = false;
                }
                if (!polygon_sharp_vertex_indices.count(end)) {
                    const auto closest_polygon_vert = FindClosestPolygonVertex(end);
                    polygon_sharp_vertex_indices[end] = closest_polygon_vert;
                    is_vertex_active_[closest_polygon_vert] = false;
                }
            }
            for (const auto& vert : sharp_vertices_) {
                if (!polygon_sharp_vertex_indices.count(vert)) {
                    const auto closest_polygon_vert = FindClosestPolygonVertex(vert);
                    polygon_sharp_vertex_indices[vert] = closest_polygon_vert;
                    is_vertex_active_[closest_polygon_vert] = false;
                }
            }

            std::vector<bool> valid_vertices(vertex_number, true);
            MeshShortestPath::DijkstraSolver dijkstra(polygon_mesh_, valid_vertices);
            for (const auto& [start, end] : sharp_edges_) {
                const auto& path_indices = dijkstra.ComputePath(polygon_sharp_vertex_indices[start],
                    polygon_sharp_vertex_indices[end],
                    true);
                for (const auto& idx : path_indices) {
                    is_vertex_active_[idx] = false;
                }
                polygon_sharp_edge_vertices_.push_back(path_indices);
            }
        }

        MeshTypes::Mesh3D ComputeProjection() {
            const size_t vertex_number = polygon_mesh_.vertices.size();
            std::vector<double> projection_distances(vertex_number, INF);
            MeshTypes::Mesh3D projection_result;
            projection_result.vertices.resize(vertex_number);
            std::vector<bool> is_vertex_active_aux = is_vertex_active_;

#pragma omp parallel for
            for (long long i = 0; i < sharp_edges_.size(); ++i) {
                const auto start = sharp_edges_[i].first;
                const auto end = sharp_edges_[i].second;
                projection_result.vertices[polygon_sharp_vertex_indices[start]] = triangle_mesh_.vertices[start];
                MeshTypes::Vec3 vec = triangle_mesh_.vertices[start] - polygon_mesh_.vertices[polygon_sharp_vertex_indices[start]];
                projection_distances[polygon_sharp_vertex_indices[start]] = vec.dot(vec);
                double total_length = 0;
                std::vector<double> accumulate_length;
                for (size_t j = 0; j < polygon_sharp_edge_vertices_[i].size() - 1; ++j) {
                    vec = polygon_mesh_.vertices[polygon_sharp_edge_vertices_[i][j + 1]] - polygon_mesh_.vertices[polygon_sharp_edge_vertices_[i][j]];
                    total_length += vec.norm();
                    accumulate_length.push_back(total_length);
                }
                double portion;
                size_t idx;
                for (size_t j = 0; j < polygon_sharp_edge_vertices_[i].size() - 1; ++j) {
                    idx = polygon_sharp_edge_vertices_[i][j + 1];
                    portion = accumulate_length[j] / total_length;
                    projection_result.vertices[idx] = portion * triangle_mesh_.vertices[end] + (1 - portion) * triangle_mesh_.vertices[start];
                    vec = projection_result.vertices[idx] - polygon_mesh_.vertices[idx];
                    projection_distances[idx] = vec.dot(vec);
                }
            }
#pragma omp parallel for
            for (long long i = 0; i < sharp_vertices_.size(); ++i) {
                const auto idx = sharp_vertices_[i];
                projection_result.vertices[polygon_sharp_vertex_indices[idx]] = triangle_mesh_.vertices[idx];
                MeshTypes::Vec3 vec = triangle_mesh_.vertices[idx] - polygon_mesh_.vertices[polygon_sharp_vertex_indices[idx]];
                projection_distances[polygon_sharp_vertex_indices[idx]] = vec.dot(vec);
            }

#pragma omp parallel for
            for (long long i = 0; i < vertex_number; ++i) {
                if (!is_vertex_active_aux[i]) {
                    continue;
                }
                const auto& point = polygon_mesh_.vertices[i];
                MeshTypes::Vec3 closest_point;
                MeshTypes::Vec2 closest_UV;
                size_t best_triangle_index = 0;

                for (size_t triangle_index = 0; triangle_index < triangles_.size(); ++triangle_index) {
                    const auto& triangle = triangles_[triangle_index];
                    MeshTypes::Vec3 projected_point;
                    double distance_squared;
                    MeshTypes::Vec2 uv;
                    PointTriangleDistance(point, triangle, projected_point, distance_squared, uv);

                    if (distance_squared < projection_distances[i]) {
                        projection_distances[i] = distance_squared;
                        closest_point = projected_point;
                        closest_UV = uv;
                        best_triangle_index = triangle_index;
                    }
                }
                projection_result.vertices[i] = closest_point;
            }

            std::list<size_t> sorted_vertex_indices;
            {
                std::vector<size_t> temp(vertex_number);
                std::iota(temp.begin(), temp.end(), 0);
                std::sort(temp.begin(), temp.end(), [&projection_distances](size_t a, size_t b) {
                    return projection_distances[a] > projection_distances[b];
                    });
                sorted_vertex_indices.assign(temp.begin(), temp.end());
            }

            std::unordered_map<size_t, Patch> current_patches;
            std::vector<std::unordered_set<size_t>> vertex_to_patch_IDs(vertex_number);
            for (size_t face_index = 0; face_index < polygon_mesh_.faces.size(); ++face_index) {
                current_patches[face_index].boundary = polygon_mesh_.faces[face_index];
                current_patches[face_index].original_faces.insert(face_index);
                for (const auto& vertexIndex : polygon_mesh_.faces[face_index]) {
                    vertex_to_patch_IDs[vertexIndex].insert(face_index);
                }
            }

            size_t next_patch_ID = polygon_mesh_.faces.size();
            bool merge_happened;
            std::unordered_map<std::pair<size_t, size_t>, size_t, HashEdge> edge_counts;
            std::unordered_set<size_t> unique_vertices;
            std::unordered_map<size_t, std::vector<size_t>> adjacency;
            std::vector<std::pair<size_t, size_t>> boundary_edges;
            std::unordered_set<size_t> affected_vertices;
            do {
                merge_happened = false;
                for (auto it = sorted_vertex_indices.begin(); it != sorted_vertex_indices.end();) {
                    const auto vertex_index = *it;
                    if (!is_vertex_active_aux[vertex_index]) {
                        it = sorted_vertex_indices.erase(it);
                        continue;
                    }
                    const auto candidate_patch_IDs = vertex_to_patch_IDs[vertex_index];
                    edge_counts.clear();
                    unique_vertices.clear();
                    adjacency.clear();
                    boundary_edges.clear();
                    for (const auto& patchID : candidate_patch_IDs) {
                        const auto& boundary = current_patches[patchID].boundary;
                        const auto n = boundary.size();
                        for (size_t i = 0; i < n; ++i) {
                            auto u = boundary[i];
                            unique_vertices.insert(u);
                            auto v = boundary[(i + 1) % n];
                            if (u > v) {
                                std::swap(u, v);
                            }
                            ++edge_counts[{u, v}];
                        }
                    }
                    for (const auto& kv : edge_counts) {
                        if (kv.second == 1) {
                            boundary_edges.emplace_back(kv.first);
                        }
                    }
                    if (unique_vertices.size() - edge_counts.size() + candidate_patch_IDs.size() != 1) {
                        ++it;
                        continue;
                    }
                    for (const auto& edge : boundary_edges) {
                        adjacency[edge.first].push_back(edge.second);
                        adjacency[edge.second].push_back(edge.first);
                    }
                    Patch new_patch;
                    size_t start_vertex = adjacency.begin()->first;
                    size_t current_vertex = start_vertex;
                    size_t previous_vertex = adjacency[current_vertex][0];
                    do {
                        new_patch.boundary.push_back(current_vertex);
                        const auto& neighbors = adjacency[current_vertex];
                        const auto next_vertex = (neighbors[0] == previous_vertex ? neighbors[1] : neighbors[0]);
                        previous_vertex = current_vertex;
                        current_vertex = next_vertex;
                    } while (current_vertex != start_vertex);
                    for (const auto& pid : candidate_patch_IDs) {
                        new_patch.original_faces.insert(current_patches[pid].original_faces.begin(),
                            current_patches[pid].original_faces.end());
                    }
                    current_patches[next_patch_ID] = std::move(new_patch);
                    affected_vertices.clear();
                    for (const auto& pid : candidate_patch_IDs) {
                        for (const auto& v : current_patches[pid].boundary) {
                            affected_vertices.insert(v);
                        }
                        current_patches.erase(pid);
                    }
                    for (const auto& v : affected_vertices) {
                        auto& patches = vertex_to_patch_IDs[v];
                        for (const auto& pid : candidate_patch_IDs) {
                            patches.erase(pid);
                        }
                        patches.insert(next_patch_ID);
                        if (patches.size() == 1) {
                            is_vertex_active_aux[v] = false;
                        }
                    }
                    merge_happened = true;
                    ++next_patch_ID;
                    it = sorted_vertex_indices.erase(it);
                }
            } while (merge_happened);

            // Step 3: 将最终的patch（每个patch的origFaces非空）写入文件a.txt
            std::ofstream fout("a.txt");
            for (const auto& [patchID, patch] : current_patches) {
                for (auto faceID : patch.original_faces)
                    fout << faceID << ",";
                fout << "\n\n";
            }
            fout.close();
            
            return projection_result;
        }

    private:
        static constexpr double INF = std::numeric_limits<double>::max();

        const MeshTypes::Mesh3D& triangle_mesh_;
        const MeshTypes::Mesh3D& polygon_mesh_;

        const std::vector<std::pair<size_t, size_t>>& sharp_edges_;
        const std::vector<size_t>& sharp_vertices_;

        std::vector<MeshTypes::Triangle> triangles_;

        std::vector<bool> is_vertex_active_;

        std::vector<std::vector<size_t>> polygon_sharp_edge_vertices_;

        std::unordered_map<size_t, size_t> polygon_sharp_vertex_indices;

        struct Patch {
            MeshTypes::VecNI boundary;
            std::unordered_set<size_t> original_faces;
        };

        struct HashEdge {
            std::size_t operator()(const std::pair<size_t, size_t>& p) const {
                return HASH_EDGE(p.first, p.second);
            }
        };

        void PointTriangleDistance(const MeshTypes::Vec3& p, const MeshTypes::Triangle& tri,
            MeshTypes::Vec3& q, double& dist_sq, MeshTypes::Vec2& uv) {
            const MeshTypes::Vec3 ap = p - tri.a;
            const double d1 = tri.ab.dot(ap);
            const double d2 = tri.ac.dot(ap);
            MeshTypes::Vec3 pq;
            if (d1 <= 0 && d2 <= 0) {
                q = tri.a;
                uv = { 0, 0 };
                pq = p - q;
                dist_sq = pq.dot(pq);
                return;
            }

            const MeshTypes::Vec3 bp = p - (tri.a + tri.ab);
            const double d3 = tri.ab.dot(bp);
            const double d4 = tri.ac.dot(bp);
            if (d3 >= 0 && d4 <= d3) {
                q = tri.a + tri.ab;
                uv = { 1, 0 };
                pq = p - q;
                dist_sq = pq.dot(pq);
                return;
            }

            const MeshTypes::Vec3 cp = p - (tri.a + tri.ac);
            const double d5 = tri.ab.dot(cp);
            const double d6 = tri.ac.dot(cp);
            if (d6 >= 0 && d5 <= d6) {
                q = tri.a + tri.ac;
                uv = { 0, 1 };
                pq = p - q;
                dist_sq = pq.dot(pq);
                return;
            }

            const double vc = d1 * d4 - d3 * d2;
            if (vc <= 0 && d1 >= 0 && d3 <= 0) {
                const double v = d1 / (d1 - d3);
                q = tri.a + tri.ab * v;
                uv = { v, 0 };
                pq = p - q;
                dist_sq = pq.dot(pq);
                return;
            }

            const double vb = d5 * d2 - d1 * d6;
            if (vb <= 0 && d2 >= 0 && d6 <= 0) {
                const double w = d2 / (d2 - d6);
                q = tri.a + tri.ac * w;
                uv = { 0, w };
                pq = p - q;
                dist_sq = pq.dot(pq);
                return;
            }

            const double va = d3 * d6 - d5 * d4;
            if (va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0) {
                const const double v = (d4 - d3) / ((d4 - d3) + (d5 - d6));
                q = tri.a + tri.ab + tri.bc * v;
                uv = { 1 - v, v };
                pq = p - q;
                dist_sq = pq.dot(pq);
                return;
            }

            const double denom = 1.0 / (va + vb + vc);
            const double v = vb * denom;
            const double w = vc * denom;
            q = tri.a + tri.ab * v + tri.ac * w;
            uv = { v, w };
            pq = p - q;
            dist_sq = pq.dot(pq);
        }
    };
}

#endif