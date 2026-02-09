#ifndef MESH_EMBEDDING_H
#define MESH_EMBEDDING_H

#include "Mesh_Types.h"
#include <algorithm>
#include <omp.h>
#include <queue>
#include <unordered_map>
#include <unordered_set>

namespace MeshEmbedding {
    class MeshProcessor {
    public:
        static void UnifyTriangleOrientations(MeshTypes::Mesh3D& mesh) {
            const size_t num_triangles = mesh.faces.size();
            std::unordered_set<size_t> remaining_triangles;
            for (size_t i = 0; i < num_triangles; ++i) {
                remaining_triangles.insert(i);
            }

            std::unordered_map<size_t, std::vector<size_t>> edge_to_triangles;
            edge_to_triangles.reserve(num_triangles * 3);
            for (size_t tri_idx = 0; tri_idx < num_triangles; ++tri_idx) {
                const MeshTypes::VecNI& tri = mesh.faces[tri_idx];
                edge_to_triangles[HASH_SORTED_EDGE(tri[0], tri[1])].push_back(tri_idx);
                edge_to_triangles[HASH_SORTED_EDGE(tri[1], tri[2])].push_back(tri_idx);
                edge_to_triangles[HASH_SORTED_EDGE(tri[2], tri[0])].push_back(tri_idx);
            }

            std::queue<size_t> processing_queue;
            while (!remaining_triangles.empty()) {
                processing_queue.push(*remaining_triangles.begin());
                remaining_triangles.erase(remaining_triangles.begin());

                while (!processing_queue.empty()) {
                    size_t current_idx = processing_queue.front();
                    processing_queue.pop();
                    MeshTypes::VecNI& current_tri = mesh.faces[current_idx];

                    for (size_t edge_num = 0; edge_num < 3; ++edge_num) {
                        const size_t v_start = current_tri[edge_num];
                        const size_t v_end = current_tri[(edge_num + 1) % 3];
                        const size_t edge_hash = HASH_SORTED_EDGE(v_start, v_end);

                        for (const auto& neighbor_idx : edge_to_triangles[edge_hash]) {
                            if (!remaining_triangles.count(neighbor_idx)) {
                                continue;
                            }

                            MeshTypes::VecNI& neighbor_tri = mesh.faces[neighbor_idx];
                            for (size_t i = 0; i < 3; ++i) {
                                size_t& a = neighbor_tri[i];
                                size_t& b = neighbor_tri[(i + 1) % 3];
                                if (a == v_start && b == v_end) {
                                    std::swap(a, b);
                                    break;
                                }
                                if (a == v_end && b == v_start) {
                                    break;
                                }
                            }
                            processing_queue.push(neighbor_idx);
                            remaining_triangles.erase(neighbor_idx);
                        }
                    }
                }
            }
        }
    };

    class MeshValidator {
    public:
        static bool HasSelfIntersections(const MeshTypes::Mesh2D& mesh) {
            bool positive = false;
            bool negative = false;

#pragma omp parallel for
            for (long long i = 0; i < mesh.faces.size(); ++i) {
                const MeshTypes::Vec2& v0 = mesh.vertices[mesh.faces[i][0]];
                const MeshTypes::Vec2& v1 = mesh.vertices[mesh.faces[i][1]];
                const MeshTypes::Vec2& v2 = mesh.vertices[mesh.faces[i][2]];
                const double area = (v1 - v0).cross(v2 - v0);

                if (area > 0) {
                    positive = true;
                    if (negative) {
                        break;
                    }
                }
                else if (area < 0) {
                    negative = true;
                    if (positive) {
                        break;
                    }
                }
                else {
                    positive = true;
                    negative = true;
                    break;
                }
            }
            return (positive && negative);
        }
    };

    class LinearSolver {
    public:
        struct SolverConfig {
            size_t check_interval = 100;
            size_t max_iterations = 100000000;
        };

        static void SolveDualSystems(
            const std::vector<double>& diagonal,
            const std::vector<std::vector<std::pair<size_t, double>>>& off_diagonal,
            const std::vector<double>& rhs_x,
            const std::vector<double>& rhs_y,
            const std::vector<size_t>& internal_vertices,
            MeshTypes::Mesh2D& parameterized_mesh,
            const SolverConfig& config = {}) {
            const size_t system_size = rhs_x.size();
            std::vector<double> x(system_size, 0.0);
            std::vector<double> y(system_size, 0.0);
            std::vector<double> x_new(system_size);
            std::vector<double> y_new(system_size);

            for (size_t iteration = 0; iteration < config.max_iterations; ++iteration) {
#pragma omp parallel for
                for (long long i = 0; i < system_size; ++i) {
                    double sum_x = 0.0, sum_y = 0.0;
                    for (const auto& entry : off_diagonal[i]) {
                        sum_x += entry.second * x[entry.first];
                        sum_y += entry.second * y[entry.first];
                    }
                    x_new[i] = (rhs_x[i] - sum_x) / diagonal[i];
                    y_new[i] = (rhs_y[i] - sum_y) / diagonal[i];
                }

                x.swap(x_new);
                y.swap(y_new);

                if (iteration % config.check_interval == 0) {
#pragma omp parallel for
                    for (long long i = 0; i < system_size; ++i) {
                        parameterized_mesh.vertices[internal_vertices[i]] = MeshTypes::Vec2(x[i], y[i]);
                    }
                    if (!MeshValidator::HasSelfIntersections(parameterized_mesh)) {
                        break;
                    }
                }
            }
        }
    };

    class TutteEmbedder {
    public:
        struct EmbeddingConfig {
            size_t weight_type = 1;
            LinearSolver::SolverConfig solver_config;
        };

        explicit TutteEmbedder(const MeshTypes::Mesh3D& input_mesh)
            : source_mesh_(input_mesh) {
        }

        MeshTypes::Mesh2D ComputeEmbedding(
            const std::vector<std::vector<size_t>>& boundary_chains,
            const EmbeddingConfig& config = {}) const {
            const size_t num_vertices = source_mesh_.vertices.size();
            MeshTypes::Mesh2D parameterized_mesh;
            parameterized_mesh.vertices.resize(num_vertices);
            parameterized_mesh.faces = source_mesh_.faces;

            const size_t n = boundary_chains.size();
            const double two_pi = 6.28318530717958647692;

            for (size_t i = 0; i < n; ++i) {
                const auto& chain = boundary_chains[i];
                const size_t corner = chain[0];
                const double angle = two_pi * i / n;
                const MeshTypes::Vec2 start(std::cos(angle), std::sin(angle));
                parameterized_mesh.vertices[corner] = start;

                const size_t next_i = (i + 1) % n;
                const double next_angle = two_pi * next_i / n;
                const MeshTypes::Vec2 end(std::cos(next_angle), std::sin(next_angle));

                auto InterpolateBoundary = [&](const std::vector<size_t>& chain) {
                    std::vector<double> lengths(chain.size(), 0.0);
                    for (size_t j = 1; j < chain.size(); ++j) {
                        const MeshTypes::Vec3& prev = source_mesh_.vertices[chain[j - 1]];
                        const MeshTypes::Vec3& curr = source_mesh_.vertices[chain[j]];
                        lengths[j] = lengths[j - 1] + (curr - prev).norm();
                    }
                    const double total_length = lengths.back();
                    for (size_t j = 1; j < chain.size() - 1; ++j) {
                        const double t = lengths[j] / total_length;
                        parameterized_mesh.vertices[chain[j]] = MeshTypes::Vec2(
                            start.x + t * (end.x - start.x),
                            start.y + t * (end.y - start.y)
                        );
                    }
                    };

                InterpolateBoundary(boundary_chains[i]);
            }

            std::unordered_set<size_t> boundary_set;
            for (const auto& chain : boundary_chains) {
                for (auto idx : chain) {
                    boundary_set.insert(idx);
                }
            }

            std::vector<std::unordered_set<size_t>> adjacency(num_vertices);
            for (const auto& tri : source_mesh_.faces) {
                for (size_t e = 0; e < 3; ++e) {
                    size_t a = tri[e];
                    size_t b = tri[(e + 1) % 3];
                    adjacency[a].insert(b);
                    adjacency[b].insert(a);
                }
            }

            std::vector<size_t> internal_vertices;
            internal_vertices.reserve(num_vertices - boundary_set.size());
            std::unordered_map<size_t, size_t> vertex_to_index;
            vertex_to_index.reserve(num_vertices - boundary_set.size());
            for (size_t i = 0; i < num_vertices; ++i) {
                if (!boundary_set.count(i) && !adjacency[i].empty()) {
                    vertex_to_index[i] = internal_vertices.size();
                    internal_vertices.push_back(i);
                }
            }

            std::unordered_map<size_t, std::vector<size_t>> edge_opposite_vertices;
            for (const auto& tri : source_mesh_.faces) {
                edge_opposite_vertices[HASH_SORTED_EDGE(tri[0], tri[1])].push_back(tri[2]);
                edge_opposite_vertices[HASH_SORTED_EDGE(tri[1], tri[2])].push_back(tri[0]);
                edge_opposite_vertices[HASH_SORTED_EDGE(tri[2], tri[0])].push_back(tri[1]);
            }

            auto ComputeWeight = [&](size_t i, size_t j) {
                if (config.weight_type == 0) {
                    return 1.0;
                }
                else if (config.weight_type == 1) {
                    const MeshTypes::Vec3& vi = source_mesh_.vertices[i];
                    const MeshTypes::Vec3& vj = source_mesh_.vertices[j];
                    const MeshTypes::Vec3 vec_ij = vj - vi;
                    const double length_ij = vec_ij.norm();
                    double total_tangent = 0.0;

                    const auto& opposite_verts = edge_opposite_vertices[HASH_SORTED_EDGE(i, j)];
                    for (const auto& k : opposite_verts) {
                        const MeshTypes::Vec3 vec_ik = source_mesh_.vertices[k] - vi;
                        const double length_ik = vec_ik.norm();
                        const double cos_theta = vec_ij.dot(vec_ik) / (length_ij * length_ik);
                        total_tangent += tan(acos(std::clamp(cos_theta, -1.0, 1.0)) / 2.0);
                    }
                    return total_tangent / length_ij;
                }
                };

            const size_t num_internal = internal_vertices.size();
            std::vector<double> diagonal(num_internal, 0.0);
            std::vector<std::vector<std::pair<size_t, double>>> off_diagonal(num_internal);
            std::vector<double> rhs_x(num_internal, 0.0);
            std::vector<double> rhs_y(num_internal, 0.0);

            for (const auto& vertex : internal_vertices) {
                const size_t row = vertex_to_index[vertex];
                double weight_sum = 0.0;

                for (const auto& neighbor : adjacency[vertex]) {
                    const double weight = ComputeWeight(vertex, neighbor);
                    weight_sum += weight;

                    if (boundary_set.count(neighbor)) {
                        rhs_x[row] -= weight * parameterized_mesh.vertices[neighbor].x;
                        rhs_y[row] -= weight * parameterized_mesh.vertices[neighbor].y;
                    }
                    else {
                        const size_t neighbor_index = vertex_to_index[neighbor];
                        off_diagonal[row].emplace_back(neighbor_index, weight);
                    }
                }
                diagonal[row] = -weight_sum;
            }

            LinearSolver::SolveDualSystems(diagonal, off_diagonal, rhs_x, rhs_y,
                internal_vertices, parameterized_mesh,
                config.solver_config);

            return parameterized_mesh;
        }

    private:
        const MeshTypes::Mesh3D& source_mesh_;
    };
}

#endif