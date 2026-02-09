#ifndef MESH_SHORTEST_PATH_H
#define MESH_SHORTEST_PATH_H

#include "Mesh_Types.h"
#include <queue>
#include <unordered_map>

namespace MeshShortestPath {
    class DijkstraSolver {
    public:
        DijkstraSolver(const MeshTypes::Mesh3D& mesh,
            std::vector<bool>& valid_vertices)
            : mesh_(mesh),
            valid_vertices_(valid_vertices),
            adjacency_list_(mesh.vertices.size()) {
            buildAdjacencyList();
        }

        std::vector<size_t> ComputePath(const size_t start, const size_t end, const bool track_change = false) {
            if (start == end) {
                return { start, end };
            }
            const size_t num_vertices = mesh_.vertices.size();
            std::vector<double> dist(num_vertices, INF);
            std::vector<size_t> predecessors(num_vertices);
            std::priority_queue<Node, std::vector<Node>, std::greater<Node>> priority_queue;

            dist[start] = 0.0;
            priority_queue.emplace(start, 0.0);

            while (!priority_queue.empty()) {
                const Node current = priority_queue.top();
                priority_queue.pop();

                if (current.index == end) {
                    break;
                }
                if (current.distance > dist[current.index]) {
                    continue;
                }

                for (const auto& neighbor : adjacency_list_[current.index]) {
                    if (!valid_vertices_[neighbor.first] && neighbor.first != end) {
                        continue;
                    }

                    const double new_dist = dist[current.index] + neighbor.second;
                    if (new_dist < dist[neighbor.first]) {
                        dist[neighbor.first] = new_dist;
                        predecessors[neighbor.first] = current.index;
                        priority_queue.emplace(neighbor.first, new_dist);
                    }
                }
            }

            if (dist[end] == INF) {
                return {};
            }

            return reconstructPath(predecessors, start, end, track_change);
        }

    private:
        struct Node {
            size_t index;
            double distance;

            Node(size_t idx, double dist) : index(idx), distance(dist) {
            }
            bool operator>(const Node& other) const {
                return distance > other.distance;
            }
        };

        static constexpr double INF = std::numeric_limits<double>::max();

        const MeshTypes::Mesh3D& mesh_;
        std::vector<bool>& valid_vertices_;
        std::vector<std::vector<std::pair<size_t, double>>> adjacency_list_;

        void buildAdjacencyList() {
            std::unordered_map<size_t, std::unordered_map<size_t, double>> adjacency_weights;

            for (const auto& face : mesh_.faces) {
                const size_t face_size = face.size();
                for (size_t i = 0; i < face_size; ++i) {
                    const size_t src = face[i];
                    const size_t dst = face[(i + 1) % face_size];

                    const double weight = (mesh_.vertices[src] - mesh_.vertices[dst]).norm();
                    adjacency_weights[src][dst] = weight;
                    adjacency_weights[dst][src] = weight;
                }
            }

            adjacency_list_.resize(mesh_.vertices.size());
            for (const auto& [src, neighbors] : adjacency_weights) {
                adjacency_list_[src].reserve(neighbors.size());
                for (const auto& [dst, weight] : neighbors) {
                    adjacency_list_[src].emplace_back(dst, weight);
                }
            }
        }

        std::vector<size_t> reconstructPath(const std::vector<size_t>& predecessors,
            const size_t start, const size_t end, const bool track_change) const {
            std::vector<size_t> path;
            size_t current = predecessors[end];
            while (current != start) {
                path.push_back(current);
                current = predecessors[current];
            }
            path.push_back(start);
            std::reverse(path.begin(), path.end());
            path.push_back(end);
            if (track_change) {
                for (const auto& i : path) {
                    valid_vertices_[i] = false;
                }
            }

            return path;
        }
    };
}

#endif