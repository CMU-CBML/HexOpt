#ifndef MESH_INPUT_OUTPUT_H
#define MESH_INPUT_OUTPUT_H

#include "Mesh_Types.h"
#include <fstream>
#include <iomanip>
#include <sstream>
#include <unordered_map>

namespace MeshInput {
    MeshTypes::Mesh3D ReadTriangleMesh(const std::string& triangle_file_name) {
        MeshTypes::Mesh3D mesh;
        std::ifstream input_file(triangle_file_name);
        if (!input_file.is_open()) {
            throw std::runtime_error("Cannot read input triangle mesh: Failed to open input file.");
        }

        std::unordered_map<size_t, size_t> chk_dup_pt, map_dup_pt;
        size_t vertex_counter = 1;
        std::string line_cache;
        char type;

        while (std::getline(input_file, line_cache)) {
            if (line_cache.empty() || (line_cache[0] != 'v' && line_cache[0] != 'V')) {
                if (line_cache[0] == 'f' || line_cache[0] == 'F') {
                    std::istringstream str_cache(line_cache);
                    size_t vert_idx_cache0, vert_idx_cache1, vert_idx_cache2;
                    if (str_cache >> type >> vert_idx_cache0 >> vert_idx_cache1 >> vert_idx_cache2) {
                        MeshTypes::VecNI face;
                        face.push_back(map_dup_pt[vert_idx_cache0]);
                        face.push_back(map_dup_pt[vert_idx_cache1]);
                        face.push_back(map_dup_pt[vert_idx_cache2]);
                        mesh.faces.push_back(face);
                    }
                    break;
                }
                continue;
            }

            std::istringstream str_cache(line_cache);
            double x, y, z;
            if (str_cache >> type >> x >> y >> z) {
                size_t hash = HASH_VERTEX(x, y, z);
                if (!chk_dup_pt.count(hash)) {
                    mesh.vertices.emplace_back(x, y, z);
                    chk_dup_pt[hash] = vertex_counter - 1;
                    map_dup_pt[vertex_counter] = vertex_counter - 1;
                }
                else {
                    map_dup_pt[vertex_counter] = chk_dup_pt[hash];
                }
                ++vertex_counter;
            }
        }

        while (std::getline(input_file, line_cache)) {
            if (line_cache.empty() || (line_cache[0] != 'f' && line_cache[0] != 'F')) {
                continue;
            }

            std::istringstream str_cache(line_cache);
            size_t vert_idx_cache0, vert_idx_cache1, vert_idx_cache2;
            if (str_cache >> type >> vert_idx_cache0 >> vert_idx_cache1 >> vert_idx_cache2) {
                MeshTypes::VecNI face;
                face.push_back(map_dup_pt[vert_idx_cache0]);
                face.push_back(map_dup_pt[vert_idx_cache1]);
                face.push_back(map_dup_pt[vert_idx_cache2]);
                mesh.faces.push_back(face);
            }
        }

        input_file.close();

        return mesh;
    }

    MeshTypes::Mesh3D ReadHexahedronMesh(const std::string& hex_file_name) {
        MeshTypes::Mesh3D mesh;
        std::ifstream input_file(hex_file_name);
        if (!input_file.is_open()) {
            throw std::runtime_error("Cannot read input hex mesh: Failed to open input file.");
        }

        std::string line;
        for (size_t i = 0; i < 4 && std::getline(input_file, line); ++i) {}

        size_t vertex_num;
        if (std::getline(input_file, line)) {
            std::istringstream iss(line);
            std::string keyword;
            if (iss >> keyword >> vertex_num) {
                for (size_t i = 0; i < vertex_num; ++i) {
                    if (std::getline(input_file, line)) {
                        std::istringstream pss(line);
                        double x, y, z;
                        if (pss >> x >> y >> z) {
                            mesh.vertices.emplace_back(x, y, z);
                        }
                    }
                }
            }
        }

        size_t element_num, dummy;
        if (std::getline(input_file, line)) {
            std::istringstream iss(line);
            std::string keyword;
            if (iss >> keyword >> element_num >> dummy) {
                for (size_t i = 0; i < element_num; ++i) {
                    if (std::getline(input_file, line)) {
                        std::istringstream cell_stream(line);
                        size_t count;
                        cell_stream >> count;
                        if (count != 8) {
                            throw std::runtime_error("Invalid cell data: Expected 8 vertices for a hexahedron.");
                        }
                        MeshTypes::VecNI cell;
                        size_t idx;
                        for (size_t j = 0; j < 8; ++j) {
                            if (cell_stream >> idx) {
                                cell.push_back(idx);
                            }
                        }
                        mesh.faces.push_back(cell);
                    }
                }
            }
        }

        input_file.close();
        return mesh;
    }

    MeshTypes::SharpFeature ReadSharpFeature(const std::string& feature_file_name) {
        MeshTypes::SharpFeature feature;
        std::ifstream input_file(feature_file_name);
        if (!input_file.is_open()) {
            throw std::runtime_error("Cannot read sharp feature file: Failed to open input file.");
        }

        std::string line;
        while (std::getline(input_file, line)) {
            if (line.empty()) {
                continue;
            }

            std::istringstream iss(line);
            std::vector<size_t> tokens;
            size_t value;
            while (iss >> value) {
                tokens.push_back(value);
            }

            if (tokens.size() == 1) {
                feature.vertices.push_back(tokens[0]);
            }
            else if (tokens.size() == 2) {
                feature.edges.emplace_back(tokens[0], tokens[1]);
            }
        }

        input_file.close();
        return feature;
    }
}

namespace MeshOutput {
    void WriteHexahedraMesh(const std::string& hexahedra_file_name,
        const MeshTypes::Mesh3D& hexahedra_mesh_) {
        std::ofstream output_file(hexahedra_file_name);
        if (!output_file.is_open()) {
            throw std::runtime_error("Cannot write output hexahedra mesh: Failed to open output file.");
        }

        output_file << std::setprecision(std::numeric_limits<double>::max_digits10);

        output_file << "# vtk DataFile Version 3.0\n";
        output_file << "Hexahedral Mesh\n";
        output_file << "ASCII\n";
        output_file << "DATASET UNSTRUCTURED_GRID\n";

        output_file << "POINTS " << hexahedra_mesh_.vertices.size() << " double\n";
        for (const auto& vertex : hexahedra_mesh_.vertices) {
            output_file << vertex.x << " " << vertex.y << " " << vertex.z << "\n";
        }

        const size_t num_hex = hexahedra_mesh_.faces.size();
        const size_t total_cell_data = num_hex * (1 + 8);

        output_file << "CELLS " << num_hex << " " << total_cell_data << "\n";
        for (const auto& hex : hexahedra_mesh_.faces) {
            output_file << "8";
            for (const auto& idx : hex) {
                output_file << " " << idx;
            }
            output_file << "\n";
        }

        output_file << "CELL_TYPES " << num_hex << "\n";
        for (size_t i = 0; i < num_hex; ++i) {
            output_file << "12\n";
        }

        output_file.close();
    }
}

#endif