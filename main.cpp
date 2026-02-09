#include "Mesh_Input_Output.h"
#include "Mesh_Optimization.h"
#include "Mesh_Shape_Diameter_Function.h"

void main() {
	MeshTypes::Mesh3D triangle_mesh = MeshInput::ReadTriangleMesh("ansys_bracket.obj");
	MeshTypes::SharpFeature sharp_feature = MeshInput::ReadSharpFeature("i08o_m8.txt");
	MeshTypes::Mesh3D hexahedra_mesh = MeshInput::ReadHexahedronMesh("ansys_bracket-hexMesh-interiorHexMesh.vtk");
	MeshOptimization::MeshOptimizer mesh_optimizer(triangle_mesh, sharp_feature, hexahedra_mesh);
	const double learning_rate = 4e-4;
	const double lambda_projection = 1e2;
	const double max_distance_threshold = 1e-10;
	MeshTypes::Mesh3D optimized_hexahedra_mesh = mesh_optimizer.Optimize(
		learning_rate,
		lambda_projection,
		max_distance_threshold);
	MeshOutput::WriteHexahedraMesh("insideHexOut.vtk", optimized_hexahedra_mesh);
}

//#include "Mesh_Embedding.h"
//#include "Mesh_Projection.h"
//#include "Mesh_Shortest_Path.h"
//#include <fstream>
//#include <sstream>
//#include <iostream>
//#include <iomanip>
//#include <chrono>
//
//int main() {
//    MeshTypes::Mesh3D triangle_mesh, polygon_mesh;
//
//    MeshTypes::Vec3 vec3;
//
//    MeshTypes::VecNI vec3i;
//    vec3i.resize(3);
//
//    // Open the polygon mesh file
//    std::ifstream inFile("cad.obj");
//
//    // A vertex counter
//    size_t vertexCounter = 0;
//
//    // Cache of a line
//    std::string lineCache;
//    std::istringstream strCache;
//
//    // Cache of a char
//    char type;
//
//    // Cache of vertex indices in a facet
//    size_t vertIdxCache0, vertIdxCache1, vertIdxCache2;
//
//    // Read vertices
//    while (std::getline(inFile, lineCache)) {
//        strCache.clear();
//        strCache.str(lineCache);
//        if (!(strCache >> type >> vec3.x >> vec3.y >> vec3.z) || type != 'v') {
//            break;
//        }
//        triangle_mesh.vertices.push_back(vec3);
//        ++vertexCounter;
//    }
//    strCache.clear();
//    strCache.str(lineCache);
//    if ((strCache >> type >> vec3i[0] >> vec3i[1] >>
//        vec3i[2])) {
//        --vec3i[0];
//        --vec3i[1];
//        --vec3i[2];
//        triangle_mesh.faces.push_back(vec3i);
//    }
//
//    while (std::getline(inFile, lineCache)) {
//        strCache.clear();
//        strCache.str(lineCache);
//        if (!(strCache >> type >> vec3i[0] >> vec3i[1] >>
//            vec3i[2])) {
//            continue;
//        }
//        --vec3i[0];
//        --vec3i[1];
//        --vec3i[2];
//        triangle_mesh.faces.push_back(vec3i);
//    }
//
//    // Close the input file
//    inFile.close(); 
//
//    inFile.open("block_surface.obj");
//
//    // A vertex counter
//    vertexCounter = 0;
//
//    // Read vertices
//    while (std::getline(inFile, lineCache)) {
//        strCache.clear();
//        strCache.str(lineCache);
//        if (!(strCache >> type >> vec3.x >> vec3.y >> vec3.z) || type != 'v') {
//            break;
//        }
//        polygon_mesh.vertices.push_back(vec3);
//        ++vertexCounter;
//    }
//    strCache.clear();
//    strCache.str(lineCache);
//    size_t tmp;
//    if ((strCache >> type >> vec3i[0] >> vec3i[1] >>
//        vec3i[2])) {
//        --vec3i[0];
//        --vec3i[1];
//        --vec3i[2];
//        while (strCache >> tmp) {
//            vec3i.push_back(tmp - 1);
//        }
//        polygon_mesh.faces.push_back(vec3i);
//    }
//
//    while (std::getline(inFile, lineCache)) {
//        vec3i.resize(3);
//        strCache.clear();
//        strCache.str(lineCache);
//        if (!(strCache >> type >> vec3i[0] >> vec3i[1] >>
//            vec3i[2])) {
//            continue;
//        }
//        --vec3i[0];
//        --vec3i[1];
//        --vec3i[2];
//        while (strCache >> tmp) {
//            vec3i.push_back(tmp - 1);
//        }
//        polygon_mesh.faces.push_back(vec3i);
//    }
//
//    // Close the input file
//    inFile.close();
//    MeshTypes::SharpFeature ftr;
//    ftr.edges = { {97,1488},{1488,466},{466,175},{175,97},{494,450},{450,226},{226,576},{576,494},{97,494},{1488,450},{466,226},{175,576} };
//    ftr.vertices = { 407 };
//
//    auto start = std::chrono::high_resolution_clock::now();
//    MeshProjection::MeshProjector proj(triangle_mesh, polygon_mesh, ftr);
//    auto result = proj.ComputeProjection();
//
//
//    auto end = std::chrono::high_resolution_clock::now();
//    std::chrono::duration<double, std::milli> duration = end - start;
//    std::cout << "Elapsed time: " << duration.count() << " ms" << std::endl;
//    std::ofstream outFile("block_projected.obj");
//    outFile << std::setprecision(std::numeric_limits<double>::max_digits10);
//    for (const auto& vertex : result.vertices) {
//        outFile << "v " << vertex.x << " " << vertex.y << " " << vertex.z << std::endl;
//    }
//
//    for (const auto& triangle : polygon_mesh.faces) {
//        outFile << "f " << triangle[0] + 1 << " " << triangle[1] + 1 << " " << triangle[2] + 1 << " "<<triangle[3]+1 << std::endl;
//    }
//    outFile.close();
//
//	//MeshTypes::Mesh3D mesh;
//
// //   MeshTypes::Vec3 vec3;
//
// //   MeshTypes::VecNI vec3i;
// //   vec3i.resize(3);
//
// //   // Open the polygon mesh file
// //   std::ifstream inFile("lionT.obj");
//
// //   // A vertex counter
// //   size_t vertexCounter = 0;
//
// //   // Cache of a line
// //   std::string lineCache;
// //   std::istringstream strCache;
//
// //   // Cache of a char
// //   char type;
//
// //   // Cache of vertex indices in a facet
// //   size_t vertIdxCache0, vertIdxCache1, vertIdxCache2;
//
// //   // Read vertices
// //   while (std::getline(inFile, lineCache)) {
// //       strCache.clear();
// //       strCache.str(lineCache);
// //       if (!(strCache >> type >> vec3.x >> vec3.y >> vec3.z) || type != 'v') {
// //           break;
// //       }
// //       mesh.vertices.push_back(vec3);
// //       ++vertexCounter;
// //   }
// //   strCache.clear();
// //   strCache.str(lineCache);
// //   if ((strCache >> type >> vec3i[0] >> vec3i[1] >>
// //       vec3i[2])) {
// //       --vec3i[0];
// //       --vec3i[1];
// //       --vec3i[2];
// //       mesh.faces.push_back(vec3i);
// //   }
//
// //   while (std::getline(inFile, lineCache)) {
// //       strCache.clear();
// //       strCache.str(lineCache);
// //       if (!(strCache >> type >> vec3i[0] >> vec3i[1] >>
// //           vec3i[2])) {
// //           continue;
// //       }
// //       --vec3i[0];
// //       --vec3i[1];
// //       --vec3i[2];
// //       mesh.faces.push_back(vec3i);
// //   }
//
// //   // Close the input file
// //   inFile.close();
//    //////TEST SHORTEST PATH
//    /*std::vector<bool> valid_flags(mesh.vertices.size(), 1);
//    valid_flags[7293] = 0;
//    valid_flags[7345] = 0;
//    valid_flags[7296] = 0;
//    MeshShortestPath::DijkstraSolver solver(mesh, valid_flags);
//    auto path = solver.ComputePath(10504, 7705);
//    for (const auto& i : path) {
//        std::cout << i << ",";
//    }*/
//    //////TEST PARAMETERIZATION
//    /*MeshEmbedding::MeshProcessor::UnifyTriangleOrientations(mesh);
//
//    std::vector<std::vector<size_t>> boundaryIndices(13);
//    boundaryIndices[0] = { 11640,11585 };
//    boundaryIndices[1] = { 11585,11624,11623 };
//    boundaryIndices[2] = { 11623,11615,11616,11625 };
//    boundaryIndices[3] = { 11625,11626,11628,11637,11636 };
//    boundaryIndices[4] = { 11636,11635,11631,11632,11649,8219,8220 };
//    boundaryIndices[5] = { 8220,8221,8222,8224,8218,8122,8120,8209 };
//    boundaryIndices[6] = { 8209,8208,8203,8202,8195,8194,8196,8197,8982 };
//    boundaryIndices[7] = { 8982,8993,8994,8995,8923,8922,8921,8929,8872,8869 };
//    boundaryIndices[8] = { 8869,8868,8864,8865,12590,12588,12589,12664,12663,12662,12752 };
//    boundaryIndices[9] = { 12752,12750,12749,12711,12709,12708,12707,12706,12676,12674,12671,12670 };
//    boundaryIndices[10] = { 12670,12673,12763,5741,5739,5738,5737,5742,5743,5725,5724,5723,5777 };
//    boundaryIndices[11] = { 5777,5776,5775,5774,4784,4769,4771,4770,4777,4774,4775,4749,4748 };
//    boundaryIndices[12] = { 4748,4736,4735,4734,4726,4728,4727,4731,4732,4768,11542,11642,11641,11640 };
//    auto start = std::chrono::high_resolution_clock::now();
//
//    MeshEmbedding::TutteEmbedder embedder(mesh);
//    MeshTypes::Mesh2D parameterized_mesh = embedder.ComputeEmbedding(boundaryIndices);
//
//    auto end = std::chrono::high_resolution_clock::now();
//
//    std::chrono::duration<double, std::milli> duration = end - start;
//    std::cout << "Elapsed time: " << duration.count() << " ms" << std::endl;
//    std::ofstream outFile("lion2.obj");
//    outFile << std::setprecision(std::numeric_limits<double>::max_digits10);
//    for (const auto& vertex : parameterized_mesh.vertices) {
//        outFile << "v " << vertex.x << " " << vertex.y << " 0" << std::endl;
//    }
//    
//    for (const auto& triangle : mesh.faces) {
//        outFile << "f " << triangle[0] + 1 << " " << triangle[1] + 1 << " " << triangle[2] + 1 << std::endl;
//    }
//    outFile.close();*/
//    return 0;
//}