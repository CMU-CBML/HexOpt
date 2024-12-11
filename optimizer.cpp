#include <ctime>
#include <iostream>
#include <omp.h>
#include <random>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include "constants.h"
#include "geometry.h"
#include "meshQuality.h"

int main() {
    // parallel computing
    omp_set_dynamic(1);

    // filename
    const char* triFileName = "isidore_horse.obj";
    const char* hexFileName = "isidore_horse.vtk";
    const char* ftrFileName = "mid2Fem.txt";
    const char* outputHexFileName = "isidore_horse2Opt.vtk";
    // declare verts and tris
    std::vector<std::vector<int>> tri;
    std::vector<std::vector<double>> triX;
    std::vector<std::vector<std::vector<double>>> triEdgeX;

    // read tris
    FILE* dataFile = fopen(triFileName, "r");
    int triPNum = 1, triENum = 0, ptIdx[8];
    double pX, pY, pZ;
    char line[256];
    while (!(fgets(line, sizeof(line), dataFile) && sscanf(line, "v %lf %lf %lf", &pX, &pY, &pZ) == 3)) {}
    triX.push_back({ pX, pY, pZ });
    while (fgets(line, sizeof(line), dataFile) && sscanf(line, "v %lf %lf %lf", &pX, &pY, &pZ) == 3) {
        triX.push_back({ pX, pY, pZ });
        ++triPNum;
    }
    while (!(sscanf(line, "f %zu %zu %zu", &ptIdx[0], &ptIdx[1], &ptIdx[2]) == 3))
        fgets(line, sizeof(line), dataFile);
    while (sscanf(line, "f %d %d %d", &ptIdx[0], &ptIdx[1], &ptIdx[2]) == 3) {
        --ptIdx[0];
        --ptIdx[1];
        --ptIdx[2];
        tri.push_back({ ptIdx[0], ptIdx[1], ptIdx[2] });
        triEdgeX.push_back({ { triX[ptIdx[1]][0] - triX[ptIdx[0]][0], triX[ptIdx[1]][1] - triX[ptIdx[0]][1], triX[ptIdx[1]][2] - triX[ptIdx[0]][2] },
            { triX[ptIdx[2]][0] - triX[ptIdx[0]][0], triX[ptIdx[2]][1] - triX[ptIdx[0]][1], triX[ptIdx[2]][2] - triX[ptIdx[0]][2] },
            { triX[ptIdx[2]][0] - triX[ptIdx[1]][0], triX[ptIdx[2]][1] - triX[ptIdx[1]][1], triX[ptIdx[2]][2] - triX[ptIdx[1]][2] } });
        ++triENum;
        if (!fgets(line, sizeof(line), dataFile)) break;
    }
    fclose(dataFile);
    
    // declare vars
    int i, j;

    // declare verts and hexes
    std::vector<double> x;
    std::vector<std::vector<int>> hex, adjPts;

    // read hex
    dataFile = fopen(hexFileName, "r");
    int pNum, eNum;
    for (i = 0; i < 4 && fgets(line, sizeof(line), dataFile); ++i) {}
    if (fgets(line, sizeof(line), dataFile) && sscanf(line, "POINTS %d double", &pNum) == 1) {
        for (i = 0; i < pNum; ++i) {
            fgets(line, sizeof(line), dataFile);
            sscanf(line, "%lf %lf %lf", &pX, &pY, &pZ);
            x.insert(x.end(), { pX, pY, pZ });
        }
        if (fgets(line, sizeof(line), dataFile) && sscanf(line, "CELLS %d %d", &eNum, &i) == 2) {
            adjPts.resize(pNum);
            std::vector<std::unordered_set<size_t>> adjPtsD3HashTable(pNum);
            for (i = 0; i < eNum; ++i) {
                fgets(line, sizeof(line), dataFile);
                sscanf(line, "%d %d %d %d %d %d %d %d %d", &j, &ptIdx[0], &ptIdx[1], &ptIdx[2], &ptIdx[3], &ptIdx[4], &ptIdx[5], &ptIdx[6], &ptIdx[7]);
                hex.push_back({ 3 * ptIdx[0], 3 * ptIdx[1], 3 * ptIdx[2], 3 * ptIdx[3], 3 * ptIdx[4], 3 * ptIdx[5], 3 * ptIdx[6], 3 * ptIdx[7] });
                for (j = 0; j < 8; ++j)
                    for (const auto& k : adjPt[j])
                        if (adjPtsD3HashTable[ptIdx[j]].insert(ptIdx[k]).second)
                            adjPts[ptIdx[j]].push_back(ptIdx[k] * 3);
            }
        }
    }
    fclose(dataFile);
    
    // declare ftrs
    std::vector<std::vector<double>> edgeFtrX;
    std::vector<std::vector<int>> hexPtType(pNum, std::vector<int>(1, 0));
    double dist[3];

    // read ftrs
    dataFile = fopen(ftrFileName, "r");
    fgets(line, sizeof(line), dataFile);
    while (fgets(line, sizeof(line), dataFile) && sscanf(line, "%d %d", &i, &j) == 2) {
        hexPtType[j] = { -2, i };
    }
    while (fgets(line, sizeof(line), dataFile) && sscanf(line, "%d %d", &i, &j) == 2) {
        dist[0] = triX[j][0] - triX[i][0];
        dist[1] = triX[j][1] - triX[i][1];
        dist[2] = triX[j][2] - triX[i][2];
        edgeFtrX.push_back({ triX[i][0], triX[i][1], triX[i][2], triX[j][0], triX[j][1], triX[j][2], dist[0], dist[1], dist[2], DOT(dist, dist)});
    }
    char* ptr;
    while (fgets(line, sizeof(line), dataFile)) {
        ptr = line;
        sscanf(ptr, "%d", &i);
        hexPtType[i] = { -1 };
        while (*ptr != '\0' && *ptr != ' ') ++ptr;
        while (*ptr == ' ') ++ptr;
        while (sscanf(ptr, "%d", &j) == 1) {
            hexPtType[i].push_back(j);
            while (*ptr != '\0' && *ptr != ' ') ++ptr;
            while (*ptr == ' ') ++ptr;
        }
    }
    fclose(dataFile);
    
    // get surf pts
    std::unordered_map<size_t, std::vector<int>> faceHashTable;// isOnFace
    faceHashTable.reserve(eNum);
    size_t hashVal;
    //std::clock_t start = std::clock();
    for (i = 0; i < eNum; ++i)
        for (j = 0; j < 6; ++j) {
            hashVal = hashFace(hex[i][facePt[j][0]], hex[i][facePt[j][1]], hex[i][facePt[j][2]], hex[i][facePt[j][3]]);
            if (faceHashTable.find(hashVal) == faceHashTable.end())
                faceHashTable[hashVal] = { true, hex[i][facePt[j][0]], hex[i][facePt[j][1]], hex[i][facePt[j][2]], hex[i][facePt[j][3]] };
            else
                faceHashTable[hashVal][0] = false;
        }
    //std::cout << ((double)(clock() - start)) * 1000.0 / CLOCKS_PER_SEC << " ";
    std::vector<int> surf;
    std::unordered_set<int> facePtHexHashTable;
    surf.reserve(faceHashTable.size() * 4);
    for (const auto& face : faceHashTable) {
        const auto& faceInfo = face.second;
        if (faceInfo[0]) {
            for (i = 1; i < 5; ++i)
                if (facePtHexHashTable.insert(faceInfo[i]).second)
                    surf.push_back(faceInfo[i]);
        }
    }
    int sPNum = surf.size();
    std::vector<int> sPNumM3(sPNum), surfD3(sPNum), isInside(pNum, -1);//-1: inside; > -1: surf index
    std::vector<double> surfX(3 * sPNum, DBL_MAX);
#pragma omp parallel for
    for (i = 0; i < sPNum; ++i) {
        sPNumM3[i] = 3 * i;
        surfD3[i] = surf[i] / 3;
        isInside[surfD3[i]] = i;
    }

    // declare optimization vars
    int iter = 0, pNum3 = 3 * pNum;
    bool eFSJ = false, fitting;
    double sJThres = -.01;
    std::vector<int> pNumM3(pNum);
#pragma omp parallel for
    for (i = 0; i < pNum; ++i)
        pNumM3[i] = 3 * i;
    std::vector<double> totGrad(pNum3, 1.0);
    std::vector<int> chkElem(eNum);


    /*double dist1, minDist;
    double tmp[3];
    for (int i = 0; i < sPNum; ++i) {
        minDist = DBL_MAX;
        for (const auto& j : edgeFtrX) {
            ptLnDist({ x[surf[i]], x[surf[i] + 1], x[surf[i] + 2] },
                { j[0], j[1], j[2] },
                { j[3], j[4], j[5] },
                { j[6], j[7], j[8], j[9] },
                tmp, dist1);
            if (dist1 < minDist) {
                minDist = dist1;
            }
        }
        if (minDist < 0.005) std::cout << surfD3[i] << std::endl;
    }
    exit(-1);*/
    /*std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-0.005, 0.005);
#pragma omp parallel for
    for (int i = 0; i < pNum; ++i) {
        if (isInside[i] == -1) {
            x[3 * i] += dis(gen);
            x[3 * i + 1] += dis(gen);
            x[3 * i + 2] += dis(gen);
        }
    }
    dataFile = fopen(outputHexFileName, "w");
    fprintf(dataFile, "# vtk DataFile Version 2.0\n");
    fprintf(dataFile, "HexOpt\n");
    fprintf(dataFile, "ASCII\n");
    fprintf(dataFile, "DATASET UNSTRUCTURED_GRID\n");
    fprintf(dataFile, "POINTS %i double\n", pNum);
    for (i = 0; i < pNum; ++i) fprintf(dataFile, "%lf %lf %lf\n", x[pNumM3[i]], x[pNumM3[i] + 1], x[pNumM3[i] + 2]);
    fprintf(dataFile, "CELLS %i %i\n", eNum, eNum * 9);
    for (i = 0; i < eNum; ++i)
        fprintf(dataFile, "8 %i %i %i %i %i %i %i %i\n", hex[i][0] / 3, hex[i][1] / 3, hex[i][2] / 3, hex[i][3] / 3, hex[i][4] / 3, hex[i][5] / 3, hex[i][6] / 3, hex[i][7] / 3);
    fprintf(dataFile, "CELL_TYPES %i\n", eNum);
    for (i = 0; i < eNum; ++i) fprintf(dataFile, "%i\n", 12);
    fclose(dataFile);
    exit(-1);*/

    // optimization
    while (true) {
        // check finish
        if (eFSJ) {
            if (fitting) {// uplevel sJ
                //std::cout << "Converged with min SJ " << -sJThres << std::endl;
                std::fill(totGrad.begin(), totGrad.end(), 1.0);

                // write to vtk file
                dataFile = fopen(outputHexFileName, "w");
                fprintf(dataFile, "# vtk DataFile Version 2.0\n");
                fprintf(dataFile, "HexOpt\n");
                fprintf(dataFile, "ASCII\n");
                fprintf(dataFile, "DATASET UNSTRUCTURED_GRID\n");
                fprintf(dataFile, "POINTS %i double\n", pNum);
                for (i = 0; i < pNum; ++i) fprintf(dataFile, "%lf %lf %lf\n", x[pNumM3[i]], x[pNumM3[i] + 1], x[pNumM3[i] + 2]);
                fprintf(dataFile, "CELLS %i %i\n", eNum, eNum * 9);
                for (i = 0; i < eNum; ++i)
                    fprintf(dataFile, "8 %i %i %i %i %i %i %i %i\n", hex[i][0] / 3, hex[i][1] / 3, hex[i][2] / 3, hex[i][3] / 3, hex[i][4] / 3, hex[i][5] / 3, hex[i][6] / 3, hex[i][7] / 3);
                fprintf(dataFile, "CELL_TYPES %i\n", eNum);
                for (i = 0; i < eNum; ++i) fprintf(dataFile, "%i\n", 12);
                fclose(dataFile);
                lapSmth(x, adjPts, pNum, pNumM3);
                std::cout << iter << ","<<-sJThres<<std::endl;
                sJThres -= .01;
            }
        }

        // get gradient
        eFSJ = gradient(tri, triX, triEdgeX, triENum,
            edgeFtrX, hexPtType,
            x, hex, eNum, pNum3,
            surf, surfD3, sPNum, sPNumM3, surfX,
            totGrad, chkElem,
            sJThres, fitting);
        
        // print
        ++iter;
        if (iter % 1000 == 0) {
            std::cout << "Iteration " << iter;
            for (i = 0; i < (16 - (int)std::to_string(iter).size()); ++i)
                std::cout << ' ';
            std::cout << "Min SJ " << eFSJ << "\t\tFitting " << fitting << std::endl;
        }
    }
    return 0;
}