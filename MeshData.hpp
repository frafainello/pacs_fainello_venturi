// MeshData.hpp
#ifndef MESHDATA_HPP
#define MESHDATA_HPP

#include "Eikonal_traits.hpp"
#include <Eigen/Dense>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <iostream>

template<std::size_t PHDIM, std::size_t INTRINSIC_DIM=PHDIM>
class MeshData {

public:
    using Traits = Eikonal::Eikonal_traits<PHDIM, INTRINSIC_DIM>;
    using Index = typename Traits::Index;
    using Indexes = typename Traits::Indexes;
    using Nodes = typename Traits::Nodes;
    using Elements = typename Traits::Elements;


    MeshData() : numElements(0), numNodes(0) {}

    // Getters
    const Nodes& getNodes() const {
        return nodes;
    }

    const Elements& getConnectivity() const {
        return connectivity;
    }

    int getNumElements() const {
        return numElements;
    }

    int getNumNodes() const {
        return numNodes;
    }

    const Indexes& getBoundaryNodes() const {
        return boundaryNodes;
    }

    // // Setters
    // void setNodes(const Nodes& nodes) {
    //     this->nodes = nodes;
    // }

    // void setConnectivity(const Elements& connectivity) {
    //     this->connectivity = connectivity;
    // }

    // void setNumElements(int numElements) {
    //     this->numElements = numElements;
    // }

    // void setNumNodes(int numNodes) {
    //     this->numNodes = numNodes;
    // }

    void readMesh(const std::string& filename) {
        // Read mesh from file
        std::ifstream file(filename);
        if (!file) {
          std::cerr << "Failed to open file: " << filename << std::endl;
          return;
        }

        // Read element connectivity and nodes
        // Nodes nodes; // 3 rows --> x, y, z, num_cols = numNodes
        // Elements connectivity; // 4 rows --> 4 nodes of the tetrahedron, num_cols = numElements
        // int numElements;
        // int numNodes;
        std::string line;
        std::istringstream ss;
        while (std::getline(file, line)) {
          if (line.find("CELLS") != std::string::npos) {
            // extract the substring starting from the first digit
            std::string numElementsStr =
                line.substr(line.find_first_of("0123456789"));
            // put the string in a stringstream
            ss.str(numElementsStr);
            ss.seekg(0, std::ios::beg);
            // read the number of elements
            ss >> numElements;
            connectivity.resize(4, numElements);
            for (int i = 0; i < numElements; i++) {
              int dummy;
              std::getline(file, line);
              ss.str(line);
              ss.seekg(0, std::ios::beg);
              ss >> dummy;
              for (int j = 0; j < 4; j++) {
                ss >> connectivity(j, i);
              }
            }
          } else if (line.find("POINTS") != std::string::npos) {
            std::string numNodesStr = line.substr(line.find_first_of("0123456789"));
            ss.str(numNodesStr);
            ss.seekg(0, std::ios::beg);
            ss >> numNodes;
            nodes.resize(3, numNodes);
            for (int i = 0; i < numNodes; i++) {
              std::getline(file, line);
              ss.str(line);
              ss.seekg(0, std::ios::beg);
              ss >> nodes(0, i) >> nodes(1, i) >> nodes(2, i);
            }
          }
        }
        // std::cout << "nodes: \n" << nodes << "\n" << std::endl;
        // std::cout << "connectivity: \n" << connectivity << "\n" << std::endl;
        // std::cout << "numElements: \n" << numElements << "\n" << std::endl;
        // std::cout << "numNodes: \n" << numNodes << "\n" << std::endl;
    }

    void convertVTK(const std::string& inputFilePath, std::string outputFilePath) {
      std::cout << "Entered in convertVTK:" << std::endl;
      std::ifstream inputFile(inputFilePath);
      if (!inputFile.is_open()) {
          std::cerr << "Error: Could not open input file " << inputFilePath << std::endl;
          return;
      }

      std::ofstream outputFile(outputFilePath);
      if (!outputFile.is_open()) {
          // outputFilePath = createOutputFileName(inputFilePath);
          outputFile.open(outputFilePath);
          if (!outputFile.is_open()) {
              std::cerr << "Error: Could not create output file " << outputFilePath << std::endl;
              return;
          }
      }

      std::string line;
      std::vector<std::vector<double>> points;
      Indexes cells;
      Indexes cellTypes;
      Indexes connectivity;

      bool readingPoints = false;
      bool readingCells = false;
      bool readingCellTypes = false;
      bool readingConnectivity = false;

      int numPoints = 0;
      int numCells = 0;
      int numDuplicatedNodes = 0;
      int numCellsTypes = 0;

      while (std::getline(inputFile, line)) {
          if (line.find("POINTS") != std::string::npos) {
              readingPoints = true;
              readingCells = false;
              readingCellTypes = false;
              readingConnectivity = false;
              // outputFile << line << std::endl;
              numPoints = std::stoi(line.substr(7)); // Assumes format "POINTS <numPoints> <dataType>"
              points.reserve(numPoints);

          } else if (line.find("CELLS") != std::string::npos) {
              readingPoints = false;
              readingCells = true;
              readingCellTypes = false;
              readingConnectivity = false;
              // outputFile << line << std::endl;
              std::istringstream ss(line);
              std::string dummy;
              ss >> dummy >> numCells >> numDuplicatedNodes;
              cells.reserve(numCells);

          } else if (line.find("CELL_TYPES") != std::string::npos) {
              readingPoints = false;
              readingCells = false;
              readingCellTypes = true;
              readingConnectivity = false;
              // outputFile << line << std::endl;
              std::istringstream ss(line);
              std::string dummy;
              ss >> dummy >> numCellsTypes;
              cellTypes.reserve(numCellsTypes);

          } else if (line.find("CONNECTIVITY") != std::string::npos) {
              readingPoints = false;
              readingCells = false;
              readingCellTypes = false;
              readingConnectivity = true;
            //   outputFile << line << std::endl;
              connectivity.reserve(numDuplicatedNodes);

          } else if (isalpha(line[0])) {
            // Skip lines starting with non-number literals that are not among the keywords
            continue;

          } else {
              std::istringstream ss(line);
              if (readingPoints) {
                  double x, y, z;
                  while (ss >> x >> y >> z) {
                      points.push_back({x, y, z});
                      // std::cout << "points: " << points.size() << std::endl;
                  // points.push_back(line);
                  // std::cout << "points: " << points.size() << std::endl;
                  }
              } else if (readingCells) {
                  int n;
                  while (ss >> n) {
                      cells.push_back(n);
                      // std::cout << "cells: " << cells.size() << std::endl;
                  // cells.push_back(line);
                  // std::cout << "cells: " << cells.size() << std::endl;
                  }
              } else if (readingConnectivity) {
                  int n;
                  while (ss >> n) {
                      connectivity.push_back(n);
                      // std::cout << "cellTypes: " << cellTypes.size() << std::endl;
                  // connectivity.push_back(line);
                  // std::cout << "cellTypes: " << cellTypes.size() << std::endl;
                  }
              } else if (readingCellTypes) {
                  int n;
                  while (ss >> n) {
                      cellTypes.push_back(n);
                      // std::cout << "cellTypes: " << cellTypes.size() << std::endl;
                  // cellTypes.push_back(line);
                  // std::cout << "cellTypes: " << cellTypes.size() << std::endl;
                  }
              }
            }
          }

      std::cout << "numPoints: " << numPoints << std::endl;
      std::cout << "numCells: " << numCells << std::endl;
      std::cout << "numDuplicatedNodes: " << numDuplicatedNodes << std::endl;
      std::cout << "numCellsTypes: " << numCellsTypes << std::endl;

      std::cout << "points.size(): " << points.size() << std::endl;
      std::cout << "cells.size(): " << cells.size() << std::endl;
      std::cout << "cellTypes.size(): " << cellTypes.size() << std::endl;
      std::cout << "connectivity.size(): " << connectivity.size() << std::endl;

      inputFile.close();

      outputFile << "# vtk DataFile Version 3.0" << std::endl;
      outputFile << "description of the mesh" << std::endl;
      outputFile << "ASCII" << std::endl;
      outputFile << "DATASET UNSTRUCTURED_GRID" << std::endl;

      outputFile << "POINTS " << numPoints << " float" << std::endl;
      for (const auto& point : points) {
        outputFile << point[0] << " " << point[1] << " " << point[2] << std::endl;
      }

      outputFile << "CELLS " << numCells-1 << " " << numDuplicatedNodes << std::endl;
      for (int i = 0; i < cells.size()-1; ++i) {

          int numCurrNodes = cells[i+1] - cells[i];
          outputFile << numCurrNodes << " ";
          for (int j = 0; j < numCurrNodes; ++j) {
              outputFile << connectivity[cells[i]+j] << " ";
          }
          outputFile << std::endl;
      }

      outputFile << "CELL_TYPES " << numCells-1 << std::endl;
      for (const auto& cellType : cellTypes) {
          outputFile << cellType << std::endl;
      }

      outputFile.close();
  }


    void updateBoundaryNodes() {
        std::map<std::set<int>, int> faceCount;
        for (int i = 0; i < numElements; ++i) {
            // Each element has 4 faces
            for (int j = 0; j < 4; ++j) {
                std::set<int> face;
                for (int k = 0; k < 4; ++k) {
                    if (k != j) {
                        face.insert(connectivity(k, i));
                    }
                }
                faceCount[face]++;
            }
        }

        std::set<int> boundaryNodesSet;
        for (const auto& entry : faceCount) {
            if (entry.second == 1) { // Boundary face
                for (int node : entry.first) {
                    boundaryNodesSet.insert(node);
                }
            }
        }

        // TAKE ALL BOUNDARY NODES
        // boundaryNodes.assign(boundaryNodesSet.begin(), boundaryNodesSet.end());

        // OR

        // TAKE ONLY ONE FACE (z=0)
        boundaryNodes.clear();
        for (int node : boundaryNodesSet) {
            if (nodes(2, node) == 0) {
                boundaryNodes.push_back(node);
            }
        }
    }

private:
    Nodes nodes;              // Node coordinates (x, y, z)
    Elements connectivity;    // Element connectivity
    int numElements;
    int numNodes;
    Indexes boundaryNodes;    // Boundary nodes indices
};

#endif // MESHDATA_HPP
