#include "LinearFiniteElement.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

int main() {
  // Read mesh from file
  std::string filename = "mesh.vtk";
  std::ifstream file(filename);
  if (!file) {
    std::cerr << "Failed to open file: " << filename << std::endl;
    return 1;
  }

  // Read element connectivity and nodes
  Eigen::Matrix<double, 3, Eigen::Dynamic> nodes; // 3 rows --> x, y, z, num_cols = numNodes
  Eigen::Matrix<Eigen::Index, 4, Eigen::Dynamic> connectivity; // 4 rows --> 4 nodes of the tetrahedron, num_cols = numElements
  int numElements;
  int numNodes;
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
  std::cout << "nodes: \n" << nodes << "\n" << std::endl;
  std::cout << "connectivity: \n" << connectivity << "\n" << std::endl;
  std::cout << "numElements: \n" << numElements << "\n" << std::endl;
  std::cout << "numNodes: \n" << numNodes << "\n" << std::endl;
  
  // prepare global matrices
  Eigen::SparseMatrix<double> stiffnessMatrix(numNodes, numNodes);
  Eigen::SparseMatrix<double> massMatrix(numNodes, numNodes);
  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(numElements);
  // for (int k = 0; k < numElements; ++k) {
  //   {
  //     for (auto i = 0; i < 4; ++i) {
  //       for (auto j = i; j < 4; ++j) {
  //         stiffnessMatrix.coeffRef(connectivity(i, k), connectivity(j, k)) = 0;
  //         stiffnessMatrix.coeffRef(connectivity(j, k), connectivity(i, k)) = 0;
  //         massMatrix.coeffRef(connectivity(i, k), connectivity(j, k)) = 0;
  //         massMatrix.coeffRef(connectivity(j, k), connectivity(i, k)) = 0;
  //       }
  //     }
  //   }
  // }
  // Create a linear finite element
  //unsigned int N = 3; // dimension of the space
  apsc::LinearFiniteElement<3> linearFiniteElement;
  using Nodes = apsc::LinearFiniteElement<3>::Nodes;
  using Indexes = apsc::LinearFiniteElement<3>::Indexes;

  Nodes localNodes;
  Indexes globalNodeNumbers;
  for (auto k = 0; k < numElements; ++k) {
    // extract element data
    for (auto i = 0; i < 4; ++i) // node numbering
    {
      globalNodeNumbers(i) = connectivity(i, k);
      for (auto j = 0; j < 3; ++j) // local node coordinates
      {
        localNodes(j, i) = nodes(j, connectivity(i, k)); // localNodes(j, i) = nodes(j, globalNodeNumbers(i));
      }
    }
    linearFiniteElement.update(localNodes);
    linearFiniteElement.updateGlobalNodeNumbers(globalNodeNumbers);
    // Compute the local stiffness matrix
    linearFiniteElement.computeLocalStiffness();
    // Compute the local mass matrix
    linearFiniteElement.computeLocalMass();
    // Add the local matrices to the global matrices
    linearFiniteElement.updateGlobalStiffnessMatrix(stiffnessMatrix);
    linearFiniteElement.updateGlobalMassMatrix(massMatrix);
  }

  // Print stiffness matrix
  std::cout << "Stiffness Matrix:" << std::endl;
  std::cout << "Rows:" << stiffnessMatrix.rows() << std::endl;
  std::cout << "Cols:" << stiffnessMatrix.cols() << std::endl;
  std::cout << stiffnessMatrix << std::endl;

  // Print mass matrix
  std::cout << "Mass Matrix:" << std::endl;
  std::cout << "Rows:" << massMatrix.rows() << std::endl;
  std::cout << "Cols:" << massMatrix.cols() << std::endl;
  std::cout << massMatrix << std::endl;

  return 0;
}

// cd /home/jammy/shared-folder/project/forStudents
// make
// ./mainFem