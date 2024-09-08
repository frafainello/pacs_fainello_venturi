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

/**
 * @brief A class to represent mesh data in finite element analysis.
 * 
 * @tparam PHDIM Physical dimension of the mesh.
 * @tparam INTRINSIC_DIM Intrinsic dimension of the mesh (default: PHDIM).
 */

template<std::size_t PHDIM, std::size_t INTRINSIC_DIM=PHDIM>
class MeshData {

public:

    using Traits = Eikonal::Eikonal_traits<PHDIM, INTRINSIC_DIM>;
    using Index = typename Traits::Index;
    using Indexes = typename Traits::Indexes;
    using Nodes = typename Traits::Nodes;
    using Elements = typename Traits::Elements;
    using Values = typename Traits::Values;

    /**
     * @brief Default constructor initializing the mesh data.
     */
    MeshData() : numElements(0), numNodes(0) {}

    // Getters
    /**
     * @brief Gets the nodes of the mesh.
     * @return The nodes matrix.
     */
    Nodes getNodes() const {
        return nodes;
    }

    /**
     * @brief Gets the element connectivity of the mesh.
     * @return Reference to the element connectivity matrix.
     */
    const Elements& getConnectivity() const {
        return connectivity;
    }

    /**
     * @brief Gets the number of elements in the mesh.
     * @return The number of elements.
     */
    int getNumElements() const {
        return numElements;
    }

    /**
     * @brief Gets the number of nodes in the mesh.
     * @return The number of nodes.
     */
    int getNumNodes() const {
        return numNodes;
    }

    /**
     * @brief Gets the indices of the boundary nodes.
     * @return Reference to the vector of boundary node indices.
     */
    const Indexes& getBoundaryNodes() const {
        return boundaryNodes;
    }

    /**
     * @brief Gets the gradient coefficient matrices.
     * @return A constant reference to the vector of gradient coefficient matrices.
     */
    const std::vector<Eigen::Matrix<double, PHDIM, PHDIM>>& getGradientCoeff() const {
        return gradientCoeff;
    }

    /**
     * @brief Gets the local stiffness matrices.
     * @return A constant reference to the vector of local stiffness matrices.
     */
    const std::vector<Eigen::Matrix<double, PHDIM+1, PHDIM+1>>& getLocalStiffnessMatrices() const {
        return localStiffnessMatrices;
    }

    /**
     * @brief Gets the local reaction matrices.
     * @return A constant reference to the vector of local reaction matrices.
     */
    const std::vector<std::vector<std::vector<Eigen::Matrix<double, PHDIM, 1>>>>& getLocalReactionMatrices() const {
        return localReactionMatrices;
    }


    /**
     * @brief Reads mesh data from a file.
     * @param filename The name of the file containing mesh data.
     */
    void readMesh(const std::string& filename) {
        // Read mesh from file
        std::ifstream file(filename);
        if (!file) {
          std::cerr << "Failed to open file: " << filename << std::endl;
          return;
        }

        // Read element connectivity and nodes
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

        // Reserve memory for vectors
        gradientCoeff.reserve(numElements);
        localStiffnessMatrices.reserve(numElements);
        localReactionMatrices.reserve(numElements);
    }

    /**
     * @brief Converts a VTK file format to a different format.
     * @param inputFilePath The path to the input VTK file.
     * @param outputFilePath The path to the output file.
     */
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

        inputFile.close();

        outputFile << "# vtk DataFile Version 3.0" << std::endl;
        outputFile << "description of the mesh" << std::endl;
        outputFile << "ASCII" << std::endl;
        outputFile << "DATASET UNSTRUCTURED_GRID" << std::endl;

        outputFile << "POINTS " << numPoints << " float" << std::endl;
        for (const auto& point : points) {
            outputFile << point[0] << " " << point[1] << " " << point[2] << std::endl;
        }

        outputFile << "CELLS " << numCells-1 << " " << numDuplicatedNodes+numCellsTypes << std::endl;
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

    /**
     * @brief Fills global variables used in finite element computations.
     * @param stiffnessMatrix The global stiffness matrix to be filled.
     * @param reactionMatrix The global reaction matrix to be filled.
     * @param globalIntegrals The global integrals to be computed and filled.
     */
    void fillGlobalVariables(Eigen::SparseMatrix<double>& stiffnessMatrix, 
                        std::vector<std::vector<Eigen::Matrix<double, PHDIM, 1>>>& reactionMatrix,
                        Values& globalIntegrals) {
                        
        Nodes localNodes(PHDIM, PHDIM+1);
        Indexes globalNodeNumbers(PHDIM+1);

        apsc::LinearFiniteElement<PHDIM> linearFiniteElement;

        for (int k = 0; k < getNumElements(); ++k) {
            for (int i = 0; i < PHDIM+1; ++i) {  // Node numbering
                globalNodeNumbers[i] = getConnectivity()(i, k);
                for (int j = 0; j < PHDIM; ++j) {  // Local node coordinates
                    localNodes(j, i) = getNodes()(j, getConnectivity()(i, k));
                }
            }
            
            // Compute local nodes and update global node numbers
            linearFiniteElement.update(localNodes);
            linearFiniteElement.updateGlobalNodeNumbers(globalNodeNumbers);

            // Compute integrals and update global matrices
            linearFiniteElement.computeLocalIntegral();
            linearFiniteElement.updateGlobalIntegrals(globalIntegrals);
            // std::cout << "Global integrals: " << std::endl << globalIntegrals << std::endl;

            // Compute stiffness and update global stiffness matrix
            linearFiniteElement.computeLocalStiffness();
            linearFiniteElement.updateGlobalStiffnessMatrix(stiffnessMatrix);

            // Compute the local reaction matrix and update the global reaction matrix
            linearFiniteElement.computeLocalReaction();
            linearFiniteElement.updateGlobalReactionMatrix(reactionMatrix);

            // Compute gradient coefficients
            gradientCoeff.push_back(linearFiniteElement.computeGradientCoeff());

            // Compute local matrices
            localStiffnessMatrices.push_back(linearFiniteElement.computeLocalStiffness());
            localReactionMatrices.push_back(linearFiniteElement.computeLocalReaction());
        }

        // Impose Dirichlet boundary conditions on the stiffness matrix
        linearFiniteElement.updateMatrixWithDirichletBoundary(stiffnessMatrix, getBoundaryNodes());
    }


    /**
     * @brief Adds a scalar field to the mesh data and writes it to a VTK file.
     * @param values The scalar field values to be added.
     * @param inputFilePath The input VTK file path.
     * @param outputFilePath The output VTK file path.
     * @throws std::invalid_argument if the size of values does not match the number of nodes.
     */
    void addScalarField(const Eigen::Matrix<double, Eigen::Dynamic, 1>& values, 
                        const std::string& inputFilePath,
                        const std::string& outputFilePath) {

        if (values.size() != numNodes) {
            throw std::invalid_argument("The size of values must match the number of nodes in the mesh.");
        }

        // std::string outputFilePath = inputFilePath.substr(0, inputFilePath.find_last_of('.')) + "_" + bc + "_" + ic + iterativeMethod + ".vtk";

        // Open the input file and create the output file
        std::ifstream inputFile(inputFilePath);
        std::ofstream outputFile(outputFilePath);

        if (!inputFile.is_open()) {
            std::cerr << "Error: Could not open input file " << inputFilePath << std::endl;
            return;
        }

        if (!outputFile.is_open()) {
            outputFile.open(outputFilePath);
            std::cerr << "Error: Could not open output file " << outputFilePath << std::endl;
            return;
        }

        // Copy the contents of the input file to the output file
        std::string line;
        while (std::getline(inputFile, line)) {
            outputFile << line << std::endl;
        }

        // Append the scalar field data to the new file
        outputFile << "POINT_DATA " << numNodes << std::endl;
        outputFile << "SCALARS sol float" << std::endl;
        outputFile << "LOOKUP_TABLE default" << std::endl;
        for (int i = 0; i < values.size(); ++i) {
            outputFile << values(i) << std::endl;
        }

        inputFile.close();
        outputFile.close();
    }

    /**
     * @brief Updates the boundary nodes based on a given boundary condition.
     * @param boundaryCondition The boundary condition to be applied.
     */
    void updateBoundaryNodes(const int& boundaryCondition) {
        std::map<std::set<int>, int> faceCount;
        for (int i = 0; i < numElements; ++i) {
            // Each element has 4 faces
            for (int j = 0; j < 4; ++j) {
                std::set<int> face;
                for (int k = 0; k < 4; ++k) {
                    if (k != j) {
                        face.insert(connectivity(k, i)); // Assuming connectivity is a 2D array or similar structure
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

        boundaryNodes.clear();
        switch (boundaryCondition) {
            case 1:
                // Include only nodes with z=0
                for (int node : boundaryNodesSet) {
                    if (nodes(2, node) == 0) { // Assuming node position z is at index 2 in nodes
                        boundaryNodes.push_back(node);
                    }
                }
                break;
            case 2:
                // Include only the center node of the face with z=0
                {
                    std::vector<int> z0Nodes;
                    for (int node : boundaryNodesSet) {
                        if (nodes(2, node) == 0) { // Check for z=0 nodes
                            z0Nodes.push_back(node);
                        }
                    }

                    if (!z0Nodes.empty()) {
                        // Calculate the geometric center by averaging coordinates
                        double centerX = 0.0, centerY = 0.0;
                        for (int node : z0Nodes) {
                            centerX += nodes(0, node); // X coordinate
                            centerY += nodes(1, node); // Y coordinate
                        }
                        centerX /= z0Nodes.size();
                        centerY /= z0Nodes.size();

                        // Find the node closest to the calculated center
                        int centerNode = z0Nodes[0];
                        double minDistance = std::numeric_limits<double>::max();
                        for (int node : z0Nodes) {
                            double distance = std::sqrt(std::pow(nodes(0, node) - centerX, 2) + std::pow(nodes(1, node) - centerY, 2));
                            if (distance < minDistance) {
                                minDistance = distance;
                                centerNode = node;
                            }
                        }
                        boundaryNodes.push_back(centerNode);
                    }
                }
                break;
            case 3:
                // Include only the vertex at the origin (0, 0, 0)
                for (int node : boundaryNodesSet) {
                    if (nodes(0, node) == 0 && nodes(1, node) == 0 && nodes(2, node) == 0) {
                        boundaryNodes.push_back(node);
                        break; // No need to search further, there is only one origin
                    }
                }
                break;
            case 4:
                // Include all boundary nodes
                boundaryNodes.assign(boundaryNodesSet.begin(), boundaryNodesSet.end());
                break;
            default:
                std::cerr << "Invalid option. Please run the program again with a valid choice.\n";
                return; // Exit with an error code
            }
    }

    /**
     * @brief Stores gradient coefficient matrices for each element.
     */
    std::vector<Eigen::Matrix<double, PHDIM, PHDIM>> gradientCoeff;
    /**
     * @brief Stores local stiffness matrices for each element.
     */
    std::vector<Eigen::Matrix<double, PHDIM+1, PHDIM+1>> localStiffnessMatrices;
    /**
     * @brief Stores local reaction matrices for each element.
     */
    std::vector<std::vector<std::vector<Eigen::Matrix<double, PHDIM, 1>>>> localReactionMatrices;

private:
    /**
     * @brief Coordinates of nodes in the mesh.
     */
    Nodes nodes;              // Node coordinates (x, y, z)
    /**
     * @brief Connectivity of the elements in the mesh.
     */
    Elements connectivity;    // Element connectivity
    /**
     * @brief Number of elements in the mesh.
     */
    int numElements;
    /**
     * @brief Number of nodes in the mesh.
     */
    int numNodes;
    /**
     * @brief Indices of the boundary nodes.
     */
    Indexes boundaryNodes;    // Boundary nodes indices
    };

#endif // MESHDATA_HPP
