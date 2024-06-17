#include "LinearFiniteElement.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <Eigen/Dense>

using Values = Eigen::Matrix<double, Eigen::Dynamic, 1>;

//  ================== MESH ==================
struct MeshData {
    Eigen::Matrix<double, 3, Eigen::Dynamic> nodes;              // Node coordinates (x, y, z)
    Eigen::Matrix<Eigen::Index, 4, Eigen::Dynamic> connectivity; // Element connectivity
    int numElements;
    int numNodes;
};

MeshData readMesh(const std::string& filename) {
    MeshData mesh;
    std::ifstream file(filename);
    std::string line;
    std::istringstream ss;

    if (!file.is_open()) {
        throw std::runtime_error("Failed to open file: " + filename);
    }

    while (std::getline(file, line)) {
        if (line.find("CELLS") != std::string::npos) {
            std::string numElementsStr = line.substr(line.find_first_of("0123456789"));
            ss.str(numElementsStr);
            ss.clear();
            ss >> mesh.numElements;
            mesh.connectivity.resize(4, mesh.numElements);

            for (int i = 0; i < mesh.numElements; i++) {
                std::getline(file, line);
                ss.str(line);
                ss.clear();
                int dummy;  // Dummy to hold the number of nodes per element (should be 4 for tetrahedra)
                ss >> dummy;
                for (int j = 0; j < 4; j++) {
                    ss >> mesh.connectivity(j, i);
                }
            }
        } else if (line.find("POINTS") != std::string::npos) {
            std::string numNodesStr = line.substr(line.find_first_of("0123456789"));
            ss.str(numNodesStr);
            ss.clear();
            ss >> mesh.numNodes;
            mesh.nodes.resize(3, mesh.numNodes);

            for (int i = 0; i < mesh.numNodes; i++) {
                std::getline(file, line);
                ss.str(line);
                ss.clear();
                ss >> mesh.nodes(0, i) >> mesh.nodes(1, i) >> mesh.nodes(2, i);
            }
        }
    }

    file.close();
    return mesh;
}

//  ================== INITIAL CONDITION ==================
Values solveHeatEquation(const Eigen::SparseMatrix<double>& stiffnessMatrix, 
                        const Values& forcingTerm) {
                        // const Mesh& mesh, 
                        // const std::vector<int>& boundaryIndices, 
                        // const Values& boundaryValues,
                        
    // std::cout << "entered in solveHeatEquation:" << std::endl;
    
    // // Print forcingTerm
    // std::cout << "forcingTerm:" << std::endl;
    // std::cout << "Rows:" << forcingTerm.rows() << std::endl;
    // std::cout << "Cols:" << forcingTerm.cols() << std::endl;
    // std::cout << forcingTerm << std::endl;
    
    Values w = Values::Zero(forcingTerm.size());
    
    // // Print w
    // std::cout << "w:" << std::endl;
    // std::cout << "Rows:" << w.rows() << std::endl;
    // std::cout << "Cols:" << w.cols() << std::endl;
    // std::cout << w << std::endl;

    // // Applying boundary conditions: MODIFY ACCORDINGLY THE STIFFNESS MATRIX
    // // Set the rows and columns for boundary nodes to zero and set the diagonal to 1
    // Eigen::SparseMatrix<double> modifiedStiffnessMatrix = stiffnessMatrix;
    // for (int idx : boundaryIndices) {
    //     for (int k = 0; k < modifiedStiffnessMatrix.outerSize(); ++k) {
    //         for (Eigen::SparseMatrix<double>::InnerIterator it(modifiedStiffnessMatrix, k); it; ++it) {
    //             if (it.row() == idx || it.col() == idx) {
    //                 it.valueRef() = (it.row() == idx && it.col() == idx) ? 1.0 : 0.0;
    //             }
    //         }
    //     }
    //     w(idx) = boundaryValues(idx);
    // }

    // Set up the solver using Eigen's Conjugate Gradient solver
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> solver;
    solver.compute(stiffnessMatrix);
    
    if (solver.info() != Eigen::Success) {
        std::cerr << "Decomposition failed!" << std::endl;
        throw std::runtime_error("Matrix decomposition failed in solveHeatEquation.");
    }

    // Solve the linear system
    w = solver.solve(forcingTerm);

    if (solver.info() != Eigen::Success) {
        std::cerr << "Solving failed!" << std::endl;
        throw std::runtime_error("Solver failed in solveHeatEquation.");
    }

    return w;
}

//  ================== INCREMENTAL SOLUTION ==================

// Function to update the solution vector based on the Eikonal equation
bool updateSolution(Values& w, 
                    const Eigen::SparseMatrix<double>& stiffnessMatrix, 
                    double gamma = 1e-3,
                    double tol = 1e-6) {
    
    std::cout << "entered in updateSolution:" << std::endl;

    // Values bilinear_form = stiffnessMatrix * w;

    // // Print stiffnessMatrix * w
    // std::cout << "stiffnessMatrix * w:" << std::endl;
    // std::cout << "Rows:" << bilinear_form.rows() << std::endl;
    // std::cout << "Cols:" << bilinear_form.cols() << std::endl;
    // std::cout << bilinear_form << std::endl;

    double norm_grad_w = 1.0; // Compute the norm of the gradient of w at each node
    std::cout << "norm_grad_w: " << norm_grad_w << std::endl;

    double coeff = (1.0 - norm_grad_w) / (norm_grad_w + gamma);

    // Compute the right-hand side for the linear problem
    Values rhs = coeff * (stiffnessMatrix * w);
    // // Print rhs
    // std::cout << "rhs:" << std::endl;
    // std::cout << "Rows:" << rhs.rows() << std::endl;
    // std::cout << "Cols:" << rhs.cols() << std::endl;
    // std::cout << rhs << std::endl;

    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> solver;
    solver.compute(stiffnessMatrix);

    if(solver.info() != Eigen::Success) {
        throw std::runtime_error("Decomposition failed!");
    }

    // Solve the lienar system
    Values z = solver.solve(rhs);
    if(solver.info() != Eigen::Success) {
        throw std::runtime_error("Solving failed!");
    }

    // Update the solution
    w += z;
    
    return (z.norm() < tol);

}

//  ================== MAIN ==================
int main() {
  
    using Nodes = apsc::LinearFiniteElement<3>::Nodes;
    using Indexes = apsc::LinearFiniteElement<3>::Indexes;
    using Values = Eigen::Matrix<double, Eigen::Dynamic, 1>; //!< The array storing values at the the nodes
    using AnisotropyM = Eigen::Matrix<double, 3, 3>; //!< The array storing the anisotropy matrix M 

//   try {
//         MeshData mesh = readMesh("path/to/your/meshfile.mesh");
//         std::cout << "Nodes: \n" << mesh.nodes << "\n" << std::endl;
//         std::cout << "Connectivity: \n" << mesh.connectivity << "\n" << std::endl;
//         std::cout << "Number of Elements: " << mesh.numElements << "\n" << std::endl;
//         std::cout << "Number of Nodes: " << mesh.numNodes << "\n" << std::endl;
//     } catch (const std::exception& e) {
//         std::cerr << "Error: " << e.what() << std::endl;
//         return -1;
//     }

    MeshData mesh = readMesh("mesh.vtk");
    
    std::cout << "nodes: \n" << mesh.nodes << "\n" << std::endl;
    std::cout << "connectivity: \n" << mesh.connectivity << "\n" << std::endl;
    std::cout << "numElements: \n" << mesh.numElements << "\n" << std::endl;
    std::cout << "numNodes: \n" << mesh.numNodes << "\n" << std::endl;
  
    // prepare global matrices
    Eigen::SparseMatrix<double> stiffnessMatrix(mesh.numNodes, mesh.numNodes);
    Eigen::SparseMatrix<double> massMatrix(mesh.numNodes, mesh.numNodes);
    
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(mesh.numElements);

    apsc::LinearFiniteElement<3> linearFiniteElement;

    Nodes localNodes;
    Indexes globalNodeNumbers;
    for (auto k = 0; k < mesh.numElements; ++k) {
        // extract element data
        for (auto i = 0; i < 4; ++i) // node numbering
        {
        globalNodeNumbers(i) = mesh.connectivity(i, k);
        for (auto j = 0; j < 3; ++j) // local node coordinates
        {
            localNodes(j, i) = mesh.nodes(j, mesh.connectivity(i, k)); // localNodes(j, i) = nodes(j, globalNodeNumbers(i));
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

    // Solve the Eikonal Equation
    Values forcingTerm = Values::Constant(mesh.numNodes, -1.0);
    Values initial_conditions = solveHeatEquation(stiffnessMatrix, forcingTerm);
    // Values initial_conditions = Values::Constant(mesh.numNodes, 5.0);
    std::cout << "initial_conditions:" << std::endl;
    std::cout << "Rows:" << initial_conditions.rows() << std::endl;
    std::cout << "Cols:" << initial_conditions.cols() << std::endl;
    std::cout << initial_conditions << std::endl;

    Values w = initial_conditions;

    bool converged = false;
    int maxIterations = 1000;

    for (int iter = 0; iter < maxIterations && !converged; ++iter) {
        converged = updateSolution(w, stiffnessMatrix);
        if (converged) {
            std::cout << "Solution converged after " << iter + 1 << " iterations." << std::endl;
        }
    }

    if (!converged) {
        std::cout << "Solution did not converge within the maximum number of iterations." << std::endl;
    }

    return 0;
}

// cd /home/jammy/shared-folder/project/forStudents
// make
// ./mainFem



