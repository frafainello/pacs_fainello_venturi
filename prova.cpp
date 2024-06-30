#include "LinearFiniteElement.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <Eigen/Dense>

using Values = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using Gradients = Eigen::Matrix<double, Eigen::Dynamic, 3>;

//  ================== MESH ==================
struct MeshData {
    Eigen::Matrix<double, 3, Eigen::Dynamic> nodes;              // Node coordinates (x, y, z)
    Eigen::Matrix<Eigen::Index, 4, Eigen::Dynamic> connectivity; // Element connectivity
    int numElements;
    int numNodes;
};

MeshData readMesh(const std::string& filename) {
    MeshData mesh;
    std::cout << "Entered in readMesh:" << std::endl;
    // std::cout << "initialization mesh.numNodes: " << mesh.numNodes << std::endl;
    // std::cout << "initialization mesh.numElements: \n" << mesh.numElements << std::endl;
    std::ifstream file(filename);
    // std::cout << "filename: " << filename << std::endl;
    std::string line;
    std::istringstream ss;

    if (!file.is_open()) {
        throw std::runtime_error("Failed to open file: " + filename);
    }
    int i = 0;
    // std::vector<int> offsets;
    while (std::getline(file, line)) {
        // std::cout << "line " << i << ": " << line << std::endl;
        if (line.find("CELLS") != std::string::npos) {
            std::cout << "line: " << line << std::endl;
            std::string numElementsStr = line.substr(line.find_first_of("0123456789"));
            // std::cout << "numElementsStr: " << numElementsStr << std::endl;
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
            std::cout << "line: " << line << std::endl;
            std::string numNodesStr = line.substr(line.find_first_of("0123456789"));
            ss.str(numNodesStr);
            ss.clear();
            ss >> mesh.numNodes;
            std::cout << "mesh.numNodes: " << mesh.numNodes << std::endl;
            mesh.nodes.resize(3, mesh.numNodes);

            int pointsRead = 0;
            std::getline(file, line); // Move to the line where points data starts, usually just after the line with "POINTS"
            ss.str(line);
            ss.clear();
            
            double x, y, z;
            while (pointsRead < mesh.numNodes) {
                // While there are still nodes left to read and data available in the stream
                while (ss >> x >> y >> z) {
                    mesh.nodes(0, pointsRead) = x;
                    mesh.nodes(1, pointsRead) = y;
                    mesh.nodes(2, pointsRead) = z;
                    pointsRead++;
                    if (pointsRead >= mesh.numNodes) break;  // Stop if we have read enough points
                }
                if (pointsRead < mesh.numNodes) {
                    // If we haven't read enough points and the stream is exhausted, read the next line
                    if (!std::getline(file, line)) {
                        throw std::runtime_error("Not enough points provided in file");
                    }
                    ss.str(line);
                    ss.clear();
                }
            }
            std::cout << "mesh.nodes: done! \n" << std::endl;
        }

        // } else if (line.find("POINTS") != std::string::npos) {
        //     std::cout << "line: " << line << std::endl;
        //     std::string numNodesStr = line.substr(line.find_first_of("0123456789"));
        //     ss.str(numNodesStr);
        //     ss.clear();
        //     ss >> mesh.numNodes;
        //     std::cout << "mesh.numNodes: " << mesh.numNodes << std::endl;
        //     mesh.nodes.resize(3, mesh.numNodes);

        //     for (int i = 0; i < mesh.numNodes; i++) {
        //         std::getline(file, line);
        //         ss.str(line);
        //         ss.clear();
        //         ss >> mesh.nodes(0, i) >> mesh.nodes(1, i) >> mesh.nodes(2, i);
        //     }
        //     std::cout << "mesh.nodes: done! \n" << std::endl;
        // }

        i += 1;
        // } else if (line.find("OFFSETS") != std::string::npos) {
        //     std::getline(file, line);
        //     ss.str(line);
        //     ss.clear();
        //     int offset;
        //     while (ss >> offset) {
        //         offsets.push_back(offset);
        //     }
        //     if (ss.fail() && !file.eof()) { // Handle potential multi-line offsets
        //         std::getline(file, line);
        //         ss.str(line);
        //         ss.clear();
        //         while (ss >> offset) {
        //             offsets.push_back(offset);
        //         }
        //     }
        // }
        
    }

    file.close();
    return mesh;
}

//  ================== INITIAL CONDITION ==================

Values HeatEquation(const Eigen::SparseMatrix<double>& stiffnessMatrix, const std::vector<int>& boundaryIndices) {
    Eigen::SparseMatrix<double> modifiedStiffnessMatrix = stiffnessMatrix;
    // memo: cambia std::vector con Eigen::qualcosa

    // Right-hand side vector of the equation, all ones
    Values rhs(stiffnessMatrix.rows());
    rhs.setOnes();

    // // Apply null Dirichlet boundary conditions
    // for (int idx : boundaryIndices) {
    //     modifiedStiffnessMatrix.row(idx).setZero();
    //     modifiedStiffnessMatrix.col(idx).setZero();

    //     modifiedStiffnessMatrix.coeffRef(idx, idx) = 1.0;

    //     rhs[idx] = 0.0;
    // }
    // Apply null Dirichlet boundary conditions
    for (int idx : boundaryIndices) {
        // Zero out the row
        for (Eigen::SparseMatrix<double>::InnerIterator it(modifiedStiffnessMatrix, idx); it; ++it) {
            it.valueRef() = 0.0;
        }

        // Zero out the column
        for (int k = 0; k < modifiedStiffnessMatrix.outerSize(); ++k) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(modifiedStiffnessMatrix, k); it; ++it) {
                if (it.row() == idx) {
                    it.valueRef() = 0.0;
                }
            }
        }

        // Set diagonal element to 1
        modifiedStiffnessMatrix.coeffRef(idx, idx) = 1.0;

        // Set corresponding rhs element to 0
        rhs[idx] = 0.0;
    }

    // Compute the solution
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> solver;
    Values solution = solver.compute(modifiedStiffnessMatrix).solve(rhs);

    // Return the computed values
    return solution;
}


//  ================== INCREMENTAL SOLUTION ==================

// Function to update the solution vector based on the Eikonal equation
bool updateSolution(Values& w, 
                    const Eigen::SparseMatrix<double>& stiffnessMatrix, 
                    const Eigen::SparseMatrix<double>& gradientMatrix,
                    const std::vector<int> &boundaryIndices,
                    double gamma = 1e-3,
                    double tol = 1e-6) {
    
    std::cout << "Entered in updateSolution:" << std::endl;

    std::cout << "w: \n" << w << std::endl;
    std::cout << "w.rows(): " << w.rows() << std::endl;
    std::cout << "w.cols(): " << w.cols() << std::endl;
    std::cout << "\n" << std::endl;

    // std::cout << "stiffnessMatrix: \n" << stiffnessMatrix << std::endl;
    // std::cout << "stiffnessMatrix.rows(): " << stiffnessMatrix.rows() << std::endl;
    // std::cout << "stiffnessMatrix.cols(): " << stiffnessMatrix.cols() << std::endl;
    // std::cout << "\n" << std::endl;

    std::cout << "gradientMatrix: \n" << gradientMatrix << std::endl;
    std::cout << "gradientMatrix.rows(): " << gradientMatrix.rows() << std::endl;
    std::cout << "gradientMatrix.cols(): " << gradientMatrix.cols() << std::endl;
    std::cout << "\n" << std::endl;

    Values bilinear_form = stiffnessMatrix * w;
    std::cout << "bilinear_form: \n" << bilinear_form << std::endl;
    std::cout << "bilinear_form.rows(): " << bilinear_form.rows() << std::endl;
    std::cout << "bilinear_form.cols(): " << bilinear_form.cols() << std::endl;
    std::cout << "\n" << std::endl;

    Gradients grad_w(gradientMatrix.rows(), 3);
    // std::cout << "grad_w: \n" << grad_w << std::endl;
    // std::cout << "grad_w.rows(): " << grad_w.rows() << std::endl;
    // std::cout << "grad_w.cols(): " << grad_w.cols() << std::endl;
    // std::cout << "\n" << std::endl;

    Gradients w_concatenated = w.replicate(1, 3);
    // std::cout << "w_concatenated: \n" << w_concatenated << std::endl;
    // std::cout << "w_concatenated.rows(): " << w_concatenated.rows() << std::endl;
    // std::cout << "w_concatenated.cols(): " << w_concatenated.cols() << std::endl;
    // std::cout << "\n" << std::endl;

    grad_w = gradientMatrix.cwiseProduct(w_concatenated);

    std::cout << "grad_w: \n" << grad_w << std::endl;
    std::cout << "grad_w.rows(): " << grad_w.rows() << std::endl;
    std::cout << "grad_w.cols(): " << grad_w.cols() << std::endl;
    std::cout << "\n" << std::endl;

    // Eigen::MatrixXd dense_grad_w = Eigen::MatrixXd(grad_w);
    Values norm_grad_w = grad_w.rowwise().norm();
    std::cout << "norm_grad_w: " << norm_grad_w << std::endl;
    std::cout << "norm_grad_w.rows(): " << norm_grad_w.rows() << std::endl;
    std::cout << "norm_grad_w.cols(): " << norm_grad_w.cols() << std::endl;
    std::cout << "\n" << std::endl;

    Values gamma_vec = Values::Constant(gradientMatrix.rows(), gamma);
    // std::cout << "gamma_vec: " << gamma_vec << std::endl;
    // std::cout << "gamma_vec.rows(): " << gamma_vec.rows() << std::endl;
    // std::cout << "gamma_vec.cols(): " << gamma_vec.cols() << std::endl;
    // std::cout << "\n" << std::endl;
    
    // Values coeffs = ((Values::Constant(gradientMatrix.cols(), 1.0) - norm_grad_w).array() / (norm_grad_w + gamma_vec).array()).matrix();
    Values coeffs = ((Values::Constant(gradientMatrix.rows(), 1.0) - norm_grad_w).array() / (norm_grad_w + gamma_vec).array()).matrix();
    // std::cout << "coeffs: " << coeffs << std::endl;
    // std::cout << "coeffs.rows(): " << coeffs.rows() << std::endl;
    // std::cout << "coeffs.cols(): " << coeffs.cols() << std::endl;
    // std::cout << "\n" << std::endl;

    // Compute the right-hand side for the linear problem
    Values rhs = coeffs * (stiffnessMatrix * w);
    // Update for Dirichlet BC. If != 0, set with value*TGV, where TGV=1e40
    for (int idx : boundaryIndices) {
        rhs(idx) = 0.0;
    }

    std::cout << "rhs: " << rhs << std::endl;
    std::cout << "rhs.rows(): " << rhs.rows() << std::endl;
    std::cout << "rhs.cols(): " << rhs.cols() << std::endl;
    std::cout << "\n" << std::endl;

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

    MeshData mesh = readMesh("mesh_cubo.vtk");
    
    std::cout << "numNodes: \n" << mesh.numNodes << "\n" << std::endl;
    // std::cout << "nodes: \n" << mesh.nodes << "\n" << std::endl;
    // std::cout << "connectivity: \n" << mesh.connectivity << "\n" << std::endl;
    std::cout << "numElements: \n" << mesh.numElements << "\n" << std::endl;
  
    // prepare global matrices
    Eigen::SparseMatrix<double> stiffnessMatrix(mesh.numNodes, mesh.numNodes);
    Eigen::SparseMatrix<double> massMatrix(mesh.numNodes, mesh.numNodes);
    Eigen::SparseMatrix<double> gradientMatrix(mesh.numNodes, 3);

    std::vector<Eigen::Triplet<double>> triplets;
    std::cout << "max allocable size: \n" << triplets.max_size() << "\n" << std::endl;
    triplets.reserve(mesh.numElements);

    apsc::LinearFiniteElement<3> linearFiniteElement;
    linearFiniteElement.initializeRefGradients();

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

        // Compute local nodes and update global node numbers
        linearFiniteElement.update(localNodes);
        linearFiniteElement.updateGlobalNodeNumbers(globalNodeNumbers);
        
        // Compute the local stiffness matrix and update the global stiffness matrix
        linearFiniteElement.computeLocalStiffness();
        linearFiniteElement.updateGlobalStiffnessMatrix(stiffnessMatrix);
        
        // Compute the local mass matrix and update the global mass matrix
        linearFiniteElement.computeLocalMass();
        linearFiniteElement.updateGlobalMassMatrix(massMatrix);
        
        // Compute the local gradient matrix and update the global gradient matrix
        linearFiniteElement.computeLocalGradient();
        linearFiniteElement.updateGlobalGradientMatrix(gradientMatrix);
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

    // Print gradient matrix
    std::cout << "Gradient Matrix:" << std::endl;
    std::cout << "Rows:" << gradientMatrix.rows() << std::endl;
    std::cout << "Cols:" << gradientMatrix.cols() << std::endl;
    std::cout << gradientMatrix << std::endl;
    // gradientMatriix = gradientMatrix.transpose();


    // Solve the Eikonal Equation
    // Solve the Heat Equation for initial conditions
    Values forcingTerm = Values::Constant(mesh.numNodes, 1.0);
    std::vector<int> boundaryIndices = {0, 2};
    Values initial_conditions = HeatEquation(stiffnessMatrix, boundaryIndices);
    
    // Values initial_conditions = Values::Constant(mesh.numNodes, 5.0);
    std::cout << "initial_conditions:" << std::endl;
    std::cout << "Rows:" << initial_conditions.rows() << std::endl;
    std::cout << "Cols:" << initial_conditions.cols() << std::endl;
    std::cout << initial_conditions << std::endl;

    Values w = initial_conditions;
    
    // Impose Dirichlet BC on stiffness matrix
    linearFiniteElement.updateMatrixWithDirichletBoundary(stiffnessMatrix, boundaryIndices);

    bool converged = false;
    int maxIterations = 1000;

    for (int iter = 0; iter < maxIterations && !converged; ++iter) {
        std::cout << "-------------- Iteration " << iter + 1 << "--------------" << std::endl;
        converged = updateSolution(w, stiffnessMatrix, gradientMatrix, boundaryIndices);
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



