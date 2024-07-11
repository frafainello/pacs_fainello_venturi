#include "Eikonal_traits.hpp"
#include "LinearFiniteElement.hpp"
#include "MeshData.hpp"
#include "InitialConditions.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <Eigen/Dense>

using Values = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using Gradients = Eigen::Matrix<double, Eigen::Dynamic, 3>;

//  ================== MESH ==================

//  ================== INITIAL CONDITION ==================

// Values HeatEquation(const Eigen::SparseMatrix<double>& stiffnessMatrix, const std::vector<int>& boundaryIndices) {
//     Eigen::SparseMatrix<double> modifiedStiffnessMatrix = stiffnessMatrix;
//     // memo: cambia std::vector con Eigen::qualcosa

//     // Right-hand side vector of the equation, all ones
//     Values rhs(stiffnessMatrix.rows());
//     rhs.setOnes();

//     // // Apply null Dirichlet boundary conditions
//     // for (int idx : boundaryIndices) {
//     //     modifiedStiffnessMatrix.row(idx).setZero();
//     //     modifiedStiffnessMatrix.col(idx).setZero();

//     //     modifiedStiffnessMatrix.coeffRef(idx, idx) = 1.0;

//     //     rhs[idx] = 0.0;
//     // }
//     // Apply null Dirichlet boundary conditions
//     for (int idx : boundaryIndices) {
//         // Zero out the row
//         for (Eigen::SparseMatrix<double>::InnerIterator it(modifiedStiffnessMatrix, idx); it; ++it) {
//             it.valueRef() = 0.0;
//         }

//         // Zero out the column
//         for (int k = 0; k < modifiedStiffnessMatrix.outerSize(); ++k) {
//             for (Eigen::SparseMatrix<double>::InnerIterator it(modifiedStiffnessMatrix, k); it; ++it) {
//                 if (it.row() == idx) {
//                     it.valueRef() = 0.0;
//                 }
//             }
//         }

//         // Set diagonal element to 1
//         modifiedStiffnessMatrix.coeffRef(idx, idx) = 1.0;

//         // Set corresponding rhs element to 0
//         rhs[idx] = 0.0;
//     }

//     // Compute the solution
//     Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> solver;
//     Values solution = solver.compute(modifiedStiffnessMatrix).solve(rhs);

//     // Return the computed values
//     return solution;
// }


//  ================== INCREMENTAL SOLUTION ==================

// Function to update the solution vector based on the Eikonal equation
bool updateSolution(Values& w, 
                    const Eigen::SparseMatrix<double>& stiffnessMatrix, 
                    const Eigen::SparseMatrix<double>& gradientMatrix,
                    const std::vector<long int> &boundaryIndices,
                    double gamma = 1e-3,
                    double tol = 1e-6) {

    Values bilinear_form = stiffnessMatrix * w;

    Gradients grad_w(gradientMatrix.rows(), 3);

    Gradients w_concatenated = w.replicate(1, 3);

    grad_w = gradientMatrix.cwiseProduct(w_concatenated);

    // Eigen::MatrixXd dense_grad_w = Eigen::MatrixXd(grad_w);
    Values norm_grad_w = grad_w.rowwise().norm();

    Values gamma_vec = Values::Constant(gradientMatrix.rows(), gamma);
    
    // Values coeffs = ((Values::Constant(gradientMatrix.cols(), 1.0) - norm_grad_w).array() / (norm_grad_w + gamma_vec).array()).matrix();
    Values coeffs = ((Values::Constant(gradientMatrix.rows(), 1.0) - norm_grad_w).array() / (norm_grad_w + gamma_vec).array()).matrix();

    // Compute the right-hand side for the linear problem
    Values rhs = coeffs * (stiffnessMatrix * w);
    // Update for Dirichlet BC. If != 0, set with value*TGV, where TGV=1e40
    for (int idx : boundaryIndices) {
        rhs(idx) = 0.0 * 1e40;
    }

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

    MeshData<3> mesh;
    // mesh.convertVTK("mesh_cubo.vtk", "mesh_cubo_3.vtk");
    mesh.readMesh("mesh_cubo_3.vtk");

    std::cout << "numNodes: \n" << mesh.getNumNodes() << "\n" << std::endl;
    std::cout << "nodes cols: \n" << mesh.getNodes().cols() << "\n" << std::endl;
    std::cout << "nodes rows: \n" << mesh.getNodes().rows() << "\n" << std::endl;
    // std::cout << "connectivity: \n" << mesh.getConnectivity() << "\n" << std::endl;
    std::cout << "numElements: \n" << mesh.getNumElements() << "\n" << std::endl;
    
    mesh.updateBoundaryNodes();
    const std::vector<long int> boundaryIndices = mesh.getBoundaryNodes();
    std::cout << "Boundary nodes count: \n" << boundaryIndices.size() << "\n" << std::endl;
    // for (int node : boundaryIndices) {
    //     std::cout << "Boundary node: " << node << std::endl;
    // }

    // prepare global matrices
    Eigen::SparseMatrix<double> stiffnessMatrix(mesh.getNumNodes(), mesh.getNumNodes());
    Eigen::SparseMatrix<double> massMatrix(mesh.getNumNodes(), mesh.getNumNodes());
    Eigen::SparseMatrix<double> gradientMatrix(mesh.getNumNodes(), 3);

    std::vector<Eigen::Triplet<double>> triplets;
    // std::cout << "max allocable size: \n" << triplets.max_size() << "\n" << std::endl;
    triplets.reserve(mesh.getNumElements());

    apsc::LinearFiniteElement<3> linearFiniteElement;
    linearFiniteElement.initializeRefGradients();

    Nodes localNodes;
    Indexes globalNodeNumbers;

    for (auto k = 0; k < mesh.getNumElements(); ++k) {
        // std::cout << k << std::endl;
        // extract element data
        for (auto i = 0; i < 4; ++i) // node numbering
        {
        globalNodeNumbers(i) = mesh.getConnectivity()(i, k);
        for (auto j = 0; j < 3; ++j) // local node coordinates
        {
            localNodes(j, i) = mesh.getNodes()(j, mesh.getConnectivity()(i, k)); // localNodes(j, i) = nodes(j, globalNodeNumbers(i));
        }
        }

        linearFiniteElement.update(localNodes);
        linearFiniteElement.updateGlobalNodeNumbers(globalNodeNumbers);
        
        linearFiniteElement.computeLocalStiffness();
        linearFiniteElement.updateGlobalStiffnessMatrix(stiffnessMatrix);
        
        linearFiniteElement.computeLocalMass();
        linearFiniteElement.updateGlobalMassMatrix(massMatrix);
        
        linearFiniteElement.computeLocalGradient();
        linearFiniteElement.updateGlobalGradientMatrix(gradientMatrix);
    }
    
    // Impose Dirichlet BC on stiffness matrix
    linearFiniteElement.updateMatrixWithDirichletBoundary(stiffnessMatrix, boundaryIndices);

    // Solve the Eikonal Equation
    // Solve the Heat Equation for initial conditions
    Values forcingTerm = Values::Constant(mesh.getNumNodes(), 1.0);
    // std::vector<int> boundaryIndices = {0};
    InitialConditions<3> initialConditions;
    Values initial_conditions = initialConditions.HeatEquation(stiffnessMatrix, boundaryIndices);
    
    // Values initial_conditions = Values::Constant(mesh.getNumNodes(), 5.0);
    std::cout << "initial_conditions:" << std::endl;
    std::cout << "Rows:" << initial_conditions.rows() << std::endl;
    std::cout << "Cols:" << initial_conditions.cols() << std::endl;
    // std::cout << initial_conditions << std::endl;

    Values w = initial_conditions;

    bool converged = false;
    int maxIterations = 1000;

    for (int iter = 0; iter < maxIterations && !converged; ++iter) {
        // std::cout << "-------------- Iteration " << iter + 1 << "--------------" << std::endl;
        converged = updateSolution(w, stiffnessMatrix, gradientMatrix, boundaryIndices);
        if (converged) {
            std::cout << "Solution converged after " << iter + 1 << " iterations." << std::endl;
            // std::cout << w << std::endl;
        }
    }

    if (!converged) {
        std::cout << "Solution did not converge within the maximum number of iterations." << std::endl;
    }

    return 0;
}

// // cd /home/jammy/shared-folder/project/forStudents
// // make
// // ./mainFem



