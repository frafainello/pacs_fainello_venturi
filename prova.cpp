#include "gradient.hpp"
#include "Eikonal_traits.hpp"
#include "LinearFiniteElement.hpp"
#include "MeshData.hpp"
#include "InitialConditions.hpp"
#include "EikonalEquation.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include<Eigen/SparseCholesky>

using Values = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using Gradients = Eigen::Matrix<double, Eigen::Dynamic, 3>;
using Elements = Eigen::Matrix<long int, 3+1, Eigen::Dynamic>;
using Nodes = apsc::LinearFiniteElement<3>::Nodes;
using Indexes = apsc::LinearFiniteElement<3>::Indexes;


// auto computeGradient(const Eigen::Matrix<double, 3, 3> &Coefficients,
//                      const Eigen::Matrix<double, 3 + 1, 1> &values)
//     -> Eigen::Matrix<double, 3, 1> {
//   // Compute the differences in values
//   Eigen::Matrix<double, 3, 1> d_values;
//   for (auto i = 0u; i < 3; ++i) {
//     d_values(i) = values(i + 1) - values(0);
//   }
//   // Compute the gradient
//   return Coefficients * d_values;
// }

// ================== INCREMENTAL SOLUTION ==================

// Function to update the solution vector based on the Eikonal equation
bool updateSolution(Values& w, 
                    const Eigen::SparseMatrix<double>& stiffnessMatrix, 
                    const std::vector<Eigen::Matrix<double, 3, 3>>& gradientCoeff,
                    const MeshData<3>& mesh,
                    const std::vector<long int> &boundaryIndices,
                    const auto& solver,
                    double gamma = 1e-3,
                    double tol = 1e-4) {
    
    std::cout << "Entered in updateSolution:" << std::endl;

    // VERSION NODEWISE (neglect by now) -----------------------------------------------------
    // Gradients grad_w(gradientMatrix.rows(), 3);
    // Gradients w_concatenated = w.replicate(1, 3);
    // grad_w = gradientMatrix.cwiseProduct(w_concatenated); // num_nodi x 3
    // Values norm_grad_w = grad_w.rowwise().norm(); // num_nodi x 1
    // Values gamma_vec = Values::Constant(gradientMatrix.rows(), gamma); // num_nodi x 1
    // Values coeffs = ((Values::Constant(gradientMatrix.rows(), 1.0) - norm_grad_w).array() / (norm_grad_w + gamma_vec).array()).matrix();
    // Values rhs = coeffs * (stiffnessMatrix * w);
    // ---------------------------------------------------------------------------------------


    // VERSION ELEMENTWISE ===================================================================
    
    // Local objects
    Values local_w(mesh.getConnectivity().rows()); //  = w[connectivity.col(i)]; // 4x1
    Eigen::SparseMatrix<double> localStiffness(mesh.getConnectivity().rows(), mesh.getConnectivity().rows()); // 4x4
    Values local_rhs(mesh.getConnectivity().rows()); // 4x1
    
    // Global objects
    Gradients grad_w(gradientCoeff.size(), 3); // num_elements x 3
    Values norm_grad_w(gradientCoeff.size(), 1); // num_elements x 1
    Values rhs = Values::Zero(stiffnessMatrix.rows()); // num_nodi x 1

    apsc::LinearFiniteElement<3> linearFiniteElement;
    
    Nodes localNodes;
    Indexes globalNodeNumbers;
    
    for(int i = 0; i < gradientCoeff.size(); i++) {
        
        for (int j = 0; j < mesh.getConnectivity().col(i).size(); j++) {
            local_w(j) = w(mesh.getConnectivity().col(i)(j));
        }
        
        grad_w.row(i) = apsc::computeGradient<3>(gradientCoeff[i], local_w);
        // std::cout << "grad_w(i): " << grad_w.row(i) << std::endl;

        double norm_grad = grad_w.row(i).norm();
        norm_grad_w.coeffRef(i) = norm_grad; // scalar
        // std::cout << "norm_grad_w: " << norm_grad << std::endl;

        double coeff = (1.0 - norm_grad) / (norm_grad + gamma);
        // std::cout << "coeff: " << coeff << std::endl;

        // for (int k = 0; k < connectivity.col(i).size(); ++k) {
        //     for (int j = 0; j < connectivity.col(i).size(); ++j) {
        //         localStiffness.coeffRef(k, j) = stiffnessMatrix.coeff(connectivity.col(i)(k), connectivity.col(i)(j));
        //     }
        // }

        for (auto k = 0; k < 4; ++k) // node numbering
            {
            globalNodeNumbers(k) = mesh.getConnectivity()(k, i);
            for (auto j = 0; j < 3; ++j) // local node coordinates
            {
                localNodes(j, k) = mesh.getNodes()(j, mesh.getConnectivity()(k, i)); // localNodes(j, i) = nodes(j, globalNodeNumbers(i));
            }
        }
        
        // Compute local nodes and update global node numbers
        linearFiniteElement.update(localNodes);
        linearFiniteElement.updateGlobalNodeNumbers(globalNodeNumbers);
        
        // Compute the local stiffness matrix and update the global stiffness matrix
        linearFiniteElement.computeLocalStiffness();
        // linearFiniteElement.updateGlobalStiffnessMatrix(stiffnessMatrix);

        local_rhs = coeff * (linearFiniteElement.getLocalStiffness() * local_w); // contributo locale di ogni elemento al rhs globale
        // std::cout << "localStiffness: " << std::endl;
        // std::cout << localStiffness << std::endl;
        // std::cout << "local_rhs: " << std::endl;
        // std::cout << local_rhs << std::endl;
        // std::cout << "\n =========" << std::endl;

        for (int j = 0; j < mesh.getConnectivity().col(i).size(); ++j) {
            rhs(mesh.getConnectivity().col(i)(j)) += local_rhs(j);
        }
    }

    // =======================================================================================

    // Update for Dirichlet BC. If == 0, set with value*TGV, where TGV=1e40
    for (int idx : boundaryIndices) {
        rhs(idx) = 0.0 * 1e40;
    }

    // std::cout << "max_grad: " << norm_grad_w.maxCoeff() << std::endl;
    // std::cout << "min_grad: " << norm_grad_w.minCoeff() << std::endl;

    // std::cout << "max_rhs: " << rhs.maxCoeff() << std::endl;
    // std::cout << "min_rhs: " << rhs.minCoeff() << std::endl;

    // Solve the lienar system
    Values z = solver.solve(rhs);
    switch (solver.info()) {
        case Eigen::Success:
            std::cout << "Solving succeeded." << std::endl;
            break;
        case Eigen::NumericalIssue:
            std::cout << "Numerical issue." << std::endl;
            break;
        case Eigen::NoConvergence:
            std::cout << "No convergence." << std::endl;
            break;
        case Eigen::InvalidInput:
            std::cout << "Invalid input." << std::endl;
            break;
        default:
            std::cout << "Unknown error." << std::endl;
    }

    // Update the solution
    w += z;
    // std::cout << "z : " << z.norm() << std::endl;
    return (z.norm() < tol);

}

//  ================== MAIN ==================
int main() {
  
    using Nodes = apsc::LinearFiniteElement<3>::Nodes;
    using Indexes = apsc::LinearFiniteElement<3>::Indexes;
    using Values = Eigen::Matrix<double, Eigen::Dynamic, 1>; //=< The array storing values at the the nodes
    using AnisotropyM = Eigen::Matrix<double, 3, 3>; //=< The array storing the anisotropy matrix M 

    // Read the mesh file and create the mesh object
    MeshData<3> mesh;
    const std::string mesh_file = "vtk_files/mesh_uniform_3.vtk";
    // mesh.convertVTK("mesh_uniform.vtk", mesh_file);
    mesh.readMesh(mesh_file);
    mesh.updateBoundaryNodes();
    const std::vector<long int> boundaryIndices = mesh.getBoundaryNodes();
    
    // Prepare global matrices
    Eigen::SparseMatrix<double> stiffnessMatrix(mesh.getNumNodes(), mesh.getNumNodes());
    
    std::vector<Eigen::Matrix<double, 3, 3>> gradientCoeff(mesh.getNumElements());
    Values globalIntegrals = Values::Constant(mesh.getNumNodes(), 0.0);

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(mesh.getNumElements());

    apsc::LinearFiniteElement<3> linearFiniteElement;
    
    Nodes localNodes;
    Indexes globalNodeNumbers;

    for (auto k = 0; k < mesh.getNumElements(); ++k) {

        for (auto i = 0; i < 4; ++i) // node numbering
            {
            globalNodeNumbers(i) = mesh.getConnectivity()(i, k);
            for (auto j = 0; j < 3; ++j) // local node coordinates
            {
                localNodes(j, i) = mesh.getNodes()(j, mesh.getConnectivity()(i, k)); // localNodes(j, i) = nodes(j, globalNodeNumbers(i));
            }
        }
        
        // Compute local nodes and update global node numbers
        linearFiniteElement.update(localNodes);
        linearFiniteElement.updateGlobalNodeNumbers(globalNodeNumbers);

        linearFiniteElement.computeLocalIntegral();
        linearFiniteElement.updateGlobalIntegrals(globalIntegrals);
        
        // Compute the local stiffness matrix and update the global stiffness matrix
        linearFiniteElement.computeLocalStiffness();
        linearFiniteElement.updateGlobalStiffnessMatrix(stiffnessMatrix);
        
        // Compute the local gradient coefficients
        gradientCoeff[k] = linearFiniteElement.computeGradientCoeff();
        
    }

    // Impose Dirichlet BC on Stiffness Matrix
    linearFiniteElement.updateMatrixWithDirichletBoundary(stiffnessMatrix, boundaryIndices);
    
    // Solve the Heat Equation for initial conditions
    Values forcingTerm = Values::Constant(mesh.getNumNodes(), 1.0);
    InitialConditions<3> initialConditions;
    Values initial_conditions = initialConditions.HeatEquation(stiffnessMatrix, 
                                                                globalIntegrals,
                                                                forcingTerm, 
                                                                boundaryIndices);
    
    std::string suffix = "ic";
    mesh.addScalarField(initial_conditions, mesh_file, suffix);

    // Solve the Eikonal Equation Iteratively
    Values w = initial_conditions;

    bool converged = false;
    int maxIterations = 10000;

    // Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> solver;
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    // Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver(stiffnessMatrix);
    solver.compute(stiffnessMatrix); // Decompose the stiffness matrix with a direct solver

    for (int iter = 0; iter < maxIterations && !converged; ++iter) {
        std::cout << "-------------- Iteration " << iter + 1 << " --------------" << std::endl;
        converged = updateSolution(w, stiffnessMatrix, 
                                    gradientCoeff, 
                                    mesh,
                                    // mesh.getConnectivity(), 
                                    boundaryIndices,
                                    solver);
        if (converged) {
            std::cout << "Solution converged after " << iter + 1 << " iterations." << std::endl;
            mesh.addScalarField(w, mesh_file, "solution");
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



