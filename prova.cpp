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
#include <Eigen/Sparse>
#include<Eigen/SparseCholesky>

using Values = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using Gradients = Eigen::Matrix<double, Eigen::Dynamic, 3>;
using Elements = Eigen::Matrix<long int, 3+1, Eigen::Dynamic>;


auto computeGradient(const Eigen::Matrix<double, 3, 3> &Coefficients,
                     const Eigen::Matrix<double, 3 + 1, 1> &values)
    -> Eigen::Matrix<double, 3, 1> {
  // Compute the differences in values
  Eigen::Matrix<double, 3, 1> d_values;
  for (auto i = 0u; i < 3; ++i) {
    d_values(i) = values(i + 1) - values(0);
  }
  // Compute the gradient
  return Coefficients * d_values;
}

//  ================== INCREMENTAL SOLUTION ==================

// Function to update the solution vector based on the Eikonal equation
bool updateSolution(Values& w, 
                    const Eigen::SparseMatrix<double>& stiffnessMatrix, 
                    const std::vector<Eigen::Matrix<double, 3, 3>>& gradientCoeff,
                    const Elements& connectivity,
                    const std::vector<long int> &boundaryIndices,
                    const auto& solver,
                    double gamma = 1e-3,
                    double tol = 1e-5) {
    
    std::cout << "Entered in updateSolution:" << std::endl;

    // VERSION NODEWISE ---------------------------------------------------------------------

    // Gradients grad_w(gradientMatrix.rows(), 3);
    
    // Gradients w_concatenated = w.replicate(1, 3);
    
    // grad_w = gradientMatrix.cwiseProduct(w_concatenated); // num_nodi x 3
    
    // Values norm_grad_w = grad_w.rowwise().norm(); // num_nodi x 1

    // Values gamma_vec = Values::Constant(gradientMatrix.rows(), gamma); // num_nodi x 1
    
    // Values coeffs = ((Values::Constant(gradientMatrix.rows(), 1.0) - norm_grad_w).array() / (norm_grad_w + gamma_vec).array()).matrix();

    // Values rhs = coeffs * (stiffnessMatrix * w);
    
    // ---------------------------------------------------------------------------------------


    // VERSION ELEMENTWISE ====================================================================
    
    // Local objects
    Values local_w(connectivity.rows()); //  = w[connectivity.col(i)]; // 4x1
    Eigen::SparseMatrix<double> localStiffness(connectivity.rows(), connectivity.rows()); // 4x4
    Values local_rhs(connectivity.rows()); // 4x1
    
    // Global objects
    Gradients grad_w(gradientCoeff.size(), 3); // num_elements x 3
    Values rhs = Values::Zero(stiffnessMatrix.rows()); // num_nodi x 1
    
    for(int i = 0; i < gradientCoeff.size(); i++) {
        
        for (int j = 0; j < connectivity.col(i).size(); j++) {
            local_w(j) = w(connectivity.col(i)(j));
        }
        
        grad_w.row(i) = computeGradient(gradientCoeff[i], local_w);

        double norm_grad_w = grad_w.row(i).norm(); // scalar
        // std::cout << "norm_grad_w: " << norm_grad_w << std::endl;

        double coeff = (1.0 - norm_grad_w) / (norm_grad_w + gamma);
        // std::cout << "coeff: " << coeff << std::endl;

        for (int k = 0; k < connectivity.col(i).size(); ++k) {
            for (int j = 0; j < connectivity.col(i).size(); ++j) {
                localStiffness.coeffRef(k, j) = stiffnessMatrix.coeff(connectivity.col(i)(k), connectivity.col(i)(j));
            }
        }

        local_rhs = coeff * (localStiffness * local_w); // matrice locale, elemento per elemento --> riasseblare a posteriori

        for (int j = 0; j < connectivity.col(i).size(); ++j) {
            rhs(connectivity.col(i)(j)) += local_rhs(j);
        }

    }

    // ========================================================================================

    // Update for Dirichlet BC. If == 0, set with value*TGV, where TGV=1e40
    for (int idx : boundaryIndices) {
        rhs(idx) = 0.0 * 1e40;
    }

    double max_rhs = rhs.maxCoeff();
    std::cout << "max_rhs: " << max_rhs << std::endl;

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
    std::cout << "z : " << z.norm() << std::endl;
    return (z.norm() < tol);

}

//  ================== MAIN ==================
int main() {
  
    using Nodes = apsc::LinearFiniteElement<3>::Nodes;
    using Indexes = apsc::LinearFiniteElement<3>::Indexes;
    using Values = Eigen::Matrix<double, Eigen::Dynamic, 1>; //=< The array storing values at the the nodes
    using AnisotropyM = Eigen::Matrix<double, 3, 3>; //=< The array storing the anisotropy matrix M 

    MeshData<3> mesh;
    const std::string mesh_file = "mesh_uniform_3.vtk";
    mesh.convertVTK("mesh_uniform.vtk", mesh_file);
    mesh.readMesh(mesh_file);

    std::cout << "numNodes: \n" << mesh.getNumNodes() << "\n" << std::endl;
    std::cout << "nodes cols: \n" << mesh.getNodes().cols() << "\n" << std::endl;
    std::cout << "nodes rows: \n" << mesh.getNodes().rows() << "\n" << std::endl;
    // std::cout << "connectivity: \n" << mesh.getConnectivity() << "\n" << std::endl;
    std::cout << "numElements: \n" << mesh.getNumElements() << "\n" << std::endl;
    
    mesh.updateBoundaryNodes();
    const std::vector<long int> boundaryIndices = mesh.getBoundaryNodes();
    std::cout << "Boundary nodes count: \n" << boundaryIndices.size() << "\n" << std::endl;

    // Prepare global matrices
    Eigen::SparseMatrix<double> stiffnessMatrix(mesh.getNumNodes(), mesh.getNumNodes());
    Eigen::SparseMatrix<double> massMatrix(mesh.getNumNodes(), mesh.getNumNodes());
    Eigen::SparseMatrix<double> gradientMatrixNodewise(mesh.getNumNodes(), 3);
    
    // std::vector<Eigen::Matrix<double, 3, 3+1>> gradientMatrixElementwise(mesh.getNumElements());
    std::vector<Eigen::Matrix<double, 3, 3>> gradientCoeff(mesh.getNumElements());

    Values globalIntegrals = Values::Constant(mesh.getNumNodes(), 0.0);

    std::vector<Eigen::Triplet<double>> triplets;
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

        // Compute local nodes and update global node numbers
        linearFiniteElement.update(localNodes);
        linearFiniteElement.updateGlobalNodeNumbers(globalNodeNumbers);

        linearFiniteElement.computeLocalIntegral();
        linearFiniteElement.updateGlobalIntegrals(globalIntegrals);
        
        // Compute the local stiffness matrix and update the global stiffness matrix
        linearFiniteElement.computeLocalStiffness();
        linearFiniteElement.updateGlobalStiffnessMatrix(stiffnessMatrix);
        
        // Compute the local mass matrix and update the global mass matrix
        linearFiniteElement.computeLocalMass();
        linearFiniteElement.updateGlobalMassMatrix(massMatrix);
        
        // Compute the local gradient matrix and update the global gradient matrix
        linearFiniteElement.computeLocalGradient();
        linearFiniteElement.updateGlobalGradientMatrix(gradientMatrixNodewise);

        gradientCoeff[k] = linearFiniteElement.computeGradientCoeff();

    }

    double max_stiffness = stiffnessMatrix.toDense().maxCoeff();
    std::cout << "max_stiffness: " << max_stiffness << std::endl;
    
    // Impose Dirichlet BC on stiffness matrix
    linearFiniteElement.updateMatrixWithDirichletBoundary(stiffnessMatrix, boundaryIndices);
    
    // -------- Solve the Eikonal Equation --------
    
    // Solve the Heat Equation for initial conditions
    Values forcingTerm = Values::Constant(mesh.getNumNodes(), 1.0);
    InitialConditions<3> initialConditions;
    Values initial_conditions = initialConditions.HeatEquation(stiffnessMatrix, 
                                                                globalIntegrals,
                                                                forcingTerm, 
                                                                boundaryIndices);
    
    mesh.addScalarField(initial_conditions, mesh_file);

    // Solve the Eikonal Equation Iteratively
    Values w = initial_conditions;

    bool converged = false;
    int maxIterations = 10000;

    // Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> solver;
    // Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver(stiffnessMatrix);

    solver.compute(stiffnessMatrix); // con metodo diretto, no decomposizione solo all'inizio
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

    for (int iter = 0; iter < maxIterations && !converged; ++iter) {
        std::cout << "-------------- Iteration " << iter + 1 << " --------------" << std::endl;
        converged = updateSolution(w, stiffnessMatrix, 
                                    gradientCoeff, 
                                    mesh.getConnectivity(), 
                                    boundaryIndices,
                                    solver);
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



