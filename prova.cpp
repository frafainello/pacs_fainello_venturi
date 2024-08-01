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
                    // const Eigen::SparseMatrix<double>& gradientMatrix,
                    // const std::vector<Eigen::Matrix<double, 3, 3+1>>& gradientMatrixElementwise,
                    const std::vector<Eigen::Matrix<double, 3, 3>>& gradientCoeffElementwise,
                    const Elements& connectivity,
                    const std::vector<long int> &boundaryIndices,
                    const auto& solver,
                    double gamma = 1e-3,
                    double tol = 1e-6) {
    
    std::cout << "Entered in updateSolution:" << std::endl;

    // std::cout << "w: \n" << w << std::endl;
    // std::cout << "w.rows(): " << w.rows() << std::endl;
    // std::cout << "w.cols(): " << w.cols() << std::endl;
    // std::cout << "\n" << std::endl;

    // std::cout << "stiffnessMatrix: \n" << stiffnessMatrix << std::endl;
    // std::cout << "stiffnessMatrix.rows(): " << stiffnessMatrix.rows() << std::endl;
    // std::cout << "stiffnessMatrix.cols(): " << stiffnessMatrix.cols() << std::endl;
    // std::cout << "\n" << std::endl;

    // std::cout << "gradientMatrix: \n" << gradientMatrix << std::endl;
    // std::cout << "gradientMatrix.rows(): " << gradientMatrix.rows() << std::endl;
    // std::cout << "gradientMatrix.cols(): " << gradientMatrix.cols() << std::endl;
    // std::cout << "\n" << std::endl;

    // Values bilinear_form = stiffnessMatrix * w;
    // std::cout << "bilinear_form: \n" << bilinear_form << std::endl;
    // std::cout << "bilinear_form.rows(): " << bilinear_form.rows() << std::endl;
    // std::cout << "bilinear_form.cols(): " << bilinear_form.cols() << std::endl;
    // std::cout << "\n" << std::endl;

    // VERSION NODEWISE ---------------------------------------------------------------------

    // Gradients grad_w(gradientMatrix.rows(), 3);
    // // std::cout << "grad_w: \n" << grad_w << std::endl;
    // // std::cout << "grad_w.rows(): " << grad_w.rows() << std::endl;
    // // std::cout << "grad_w.cols(): " << grad_w.cols() << std::endl;
    // // std::cout << "\n" << std::endl;

    // Gradients w_concatenated = w.replicate(1, 3);
    // // std::cout << "w_concatenated: \n" << w_concatenated << std::endl;
    // // std::cout << "w_concatenated.rows(): " << w_concatenated.rows() << std::endl;
    // // std::cout << "w_concatenated.cols(): " << w_concatenated.cols() << std::endl;
    // // std::cout << "\n" << std::endl;

    // grad_w = gradientMatrix.cwiseProduct(w_concatenated); // num_nodi x 3
    // // std::cout << "grad_w: \n" << grad_w << std::endl;
    // // std::cout << "grad_w.rows(): " << grad_w.rows() << std::endl;
    // // std::cout << "grad_w.cols(): " << grad_w.cols() << std::endl;
    // // std::cout << "\n" << std::endl;
    // // Eigen::MatrixXd dense_grad_w = Eigen::MatrixXd(grad_w);

    // Values norm_grad_w = grad_w.rowwise().norm(); // num_nodi x 1
    // // std::cout << "norm_grad_w: " << norm_grad_w << std::endl;
    // // std::cout << "norm_grad_w.rows(): " << norm_grad_w.rows() << std::endl;
    // // std::cout << "norm_grad_w.cols(): " << norm_grad_w.cols() << std::endl;
    // // std::cout << "\n" << std::endl;

    // Values gamma_vec = Values::Constant(gradientMatrix.rows(), gamma); // num_nodi x 1
    // // std::cout << "gamma_vec: " << gamma_vec << std::endl;
    // // std::cout << "gamma_vec.rows(): " << gamma_vec.rows() << std::endl;
    // // std::cout << "gamma_vec.cols(): " << gamma_vec.cols() << std::endl;
    // // std::cout << "\n" << std::endl;
    
    // // Values coeffs = ((Values::Constant(gradientMatrix.cols(), 1.0) - norm_grad_w).array() / (norm_grad_w + gamma_vec).array()).matrix();
    // Values coeffs = ((Values::Constant(gradientMatrix.rows(), 1.0) - norm_grad_w).array() / (norm_grad_w + gamma_vec).array()).matrix();
    // // std::cout << "coeffs: " << coeffs << std::endl;
    // // std::cout << "coeffs.rows(): " << coeffs.rows() << std::endl;
    // // std::cout << "coeffs.cols(): " << coeffs.cols() << std::endl;
    // // std::cout << "\n" << std::endl;

    // // Compute the right-hand side for the linear problem
    // Values rhs = coeffs * (stiffnessMatrix * w);
    // ---------------------------------------------------------------------------------------


    // VERSION ELEMENTWISE ====================================================================
    // Gradients grad_w_old(gradientMatrixElementwise.size(), 3); // num_elements x 3
    // for(int i = 0; i < gradientMatrixElementwise.size(); i++) {
    //     // std::cout << "i: " << i << std::endl;
    //     // std::cout << "connectivity.col(i): " << connectivity.col(i) << std::endl;
    //     // grad_w.row(i) = gradientMatrixElementwise[i] * w[connectivity.col(i)]; // 3x4 * 4x1 = 3x1
        
    //     Values local_w(connectivity.col(i).size()); //  = w[connectivity.col(i)]; // 4x1
    //     for (int j = 0; j < connectivity.col(i).size(); j++) {
    //         local_w(j) = w(connectivity.col(i)(j));
    //     }
    //     grad_w_old.row(i) = gradientMatrixElementwise[i] * local_w; // 3x4 * 4x1 = 3x1
    // }
    // std::cout << "grad_w_old:" << std::endl;
    // std::cout << "grad_w_old.rows(): " << grad_w_old.rows() << std::endl;
    // std::cout << "grad_w_old.cols(): " << grad_w_old.cols() << std::endl;
    // std::cout << grad_w_old << std::endl;
    // std::cout << "\n" << std::endl;

    // Local objects
    Values local_w(connectivity.rows()); //  = w[connectivity.col(i)]; // 4x1
    Eigen::SparseMatrix<double> localStiffness(connectivity.rows(), connectivity.rows()); // 4x4
    Values local_rhs(connectivity.rows()); // 4x1
    
    // Global objects
    Gradients grad_w(gradientCoeffElementwise.size(), 3); // num_elements x 3
    Values rhs = Values::Zero(stiffnessMatrix.rows()); // num_nodi x 1
    
    for(int i = 0; i < gradientCoeffElementwise.size(); i++) {
        // std::cout << "i: " << i << std::endl;
        // std::cout << "connectivity.col(i): " << connectivity.col(i) << std::endl;
        
        for (int j = 0; j < connectivity.col(i).size(); j++) {
            local_w(j) = w(connectivity.col(i)(j));
        }
        // std::cout << "local_w:" << std::endl;
        // std::cout << "local_w.rows(): " << local_w.rows() << std::endl;
        // std::cout << "local_w.cols(): " << local_w.cols() << std::endl;
        // std::cout << local_w << std::endl;
        // std::cout << "\n" << std::endl;

        grad_w.row(i) = computeGradient(gradientCoeffElementwise[i], local_w);
        
        std::cout << "\n" << std::endl;
        std::cout << "grad_w.row(i):" << std::endl;
        // std::cout << "grad_w.rows(): " << grad_w.rows() << std::endl;
        // std::cout << "grad_w.cols(): " << grad_w.cols() << std::endl;
        std::cout << grad_w.row(i) << std::endl;

        // if (grad_w.isApprox(grad_w_old)) {
        //     std::cout << "The matrices are approximately equal." << std::endl;
        // } else {
        //     std::cout << "The matrices are not equal." << std::endl;
        // }

        // double norm_grad_w_old = grad_w_old.norm(); // scalar
        double norm_grad_w = grad_w.row(i).norm(); // scalar
        // std::cout << "norm_grad_w_old: " << norm_grad_w_old << std::endl;
        std::cout << "norm_grad_w: " << norm_grad_w << std::endl;

        double coeff = (1.0 - norm_grad_w) / (norm_grad_w + gamma);
        std::cout << "coeff: " << coeff << std::endl;

        for (int k = 0; k < connectivity.col(i).size(); ++k) {
            for (int j = 0; j < connectivity.col(i).size(); ++j) {
                localStiffness.coeffRef(k, j) = stiffnessMatrix.coeff(connectivity.col(i)(k), connectivity.col(i)(j));
            }
        }
        // std::cout << "localStiffness:" << std::endl;
        // std::cout << "localStiffness.rows(): " << localStiffness.rows() << std::endl;
        // std::cout << "localStiffness.cols(): " << localStiffness.cols() << std::endl;
        // std::cout << localStiffness << std::endl;
        // std::cout << "\n" << std::endl;

        local_rhs = coeff * (localStiffness * local_w); // matrice locale, elemento per elemento --> riasseblare a posteriori
        // std::cout << "local_rhs:" << std::endl;
        // std::cout << "local_rhs.rows(): " << local_rhs.rows() << std::endl;
        // std::cout << "local_rhs.cols(): " << local_rhs.cols() << std::endl;
        // std::cout << local_rhs << std::endl;
        // std::cout << "\n" << std::endl;

        for (int j = 0; j < connectivity.col(i).size(); ++j) {
            rhs(connectivity.col(i)(j)) += local_rhs(j);
        }

    }



    // Values rhs = coeff * (stiffnessMatrix * w); // matrice locale, elemento per elemento --> riasseblare a posteriori
    // std::cout << "rhs: " << rhs << std::endl;
    // std::cout << "rhs.rows(): " << rhs.rows() << std::endl;
    // std::cout << "rhs.cols(): " << rhs.cols() << std::endl;
    // std::cout << "\n" << std::endl;
    // ========================================================================================

    // Update for Dirichlet BC. If == 0, set with value*TGV, where TGV=1e40
    for (int idx : boundaryIndices) {
        rhs(idx) = 0.0 * 1e40;
    }

    // std::cout << "rhs: " << rhs << std::endl;
    // std::cout << "rhs.rows(): " << rhs.rows() << std::endl;
    // std::cout << "rhs.cols(): " << rhs.cols() << std::endl;
    // std::cout << "\n" << std::endl;

    // Solve the lienar system
    Values z = solver.solve(rhs);
    if(solver.info() == Eigen::Success) {
        throw std::runtime_error("Solving failed=");
    }

    // Update the solution
    w += z;
    
    return (z.norm() < tol);

}

//  ================== MAIN ==================
int main() {
  
    using Nodes = apsc::LinearFiniteElement<3>::Nodes;
    using Indexes = apsc::LinearFiniteElement<3>::Indexes;
    using Values = Eigen::Matrix<double, Eigen::Dynamic, 1>; //=< The array storing values at the the nodes
    using AnisotropyM = Eigen::Matrix<double, 3, 3>; //=< The array storing the anisotropy matrix M 

    MeshData<3> mesh;
    const std::string mesh_file = "mesh_simple_3.vtk";
    mesh.convertVTK("mesh_simple.vtk", mesh_file);
    mesh.readMesh(mesh_file);

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
    Eigen::SparseMatrix<double> gradientMatrixNodewise(mesh.getNumNodes(), 3);
    
    // std::vector<Eigen::SparseMatrix<double>> gradientMatrixElementwise(3, 3+1);
    // gradientMatrixElementwise.reserve(mesh.getNumElements());
    std::vector<Eigen::Matrix<double, 3, 3+1>> gradientMatrixElementwise(mesh.getNumElements());
    std::vector<Eigen::Matrix<double, 3, 3>> gradientCoeffElementwise(mesh.getNumElements());

    Values globalIntegrals = Values::Constant(mesh.getNumNodes(), 0.0);

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
        // std::cout << "local gradient matrix: \n" << linearFiniteElement.getLocalGradient() << std::endl;
        gradientMatrixElementwise[k] = linearFiniteElement.getLocalGradient();
        linearFiniteElement.updateGlobalGradientMatrix(gradientMatrixNodewise);

        gradientCoeffElementwise[k] = linearFiniteElement.computeGradientCoeff();
        // std::cout << "gradientCoeffElementwise[" << k << "]: \n" << gradientCoeffElementwise[k] << std::endl;

    }

    std::cout << "globalIntegrals: \n" << globalIntegrals << std::endl;

    
    double max_value = -std::numeric_limits<double>::infinity();
    for (const auto& matrix : gradientMatrixElementwise) {
        double max_in_matrix = matrix.maxCoeff();
        if (max_in_matrix > max_value) {
            max_value = max_in_matrix;
        }
    }
    std::cout << "Maximum value: " << max_value << std::endl;

    // Print stiffness matrix
    std::cout << "Stiffness Matrix:" << std::endl;
    std::cout << "Rows:" << stiffnessMatrix.rows() << std::endl;
    std::cout << "Cols:" << stiffnessMatrix.cols() << std::endl;
    std::cout << "\n" << std::endl;
    // std::cout << stiffnessMatrix << std::endl;

    // Print mass matrix
    std::cout << "Mass Matrix:" << std::endl;
    std::cout << "Rows:" << massMatrix.rows() << std::endl;
    std::cout << "Cols:" << massMatrix.cols() << std::endl;
    std::cout << "\n" << std::endl;
    // std::cout << massMatrix << std::endl;

    // Print gradient matrix
    std::cout << "Gradient Matrix:" << std::endl;
    // std::cout << "Rows:" << gradientMatrixNodewise.rows() << std::endl;
    // std::cout << "Cols:" << gradientMatrixNodewise.cols() << std::endl;
    std::cout << "\n" << std::endl;
    // std::cout << gradientMatrix << std::endl;
    // gradientMatriix = gradientMatrix.transpose();

    // Print gradient matrix elementwise
    // std::cout << "Gradient Matrix Elementwise:" << std::endl;
    // for (size_t k = 0; k < mesh.getNumElements(); ++k) {
    //     std::cout << "Element " << k << " matrix:\n" << gradientMatrixElementwise[k] << "\n\n";
    // }
    // std::cout << "\n" << std::endl;

    // Impose Dirichlet BC on stiffness matrix
    linearFiniteElement.updateMatrixWithDirichletBoundary(stiffnessMatrix, boundaryIndices);
    // std::cout << "Stiffness Matrix after BC:" << std::endl;
    // std::cout << stiffnessMatrix << std::endl;

    // Solve the Eikonal Equation
    // Solve the Heat Equation for initial conditions
    Values forcingTerm = Values::Constant(mesh.getNumNodes(), 1.0);
    // std::vector<int> boundaryIndices = {0};
    InitialConditions<3> initialConditions;
    Values initial_conditions = initialConditions.HeatEquation(stiffnessMatrix, 
                                                                globalIntegrals,
                                                                forcingTerm, 
                                                                boundaryIndices);
    
    mesh.addScalarField(initial_conditions, mesh_file);

    // Values initial_conditions = Values::Constant(mesh.getNumNodes(), 5.0);
    std::cout << "initial_conditions:" << std::endl;
    std::cout << "Rows:" << initial_conditions.rows() << std::endl;
    std::cout << "Cols:" << initial_conditions.cols() << std::endl;
    // std::cout << initial_conditions << std::endl;

    // Values w = initial_conditions;

    // bool converged = false;
    // int maxIterations = 1000;

    // // Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> solver;
    // Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    // solver.compute(stiffnessMatrix); // con metodo diretto, no decomposizione solo all'inizio

    // // if(solver.info() == Eigen::Success) {
    // //     throw std::runtime_error("Decomposition failed=");
    // // }

    // for (int iter = 0; iter < maxIterations && !converged; ++iter) {
    //     std::cout << "-------------- Iteration " << iter + 1 << " --------------" << std::endl;
    //     converged = updateSolution(w, stiffnessMatrix, 
    //                             // gradientMatrixNodewise, gradientMatrixElementwise, 
    //                             gradientCoeffElementwise, 
    //                             mesh.getConnectivity(), 
    //                             boundaryIndices,
    //                             solver);
    //     if (converged) {
    //         std::cout << "Solution converged after " << iter + 1 << " iterations." << std::endl;
    //         // std::cout << w << std::endl;
    //     }
    // }

    // if (!converged) {
    //     std::cout << "Solution did not converge within the maximum number of iterations." << std::endl;
    // }

    return 0;
}

// // cd /home/jammy/shared-folder/project/forStudents
// // make
// // ./mainFem



