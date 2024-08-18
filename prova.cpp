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


//  ================== MAIN ==================
int main() {
    
    const std::size_t PHDIM = 3;

    using Traits = Eikonal::Eikonal_traits<PHDIM>;
    using Values = typename Traits::Values;
    using Nodes = typename Traits::Nodes;
    
    // Read the mesh file and create the mesh object
    MeshData<PHDIM> mesh;
    const std::string mesh_file = "vtk_files/mesh_uniform_3.vtk";
    // mesh.convertVTK("mesh_uniform.vtk", mesh_file);
    mesh.readMesh(mesh_file);
    mesh.updateBoundaryNodes();
    const std::vector<long int> boundaryIndices = mesh.getBoundaryNodes();
    
    // Prepare global matrices
    Eigen::SparseMatrix<double> stiffnessMatrix(mesh.getNumNodes(), mesh.getNumNodes());
    Eigen::SparseMatrix<double> massMatrix(mesh.getNumNodes(), mesh.getNumNodes());
    std::vector<std::vector<Eigen::Matrix<double, PHDIM, 1>>> reactionMatrix(mesh.getNumNodes(), std::vector<Eigen::Matrix<double, PHDIM, 1>>(mesh.getNumNodes()));
    
    std::vector<Eigen::Matrix<double, PHDIM, PHDIM>> gradientCoeff(mesh.getNumElements());
    Values globalIntegrals = Values::Constant(mesh.getNumNodes(), 0.0);

    // std::vector<Eigen::Triplet<double>> triplets;
    // triplets.reserve(mesh.getNumElements());

    mesh.fillGlobalVariables(stiffnessMatrix, massMatrix, reactionMatrix, globalIntegrals, gradientCoeff, boundaryIndices);

    // Solve the Heat Equation for initial conditions
    Values forcingTerm = Values::Constant(mesh.getNumNodes(), 1.0);
    InitialConditions<PHDIM> initialConditions;
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

    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(stiffnessMatrix); // Decompose the stiffness matrix with a direct solver
    
    StandardEikonal<PHDIM> eikonal(mesh, w, stiffnessMatrix, gradientCoeff, boundaryIndices, solver);
    
    // double r = 0.001;
    // PenaltyEikonal<PHDIM> eikonal(mesh, w, stiffnessMatrix, gradientCoeff, boundaryIndices, solver, r);

    Eigen::Matrix<double, Eigen::Dynamic, PHDIM> lagrangians = Eigen::Matrix<double, Eigen::Dynamic, PHDIM>::Constant(mesh.getNumElements(), PHDIM, 0.0);

    for (int iter = 0; iter < maxIterations && !converged; ++iter) {
        // std::cout << "-------------- Iteration " << iter + 1 << " --------------" << std::endl;
        
        converged = eikonal.updateSolution();
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



