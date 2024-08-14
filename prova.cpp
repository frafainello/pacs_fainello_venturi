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
    Indexes globalNodeNumbers(3+1);
    
    for (auto k = 0; k < mesh.getNumElements(); ++k) {

        for (auto i = 0; i < 4; ++i) // node numbering
            {
            globalNodeNumbers[i] = mesh.getConnectivity()(i, k);
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

    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(stiffnessMatrix); // Decompose the stiffness matrix with a direct solver
    
    StandardEikonal<3> eikonal(mesh, w, stiffnessMatrix, gradientCoeff, boundaryIndices, solver);

    for (int iter = 0; iter < maxIterations && !converged; ++iter) {
        std::cout << "-------------- Iteration " << iter + 1 << " --------------" << std::endl;
        
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



