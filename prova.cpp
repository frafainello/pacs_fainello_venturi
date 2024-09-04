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
    // const std::string mesh_file = "vtk_files/mesh_uniform_3.vtk";
    std::string mesh_file;
    std::cout << "Enter the mesh file name (default: mesh_uniform_3.vtk): ";
    std::getline(std::cin, mesh_file);

    // Append directory to filename if not empty, otherwise use default file
    if (!mesh_file.empty()) {
        mesh_file = "vtk_files/" + mesh_file;
    } else {
        mesh_file = "vtk_files/mesh_uniform_3.vtk";
    }
    // Check if the file can be opened
    std::ifstream file(mesh_file);
    if (!file) {
        std::cerr << "Error: File '" << mesh_file << "' not found. Please check the filename and try again.\n";
        return 1; // Exit with error code 1
    } else {
        std::cout << "Using mesh file: " << mesh_file << std::endl;
    }

    // mesh.convertVTK("mesh_uniform.vtk", mesh_file);
    mesh.readMesh(mesh_file);

    int choice = 1; // Initialize choice with the default value
    std::cout << "Enter the Dirichlet Boundary Conditions:\n";
    std::cout << "1. Dirichlet Null BCs on the bottom face of the domain (z=0)\n";
    std::cout << "2. Dirichlet Null BCs on the central point of the bottom face (z=0, x=0.5, y=0.5)\n";
    std::cout << "3. Dirichlet Null BCs on one vertex of the domain (x=0, y=0, z=0)\n";
    std::cout << "4. Dirichlet Null BCs on the whole boundary of the domain\n";
    std::cout << "Please choose an option (1-4, default 1): ";

    if (!(std::cin >> choice) || choice < 1 || choice > 4) {
        std::cin.clear(); // Clear the error flag that may be set by an invalid input
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Ignore the rest of the incorrect input
        choice = 1; // Set to default if input is invalid or not within the range
        std::cout << "Invalid input. Defaulting to option 1." << std::endl;
    }

    mesh.updateBoundaryNodes(choice);
    
    // Prepare global matrices
    Eigen::SparseMatrix<double> stiffnessMatrix(mesh.getNumNodes(), mesh.getNumNodes());
    Eigen::SparseMatrix<double> massMatrix(mesh.getNumNodes(), mesh.getNumNodes());
    std::vector<std::vector<Eigen::Matrix<double, PHDIM, 1>>> reactionMatrix(mesh.getNumNodes(), std::vector<Eigen::Matrix<double, PHDIM, 1>>(mesh.getNumNodes()));
    
    std::vector<Eigen::Matrix<double, PHDIM, PHDIM>> gradientCoeff(mesh.getNumElements());
    Values globalIntegrals = Values::Constant(mesh.getNumNodes(), 0.0);

    mesh.fillGlobalVariables(stiffnessMatrix, massMatrix, reactionMatrix, globalIntegrals, gradientCoeff);

    // Solve the Heat Equation for initial conditions
    Values forcingTerm = Values::Constant(mesh.getNumNodes(), 1.0);
    InitialConditions<PHDIM> initialConditions;
    Values initial_conditions = initialConditions.HeatEquation(stiffnessMatrix, 
                                                                globalIntegrals,
                                                                forcingTerm, 
                                                                mesh.getBoundaryNodes());
    
    mesh.addScalarField(initial_conditions, mesh_file, "heat", choice, false);

    // Solve the Eikonal Equation Iteratively
    Values w = initial_conditions;

    bool converged = false;
    int maxIterations = 10000;

    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(stiffnessMatrix); // Decompose the stiffness matrix with a direct solver
    
    // // StandardEikonal<PHDIM> eikonal(mesh, w, stiffnessMatrix, gradientCoeff, solver);
    
    // // double r = 0.001;
    // // PenaltyEikonal<PHDIM> eikonal(mesh, w, stiffnessMatrix, gradientCoeff, solver, r);

    // double r = 5;
    // Eigen::Matrix<double, Eigen::Dynamic, PHDIM> lagrangians = Eigen::Matrix<double, Eigen::Dynamic, PHDIM>::Constant(mesh.getNumElements(), PHDIM, 0.0);
    // LagrangianEikonal<PHDIM> eikonal(mesh, w, stiffnessMatrix, gradientCoeff, solver, r, lagrangians);

    int methodChoice;
    std::cout << "Select the eikonal method to use:\n";
    std::cout << "1. Standard Eikonal\n";
    std::cout << "2. Penalty Eikonal\n";
    std::cout << "3. Lagrangian Eikonal\n";
    std::cout << "Enter your choice (1-3, default: 1): ";
    std::cin >> methodChoice;

    if (std::cin.fail() || methodChoice < 1 || methodChoice > 3) {
        methodChoice = 1; // Default to Standard Eikonal if input is invalid
        std::cin.clear(); // Clear error flag
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Ignore wrong input
    }

    StandardEikonal<PHDIM> eikonal(mesh, w, stiffnessMatrix, gradientCoeff, solver);
    if (methodChoice == 2) {
        double r_penalty = 0.001;
        PenaltyEikonal<PHDIM> eikonal(mesh, w, stiffnessMatrix, gradientCoeff, solver, r_penalty);
    } else if (methodChoice == 3) {
        double r_lagrangian = 5;
        Eigen::Matrix<double, Eigen::Dynamic, PHDIM> lagrangians = Eigen::Matrix<double, Eigen::Dynamic, PHDIM>::Constant(mesh.getNumElements(), PHDIM, 0.0);
        LagrangianEikonal<PHDIM> eikonal(mesh, w, stiffnessMatrix, gradientCoeff, solver, r_lagrangian, lagrangians);
    }

    for (int iter = 0; iter < maxIterations && !converged; ++iter) {
        std::cout << "-------------- Iteration " << iter + 1 << " --------------" << std::endl;
        
        converged = eikonal.updateSolution();
        if (converged) {
            std::cout << "Solution converged after " << iter + 1 << " iterations." << std::endl;
            mesh.addScalarField(w, mesh_file, "heat", choice, true);
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



