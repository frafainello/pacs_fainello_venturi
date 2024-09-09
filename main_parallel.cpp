#include "src/gradient.hpp"
#include "src/Eikonal_traits.hpp"
#include "src/LinearFiniteElement.hpp"
#include "src/MeshData.hpp"
#include "src/InitialConditions.hpp"
#include "src/EikonalEquation.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <mpi.h>

//  ================== MAIN ==================
/**
 * @brief Entry point of the program.
 * 
 * Initializes MPI and performs the following steps:
 * 1. Reads input mesh file and boundary conditions.
 * 2. Sets up mesh and computes global matrices.
 * 3. Solves initial conditions using the heat equation.
 * 4. Broadcasts necessary data to all MPI processes.
 * 5. Initializes and executes the chosen Eikonal solver.
 * 
 * @param argc Number of command-line arguments.
 * @param argv Array of command-line arguments.
 * @return int Exit status of the program.
 */
int main(int argc, char **argv) {

    // Initialize MPI
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const std::size_t PHDIM = 3;

    using Traits = Eikonal::Eikonal_traits<PHDIM>;
    using Values = typename Traits::Values;
    using Nodes = typename Traits::Nodes;
    
    MeshData<PHDIM> mesh;
    std::string mesh_file;

    //  ================== GET INPUT INSTRUCTIONS ==================
    // Only rank 0 will handle input and output operations
    if (rank == 0) {
        std::cout << "\nEnter the mesh file name (default: mesh_10.vtk): ";
        std::getline(std::cin, mesh_file);

        if (!mesh_file.empty()) {
            mesh_file = "vtk_files/" + mesh_file;
        } else {
            mesh_file = "vtk_files/mesh_10.vtk";
        }
       
        std::ifstream file(mesh_file);
        if (!file) {
            std::cerr << "Error: File '" << mesh_file << "' not found. Please check the filename and try again.\n";
            MPI_Abort(MPI_COMM_WORLD, 1); // Abort if file is not found
        } else {
            std::cout << "Using mesh file: " << mesh_file << std::endl;
        }

        // Broadcasting the mesh file name to all processes
        int filename_length = mesh_file.size();
        MPI_Bcast(&filename_length, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(mesh_file.data(), filename_length, MPI_CHAR, 0, MPI_COMM_WORLD);
    } else {
        // Other ranks receive the mesh file name
        int filename_length;
        MPI_Bcast(&filename_length, 1, MPI_INT, 0, MPI_COMM_WORLD);
        mesh_file.resize(filename_length);
        MPI_Bcast(mesh_file.data(), filename_length, MPI_CHAR, 0, MPI_COMM_WORLD);
    }

    int choice = 1; // Default value
    std::string input;

    if (rank == 0) {
        std::cout << "\nEnter the Dirichlet Boundary Conditions:\n";
        std::cout << "1. Dirichlet Null BCs on the bottom face of the domain (z=0)\n";
        std::cout << "2. Dirichlet Null BCs on the central point of the bottom face (z=0, x=0.5, y=0.5)\n";
        std::cout << "3. Dirichlet Null BCs on one vertex of the domain (x=0, y=0, z=0)\n";
        std::cout << "4. Dirichlet Null BCs on the whole boundary of the domain\n";
        std::cout << "Please choose an option (1-4, default 1): ";
        std::getline(std::cin, input);

        std::stringstream ss(input);
        if (!(ss >> choice) || choice < 1 || choice > 4) {
            choice = 1;
            std::cout << "Defaulting to option 1." << std::endl;
        }

        std::cout << "Selected option: " << choice << std::endl;
    }

    // Broadcast the choice from rank 0 to all other ranks
    MPI_Bcast(&choice, 1, MPI_INT, 0, MPI_COMM_WORLD);

    std::string choiceString = "";
    if (choice == 1) {
        choiceString = "face";
    } else if (choice == 2) {
        choiceString = "center";
    } else if (choice == 3) {
        choiceString = "vertex";
    } else if (choice == 4) {
        choiceString = "all";
    }
    
    int methodChoice = 1; // Default value

    if (rank == 0) {
        std::cout << "\nSelect the eikonal method to use:\n";
        std::cout << "1. Standard Eikonal\n";
        std::cout << "2. Penalty Eikonal\n";
        std::cout << "3. Lagrangian Eikonal\n";
        std::cout << "Enter your choice (1-3, default: 1): ";
        std::getline(std::cin, input);

        std::stringstream sss(input);
        if (!(sss >> methodChoice) || methodChoice < 1 || methodChoice > 3) {
            methodChoice = 1;
            std::cout << "Defaulting to option 1." << std::endl;
        }
        std::cout << "Selected option: " << methodChoice << std::endl;
    }

    // Broadcast the method choice to all processes
    MPI_Bcast(&methodChoice, 1, MPI_INT, 0, MPI_COMM_WORLD);

    std::string itMethString = "";
    if (methodChoice == 1) {
        itMethString = "_standard";
    } else if (methodChoice == 2) {
        itMethString = "_penalty";
    } else if (methodChoice == 3) {
        itMethString = "_lagrangian";
    }

    //  ================== FILL MESH INFORMATION ==================
    if (rank == 0) std::cout << "\nReading mesh...\n";
    mesh.readMesh(mesh_file);
    if (rank == 0) std::cout << "Computing boundary nodes...\n";
    mesh.updateBoundaryNodes(choice);
    
    //  ================== COMPUTE GLOBAL MATRICES AND INITIAL CONDITION ==================
    // Computing global matrices and initial condition in rank 0
    Eigen::SparseMatrix<double> stiffnessMatrix(mesh.getNumNodes(), mesh.getNumNodes());
    // Eigen::SparseMatrix<double> massMatrix(mesh.getNumNodes(), mesh.getNumNodes());
    std::vector<std::vector<Eigen::Matrix<double, PHDIM, 1>>> reactionMatrix(mesh.getNumNodes(), std::vector<Eigen::Matrix<double, PHDIM, 1>>(mesh.getNumNodes()));

    Values globalIntegrals = Values::Constant(mesh.getNumNodes(), 0.0);
    Values w;
    int w_size = 0;

    const std::string bc = choiceString;
    const std::string ic="heat";
    std::string outputFilePath = mesh_file.substr(0, mesh_file.find_last_of('.')) + "/initial_conditions/" + bc + "_" + ic + "_ic.vtk";


    if (rank == 0){
        std::cout << "Filling global matrices...\n";
        mesh.fillGlobalVariables(stiffnessMatrix, reactionMatrix, globalIntegrals);

        // Solve the Heat Equation for initial conditions
        std::cout << "Computing initial condition...\n";
        Values forcingTerm = Values::Constant(mesh.getNumNodes(), 1.0);
        InitialConditions<PHDIM> initialConditions;
        Values initial_conditions = initialConditions.HeatEquation(stiffnessMatrix, 
                                                                    globalIntegrals,
                                                                    forcingTerm, 
                                                                    mesh.getBoundaryNodes());
        
        // Random
        // std::srand(42);
        // Values initial_conditions = Values(mesh.getNumNodes());
        // initial_conditions.setRandom();
        
        // "Like" BC
        // Values initial_conditions = Values::Constant(mesh.getNumNodes(), 1.0);
        // // Set specified indices to 0.0
        // for (long int index : mesh.getBoundaryNodes()) {
        //     if (index >= 0 && index < mesh.getNumNodes()) {  // Ensure index is within bounds
        //         initial_conditions(static_cast<int>(index)) = 0.0;
        //     }
        // }


        mesh.addScalarField(initial_conditions, mesh_file, outputFilePath);

        w = initial_conditions;
        w_size = w.size();

    }

    //  ================== BROADCAST INFORMATION ==================
    if (rank == 0) std::cout << "Broadcasting to all ranks...\n";
    // Broadcast the size of w
    MPI_Bcast(&w_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Resize w on all other ranks
    if (rank != 0) {
        w.resize(w_size);
    }
    // Broadcast the content of w
    MPI_Bcast(w.data(), w_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Broadcasting stiffness matrix
    // Rank 0: flatten and send
    if (rank == 0) {
        std::vector<Eigen::Triplet<double>> tripletList;
        for (int k = 0; k < stiffnessMatrix.outerSize(); ++k) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(stiffnessMatrix, k); it; ++it) {
                tripletList.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value()));
            }
        }
        int tripletSize = tripletList.size();
        MPI_Bcast(&tripletSize, 1, MPI_INT, 0, MPI_COMM_WORLD);

        std::vector<int> rows(tripletSize), cols(tripletSize);
        std::vector<double> values(tripletSize);
        for (int i = 0; i < tripletSize; ++i) {
            rows[i] = tripletList[i].row();
            cols[i] = tripletList[i].col();
            values[i] = tripletList[i].value();
        }
        MPI_Bcast(rows.data(), tripletSize, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(cols.data(), tripletSize, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(values.data(), tripletSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    } else {
        int tripletSize;
        MPI_Bcast(&tripletSize, 1, MPI_INT, 0, MPI_COMM_WORLD);

        std::vector<int> rows(tripletSize), cols(tripletSize);
        std::vector<double> values(tripletSize);
        MPI_Bcast(rows.data(), tripletSize, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(cols.data(), tripletSize, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(values.data(), tripletSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        // Reconstruct the sparse matrix on other ranks
        Eigen::SparseMatrix<double> stiffnessMatrix(mesh.getNumNodes(), mesh.getNumNodes());
        std::vector<Eigen::Triplet<double>> tripletList(tripletSize);
        for (int i = 0; i < tripletSize; ++i) {
            tripletList[i] = Eigen::Triplet<double>(rows[i], cols[i], values[i]);
        }
        stiffnessMatrix.setFromTriplets(tripletList.begin(), tripletList.end());
    }

    // Broadcasting reaction matrix
    // Flatten reactionMatrix
    if (rank == 0) {
        std::vector<double> flatReactionMatrix;
        for (const auto& outerVec : reactionMatrix) {
            for (const auto& innerVec : outerVec) {
                for (int i = 0; i < PHDIM; ++i) {
                    flatReactionMatrix.push_back(innerVec(i));
                }
            }
        }
        int flatSize = flatReactionMatrix.size();
        MPI_Bcast(&flatSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(flatReactionMatrix.data(), flatSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    } else {
        int flatSize;
        MPI_Bcast(&flatSize, 1, MPI_INT, 0, MPI_COMM_WORLD);

        std::vector<double> flatReactionMatrix(flatSize);
        MPI_Bcast(flatReactionMatrix.data(), flatSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        // Reconstruct reactionMatrix on other ranks
        int idx = 0;
        for (int i = 0; i < mesh.getNumNodes(); ++i) {
            for (int j = 0; j < mesh.getNumNodes(); ++j) {
                for (int k = 0; k < PHDIM; ++k) {
                    reactionMatrix[i][j](k) = flatReactionMatrix[idx++];
                }
            }
        }
    }

    // Broadcasting mesh members
    // Flatten gradientCoeff
    if (rank == 0) {
        std::vector<double> flatGradientCoeff;
        for (const auto& mat : mesh.gradientCoeff) {
            for (int i = 0; i < mat.size(); ++i) {
                flatGradientCoeff.push_back(mat(i));
            }
        }
        int flatSize = flatGradientCoeff.size();
        MPI_Bcast(&flatSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(flatGradientCoeff.data(), flatSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        
    } else {
        int flatSize;
        MPI_Bcast(&flatSize, 1, MPI_INT, 0, MPI_COMM_WORLD);

        std::vector<double> flatGradientCoeff(flatSize);
        MPI_Bcast(flatGradientCoeff.data(), flatSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        // Reconstruct gradientCoeff on other ranks
        int idx = 0;

        for (int i = 0; i < mesh.getNumElements(); ++i) {
            Eigen::Matrix<double, PHDIM, PHDIM> tempMatrix;  // Temporary matrix
            for (int j = 0; j < PHDIM * PHDIM; ++j) {
                tempMatrix(j % PHDIM, j / PHDIM) = flatGradientCoeff[idx++];
            }
            mesh.gradientCoeff.push_back(tempMatrix);  // Assign the filled matrix to the vector
        }
    }

    // Flatten localStiffnessMatrices (std::vector<Eigen::Matrix<double, PHDIM+1, PHDIM+1>>)
    if (rank == 0) {
        std::vector<double> flatLocalStiffnessMatrices;
        for (const auto& matrix : mesh.localStiffnessMatrices) {
            for (int i = 0; i < (PHDIM + 1) * (PHDIM + 1); ++i) {
                flatLocalStiffnessMatrices.push_back(matrix(i));
            }
        }
        int flatSize = flatLocalStiffnessMatrices.size();
        MPI_Bcast(&flatSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(flatLocalStiffnessMatrices.data(), flatSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    } else {
        int flatSize;
        MPI_Bcast(&flatSize, 1, MPI_INT, 0, MPI_COMM_WORLD);

        std::vector<double> flatLocalStiffnessMatrices(flatSize);
        MPI_Bcast(flatLocalStiffnessMatrices.data(), flatSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        // Reconstruct localStiffnessMatrices on other ranks
        int idx = 0;
        for (int i = 0; i < mesh.getNumElements(); ++i) {
            Eigen::Matrix<double, PHDIM+1, PHDIM+1> tempMatrix;  // Temporary matrix
            for (int j = 0; j < (PHDIM + 1) * (PHDIM + 1); ++j) {
                tempMatrix(j % (PHDIM + 1), j / (PHDIM + 1)) = flatLocalStiffnessMatrices[idx++];
            }
            mesh.localStiffnessMatrices.push_back(tempMatrix);  // Assign the filled matrix to the vector
        }
    }


    // Flatten localReactionMatrices (std::vector<std::vector<std::vector<Eigen::Matrix<double, PHDIM, 1>>>>)
    if (rank == 0) {
        std::vector<double> flatLocalReactionMatrices;
        for (const auto& outerVec : mesh.localReactionMatrices) {
            for (const auto& middleVec : outerVec) {
                for (const auto& matrix : middleVec) {
                    for (int i = 0; i < PHDIM; ++i) {
                        flatLocalReactionMatrices.push_back(matrix(i));
                    }
                }
            }
        }
        int flatSize = flatLocalReactionMatrices.size();
        MPI_Bcast(&flatSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(flatLocalReactionMatrices.data(), flatSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    } else {
        int flatSize;
        MPI_Bcast(&flatSize, 1, MPI_INT, 0, MPI_COMM_WORLD);

        std::vector<double> flatLocalReactionMatrices(flatSize);
        MPI_Bcast(flatLocalReactionMatrices.data(), flatSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        // Reconstruct localReactionMatrices on other ranks
        int idx = 0;
        for (int i = 0; i < mesh.getNumElements(); ++i) {
            std::vector<std::vector<Eigen::Matrix<double, PHDIM, 1>>> tempTensor; 

            for (int j = 0; j < (PHDIM + 1); ++j) {
                std::vector<Eigen::Matrix<double, PHDIM, 1>> tempIntermediate; 

                for (int k = 0; k < (PHDIM + 1); ++k) {
                    Eigen::Matrix<double, PHDIM, 1> tempVec;  // Temporary matrix
                    for (int l = 0; l < PHDIM; ++l) {
                        tempVec(l) = flatLocalReactionMatrices[idx++];
                    }
                    tempIntermediate.push_back(tempVec);
                }
                tempTensor.push_back(tempIntermediate);
            }
            mesh.localReactionMatrices.push_back(tempTensor);
        }
    }

    //  ================== SOLVER AND EIKONAL CLASS ==================
    if (rank == 0) std::cout << "Initializing solver...\n";
    int maxIterations = 10000;
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    if (rank == 0){
        solver.compute(stiffnessMatrix); // Decompose the stiffness matrix with a direct solver
    }

    std::unique_ptr<EikonalEquation<PHDIM>> eikonal = nullptr;
    if (methodChoice == 1) {
        eikonal = std::make_unique<StandardEikonal<PHDIM>>(mesh, w, stiffnessMatrix, solver);
        if (rank == 0) std::cout << "\nStandard Eikonal selected." << std::endl;
    } else if (methodChoice == 2) {
        double r_penalty = 0.1;
        eikonal = std::make_unique<PenaltyEikonal<PHDIM>>(mesh, w, stiffnessMatrix, solver, r_penalty);
        if (rank == 0) std::cout << "\nPenalty Eikonal selected." << std::endl;
    } else if (methodChoice == 3) {
        double r_lagrangian = 5;
        Eigen::Matrix<double, Eigen::Dynamic, PHDIM> lagrangians = Eigen::Matrix<double, Eigen::Dynamic, PHDIM>::Constant(mesh.getNumElements(), PHDIM, 0.0);
        eikonal = std::make_unique<LagrangianEikonal<PHDIM>>(mesh, w, stiffnessMatrix, solver, r_lagrangian, lagrangians);
        if (rank == 0) std::cout << "\nLagrangian Eikonal selected." << std::endl;
    } else {
        if (rank == 0) std::cerr << "\nInvalid choice. Exiting..." << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    
    //  ================== ITERATIONS ==================
    std::string num_ranks_str = std::to_string(size);
    outputFilePath = mesh_file.substr(0, mesh_file.find_last_of('.')) + "/metrics/np" + num_ranks_str + '_' + bc + "_" + ic + itMethString + ".txt";
    std::ofstream metricsFile(outputFilePath); // Open file for writing metrics
    const std::string HEADER = "Iteration norm_z max_grad min_grad";
    if (rank == 0) {
        // Write header to file
        metricsFile << HEADER << '\n';
    }
    
    // Now the parallelized section for all ranks
    bool converged = false;
    
    for (int iter = 0; iter < maxIterations && !converged; ++iter) {
        if (rank == 0) std::cout << "-------------- Iteration " << iter + 1 << " --------------" << std::endl;
        
        bool localConverged = eikonal->updateSolution(rank, size);
        MPI_Allreduce(&localConverged, &converged, 1, MPI_C_BOOL, MPI_LAND, MPI_COMM_WORLD);
        
        if (rank == 0){
            auto metrics = eikonal->getMetrics();

            // Write metrics to file
            metricsFile << iter + 1 << ' '; // Write iteration number
            
            metricsFile << std::get<0>(metrics) << ' '
                        << std::get<1>(metrics) << ' '
                        << std::get<2>(metrics) << '\n';
                        
            if (converged) {
                std::cout << "Solution converged after " << iter + 1 << " iterations." << std::endl;

                outputFilePath = mesh_file.substr(0, mesh_file.find_last_of('.')) + "/solutions/" + bc + "_" + ic + itMethString + ".vtk";
  
                mesh.addScalarField(w, mesh_file, outputFilePath);
            }
        }
    }

    if (!converged && rank == 0) {
        std::cout << "Solution did not converge within the maximum number of iterations." << std::endl;
    }

    MPI_Finalize();

    return 0;
};
