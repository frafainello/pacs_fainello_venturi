// InitialConditions.hpp
#ifndef INITIALCONDITIONS_HPP
#define INITIALCONDITIONS_HPP

#include "Eikonal_traits.hpp"
#include <Eigen/Dense>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <iostream>

/**
 * @class InitialConditions
 * @brief This class provides methods for setting initial conditions and solving equations for a finite element mesh.
 * 
 * The class is templated to handle dimensions PHDIM and INTRINSIC_DIM, allowing for flexibility in different physical
 * and geometric configurations. It includes methods for computing initial conditions, such as those for solving the
 * heat equation.
 * 
 * @tparam PHDIM The physical dimension of the problem.
 * @tparam INTRINSIC_DIM The intrinsic dimension of the elements in the mesh (defaults to PHDIM).
 */

template<std::size_t PHDIM, std::size_t INTRINSIC_DIM=PHDIM>
class InitialConditions {

public:
    using Traits = Eikonal::Eikonal_traits<PHDIM, INTRINSIC_DIM>;
    using Values = typename Traits::Values;
    using Indexes = typename Traits::Indexes;

    /**
     * @brief Default constructor for InitialConditions.
     */
    InitialConditions() = default;

    /**
     * @brief Solves the heat equation using a sparse linear solver.
     * 
     * This method computes the solution for a heat equation by solving a linear system of equations,
     * taking into account boundary conditions provided as input. The solution is obtained using
     * LU decomposition for sparse matrices.
     * 
     * @param stiffnessMatrix The global stiffness matrix of the mesh.
     * @param Integrals Vector of integral values for each node.
     * @param forcingTerm Vector representing the forcing term in the heat equation.
     * @param boundaryIndices Indices of the nodes where boundary conditions are applied.
     * @return The computed values of the solution at each node.
     */
    Values HeatEquation(const Eigen::SparseMatrix<double>& stiffnessMatrix, 
                        const Values& Integrals,
                        const Values& forcingTerm,
                        const Indexes& boundaryIndices) const {

        // Right-hand side vector of the equation, all ones
        Values rhs = forcingTerm * Integrals;

        // Apply null Dirichlet boundary conditions
        for (int idx : boundaryIndices) {
            rhs[idx] = 0.0 * 1e40;
        }

        // Compute the solution
        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;

        solver.compute(stiffnessMatrix);
        Values solution = solver.solve(rhs);

        // Return the computed values
        return solution;
    }

private:

};

#endif // INITIALCONDITIONS_HPP