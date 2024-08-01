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

template<std::size_t PHDIM, std::size_t INTRINSIC_DIM=PHDIM>
class InitialConditions {

public:
    using Traits = Eikonal::Eikonal_traits<PHDIM, INTRINSIC_DIM>;
    using Values = typename Traits::Values;
    using Indexes = typename Traits::Indexes;

    InitialConditions() = default;

    Values HeatEquation(const Eigen::SparseMatrix<double>& stiffnessMatrix, 
                        const Values& Integrals,
                        const Values& forcingTerm,
                        const Indexes& boundaryIndices) const {
        // Eigen::SparseMatrix<double> modifiedStiffnessMatrix = stiffnessMatrix;

        // Right-hand side vector of the equation, all ones
        Values rhs = forcingTerm * Integrals;

        // Apply null Dirichlet boundary conditions
        for (int idx : boundaryIndices) {
            // // Zero out the row
            // for (Eigen::SparseMatrix<double>::InnerIterator it(modifiedStiffnessMatrix, idx); it; ++it) {
            //     it.valueRef() = 0.0;
            // }

            // // Zero out the column
            // for (int k = 0; k < modifiedStiffnessMatrix.outerSize(); ++k) {
            //     for (Eigen::SparseMatrix<double>::InnerIterator it(modifiedStiffnessMatrix, k); it; ++it) {
            //         if (it.row() == idx) {
            //             it.valueRef() = 0.0;
            //         }
            //     }
            // }

            // // Set diagonal element to 1
            // modifiedStiffnessMatrix.coeffRef(idx, idx) = 1.0;

            // Set corresponding rhs element to 0
            rhs[idx] = 0.0 * 1e40;
        }

        // Compute the solution
        // Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> solver;
        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
        solver.compute(stiffnessMatrix);
        // Values solution = solver.compute(modifiedStiffnessMatrix).solve(rhs);
        Values solution = solver.solve(rhs);

        // Return the computed values
        return solution;
    }

private:
    // You can add other private members here if needed

};

#endif // INITIALCONDITIONS_HPP