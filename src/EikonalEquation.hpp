#ifndef EIKONAL_EQUATION_HPP
#define EIKONAL_EQUATION_HPP

#include <mpi.h>

template<std::size_t PHDIM, std::size_t INTRINSIC_DIM=PHDIM>
class EikonalEquation {

public:

    // Alias definitions for ease of use. 
    // Check Eikonal_traits.hpp for more details.
    using Traits = Eikonal::Eikonal_traits<PHDIM, INTRINSIC_DIM>;
    using Index = typename Traits::Index;
    using Indexes = typename Traits::Indexes;
    using Nodes = typename Traits::Nodes;
    using Values = typename Traits::Values;
    using AnisotropyM = typename Traits::AnisotropyM;

    // Virtual destructor to allow proper cleanu  p in derived classes
    virtual ~EikonalEquation() {}
    

    /// Distributes computation of rhs across MPI processes and updates solution.
    /// @param rank The rank of the current MPI process.
    /// @param size The total number of MPI processes.
    /// @return True if the local solution converged; otherwise, false.
    bool updateSolution(const int rank, const int size) {
        
        // Reset attributes for each update
        reset_attributes();

        // Prepare to parallelization. The last process handles the remainder
        int n = this->mesh.getNumElements();
        int local_n = n / size;
        int start = rank * local_n;
        int end = (rank == size - 1) ? n : start + local_n;
        
        // Variables for timing
        double start_time, end_time, elapsed_time;
        start_time = MPI_Wtime();

        // Compute local rhs for each element and update global rhs
        for (int i = start; i < end; i++) {
            computeLocalRhs(i);
            updateGlobalRhs(i);
        }
        
        // Sum global rhs values from all processes to rank 0
        Values global_rhs(rhs.size());
        MPI_Reduce(rhs.data(), global_rhs.data(), rhs.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        
        // Timing computation
        end_time = MPI_Wtime();
        elapsed_time = end_time - start_time;
        if (rank == 0){
            std::cout << "Time taken for rhs with " << size << " cores is " << elapsed_time << " s" << std::endl;
        }

        // Sum global grad_w values from all processes to rank 0
        Eigen::Matrix<double, Eigen::Dynamic, PHDIM> global_grad_w(grad_w.size() / PHDIM, PHDIM);
        MPI_Reduce(grad_w.data(), global_grad_w.data(), grad_w.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

        bool localConverged = true; // Initialize to true for ranks other than 0

        // rank 0 updates the solution
        if (rank == 0) { 
            rhs = global_rhs;
            grad_w = global_grad_w;
            
            // Apply boundary conditions with TGV method
            for (int idx : mesh.getBoundaryNodes()) {
                rhs(idx) = 0.0 * 1e40;
            }
            
            // Solve the linear system for z
            z = solver.solve(rhs);

            // Incremental update of the solution w
            w += z;
            std::cout << "z : " << z.norm() << std::endl;
            updateLagrangians(z);

            // Check for convergence
            localConverged = (z.norm() < tol);
        }
    
        // Broadcast the updated value of w from rank 0 to all other ranks
        MPI_Bcast(w.data(), w.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
        
        return localConverged;
    }

    /// Computes solution metrics like max/min values of z and gradient norms.
    /// This method should be called only by rank 0.
    /// @return A tuple containing the norm of z, maximum gradient norm, and minimum gradient norm.
    std::tuple<double, double, double> getMetrics() {
        // This method should be called only by rank 0
        double Z = z.norm();
        
        // Return the result as a pair
        return std::make_tuple(Z, maxGrad, minGrad);
    }


protected:
    // Constructor initializing mesh data, values, and solver
    /// @param mesh The mesh data for the problem.
    /// @param w Reference to the solution vector.
    /// @param stiffnessMatrix The global stiffness matrix.
    /// @param solver The sparse solver to be used.
    /// @param gamma Stabilization parameter.
    /// @param tol Convergence tolerance.
    EikonalEquation(const MeshData<PHDIM>& mesh, 
                    Values& w, 
                    const Eigen::SparseMatrix<double>& stiffnessMatrix, 
                    const Eigen::SparseLU<Eigen::SparseMatrix<double>>& solver,
                    double gamma = 1e-3,
                    double tol = 1e-3)
        : mesh(mesh), 
          w(w), 
          stiffnessMatrix(stiffnessMatrix), 
          solver(solver), 
          gamma(gamma), 
          tol(tol),

          z(mesh.getNumNodes()),  
          local_w(PHDIM+1),
          localStiffness(PHDIM+1, PHDIM+1),
          local_rhs(PHDIM+1),
          grad_w(this->mesh.getNumElements(), PHDIM),
          rhs(Values::Zero(mesh.getNumNodes())), 
          linearFiniteElement(),
          localNodes(PHDIM, PHDIM+1),
          globalNodeNumbers(PHDIM+1),
          M(1.0 * Eigen::Matrix<double, PHDIM, PHDIM>::Identity()),
          maxGrad(-std::numeric_limits<double>::infinity()),
          minGrad(std::numeric_limits<double>::infinity())
          {}
    
    /// Attributes related to mesh, values, and solver
    MeshData<PHDIM> mesh; ///< Mesh data for the problem.
    Values rhs; ///< Global right-hand side
    Values& w; ///< Solution
    Values z; ///< Incremental update
    const Eigen::SparseMatrix<double>& stiffnessMatrix; ///< Global stiffness matrix.
    const Eigen::SparseLU<Eigen::SparseMatrix<double>>& solver; ///< Sparse solver for the system.
    double gamma; ///< Stabilization parameter
    double tol; ///< Convergence tolerance

    // Local variables for each element of the mesh
    Values local_w; ///< Local solution vector.
    Eigen::SparseMatrix<double> localStiffness; ///< Local stiffness matrix.
    Values local_rhs; ///< Local right-hand side vector.
    Eigen::Matrix<double, Eigen::Dynamic, PHDIM> grad_w; ///< Gradient of w for each element.
    apsc::LinearFiniteElement<PHDIM> linearFiniteElement; ///< Linear finite element for local operations.
    Nodes localNodes; ///< Local nodes for each element.
    Indexes globalNodeNumbers; ///< Global node numbers corresponding to local nodes.

    // Anisotropy matrix
    AnisotropyM M; ///< Anisotropy matrix for the problem.

    /// Virtual function to compute the stiffness term for the i-th element.
    /// Must be implemented by derived classes.
    /// @param i The index of the element.
    /// @return The computed stiffness term.
    virtual Values computeStiffnessTerm(int i) = 0;
    /// Virtual function to compute the reaction term for the i-th element.
    /// Must be implemented by derived classes.
    /// @param i The index of the element.
    /// @return The computed reaction term.
    virtual Values computeReactionTerm(int i) = 0;
    /// Virtual function to update the Lagrangian multipliers for each element.
    /// Must be implemented by derived classes.
    /// @param z The incremental update to the solution.
    virtual void updateLagrangians(const Values& z) = 0;

    // Computes the anisotropic norm of a vector using the matrix M
    /// @param u The input vector.
    /// @return The anisotropic norm.   
    double anisotropicNorm(const Eigen::Matrix<double, 1, PHDIM>& u) const {
        return std::sqrt((u * M * u.transpose()).value());
    }

    // Saving for plots
    double maxGrad; ///< Maximum gradient norm.
    double minGrad; ///< Minimum gradient norm.

private:

    // Resets attributes before each update cycle
    void reset_attributes() {

        local_w.setZero();
        localStiffness.setZero();
        local_rhs.setZero();

        grad_w.setZero();
        rhs.setZero();

        apsc::LinearFiniteElement<PHDIM> linearFiniteElement;
        
        localNodes.setZero();
        globalNodeNumbers.resize(globalNodeNumbers.size(), 0);

        maxGrad = -std::numeric_limits<double>::infinity();
        minGrad = std::numeric_limits<double>::infinity();
    }

    // Computes the local right-hand side (rhs) for the i-th element
    /// @param i The index of the element.
    void computeLocalRhs(int i) {
        
        // Collect local w values
        for (int j = 0; j < PHDIM+1; j++) {
            local_w(j) = w(mesh.getConnectivity().col(i)(j));
        }
        
        // Compute local gradient
        grad_w.row(i) = apsc::computeGradient<PHDIM>(this->mesh.getGradientCoeff()[i], local_w);
        
        // Save max/min for plots
        double norm = anisotropicNorm(grad_w.row(i));

        if (norm > maxGrad) {
            maxGrad = norm;
        }
        if (norm < minGrad) {
            minGrad = norm;
        }

        // Compute stiffness term and reaction term
        Values stiffnessTerm = computeStiffnessTerm(i);
        Values reactionTerm = computeReactionTerm(i);

        // Compute local rhs
        local_rhs = stiffnessTerm + reactionTerm; 
    }

    // Updates the global right-hand side with local rhs values
    /// @param i The index of the element.
    void updateGlobalRhs(int i) {
        for (int j = 0; j < PHDIM+1; ++j) {
            // Accumulate local rhs into global rhs
            rhs(mesh.getConnectivity().col(i)(j)) += local_rhs(j);
        }
    }
};

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// Method of Laplacian Iterations (Standard Eikonal)
template<std::size_t PHDIM, std::size_t INTRINSIC_DIM=PHDIM>
class StandardEikonal : public EikonalEquation<PHDIM> {
public:
    
    using Traits = Eikonal::Eikonal_traits<PHDIM, INTRINSIC_DIM>;
    using Index = typename Traits::Index;
    using Indexes = typename Traits::Indexes;
    using Nodes = typename Traits::Nodes;
    using Values = typename Traits::Values;
    
    // The StandardEikonal uses the EikonalEquation constructor
    StandardEikonal(MeshData<PHDIM>& mesh,
                    Values& w,
                    const Eigen::SparseMatrix<double>& stiffnessMatrix,
                    const Eigen::SparseLU<Eigen::SparseMatrix<double>>& solver)
        : EikonalEquation<PHDIM>(mesh, w, stiffnessMatrix, solver) {}

    // Compute the stiffness term for the i-th element
    Values computeStiffnessTerm(int i) override {
        Eigen::Matrix<double, 1, PHDIM> grad = this->grad_w.row(i);
        double stiffnessCoeff = (1.0 - this->anisotropicNorm(grad)) / (this->anisotropicNorm(grad) + this->gamma);
        return stiffnessCoeff * (this->mesh.getLocalStiffnessMatrices()[i] * this->local_w);
    }

    // Not used in the StandardEikonal
    Values computeReactionTerm(int i) override {
        return Values::Zero(PHDIM+1);
    }

    // Not used in the StandardEikonal
    void updateLagrangians(const Values& z) override { return; }
};

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// Relaxation of the Method of Laplacian Iterations (Penalty Eikonal)
template<std::size_t PHDIM, std::size_t INTRINSIC_DIM=PHDIM>
class PenaltyEikonal : public EikonalEquation<PHDIM> {
public:
    using Traits = Eikonal::Eikonal_traits<PHDIM, INTRINSIC_DIM>;
    using Index = typename Traits::Index;
    using Indexes = typename Traits::Indexes;
    using Nodes = typename Traits::Nodes;
    using Values = typename Traits::Values;
    
    // Constructor for Penalty Eikonal
    PenaltyEikonal(MeshData<PHDIM>& mesh,
                    Values& w,
                    const Eigen::SparseMatrix<double>& stiffnessMatrix,
                    const Eigen::SparseLU<Eigen::SparseMatrix<double>>& solver,
                    double r)
        : EikonalEquation<PHDIM>(mesh, w, stiffnessMatrix, solver), r(r) {
        }

    // Compute the stiffness term for the i-th element
    Values computeStiffnessTerm(int i) override {
        Eigen::Matrix<double, 1, PHDIM> grad = this->grad_w.row(i);
        double stiffnessCoeff = (1.0 - this->anisotropicNorm(grad)) / ((1.0 + this->r) * this->anisotropicNorm(grad) + this->gamma);
        return stiffnessCoeff * (this->mesh.getLocalStiffnessMatrices()[i] * this->local_w);
    }

    // Not used in the PenaltyEikonal
    Values computeReactionTerm(int i) override {
        return Values::Zero(PHDIM+1);
    }

    // Not used in the PenaltyEikonal
    void updateLagrangians(const Values& z) override { return; }
private:
    double r;
};

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// Alternating Direction Method of Lagrange Multipliers (Lagrangian Eikonal)
template<std::size_t PHDIM, std::size_t INTRINSIC_DIM=PHDIM>
class LagrangianEikonal : public EikonalEquation<PHDIM> {
public:
    using Traits = Eikonal::Eikonal_traits<PHDIM, INTRINSIC_DIM>;
    using Index = typename Traits::Index;
    using Indexes = typename Traits::Indexes;
    using Nodes = typename Traits::Nodes;
    using Values = typename Traits::Values;

    // Constructor for Lagrangian Eikonal
    LagrangianEikonal(MeshData<PHDIM>& mesh,
                Values& w,
                const Eigen::SparseMatrix<double>& stiffnessMatrix,
                const Eigen::SparseLU<Eigen::SparseMatrix<double>>& solver,
                double r,
                Eigen::Matrix<double, Eigen::Dynamic, PHDIM>& lagrangians
                )
        : EikonalEquation<PHDIM>(mesh, w, stiffnessMatrix, solver), r(r), lagrangians(lagrangians) {}
    
    // Compute the stiffness term for the i-th element
    Values computeStiffnessTerm(int i) override{
        Eigen::Matrix<double, 1, PHDIM> grad = this->grad_w.row(i);
        Eigen::Matrix<double, 1, PHDIM> localLagrangian = lagrangians.row(i);
        double q_norm = this->anisotropicNorm(grad - localLagrangian/this->r);
        double stiffnessCoeff = (1.0 - q_norm) / ((1.0 + this->r)*q_norm + this->gamma);
        return stiffnessCoeff * (this->mesh.getLocalStiffnessMatrices()[i] * this->local_w);
    }

    // Compute the reaction term for the i-th element
    Values computeReactionTerm(int i) override {
        Eigen::Matrix<double, 1, PHDIM> grad = this->grad_w.row(i);
        Eigen::Matrix<double, 1, PHDIM> localLagrangian = lagrangians.row(i);
        double q_norm = this->anisotropicNorm(grad - localLagrangian/r);
        double reactionCoeff = (q_norm - 1.0) / ((r * (1.0+r) * q_norm) + this->gamma);
        
        Values uu(PHDIM+1);
        for (int j = 0; j < PHDIM+1; j++) {
            for (int k = 0; k < PHDIM+1; k++) {
                uu(j) += this->mesh.getLocalReactionMatrices()[i][j][k].dot(localLagrangian.transpose()); // (PHDIM x 1) x (PHDIM x 1)
            }
        }
        return reactionCoeff * uu;
    }

    // Update the Lagrangian multipliers for each element
    void updateLagrangians(const Values& z) override {
        Eigen::Matrix<double, Eigen::Dynamic, 1> q_norm(this->grad_w.rows());

        // Compute the norm of the difference between the gradient and the Lagrangian
        for (int i = 0; i < this->grad_w.rows(); ++i) {
            Eigen::Matrix<double, 1, PHDIM> diff = this->grad_w.row(i) - lagrangians.row(i) / this->r;
            q_norm(i) = this->anisotropicNorm(diff);
        }
        
        // Compute coefficients for Lagrangian and gradient updates
        Eigen::Matrix<double, Eigen::Dynamic, 1> c = (1.0 + this->r * q_norm.array()) / ((1 + this->r) * q_norm.array());
        Eigen::Matrix<double, Eigen::Dynamic, 1> lagrangianCoeff = 1 - c.array();
        Eigen::Matrix<double, Eigen::Dynamic, 1> gradzCoeff = this->r * c;
        
        // Local variables for each element
        Values local_z = Values::Zero(this->mesh.getConnectivity().rows());
        Eigen::Matrix<double, Eigen::Dynamic, PHDIM> grad_z(this->mesh.getNumElements(), PHDIM);

        // Update the Lagrangian multipliers for each element
        for(int i = 0; i < this->mesh.getGradientCoeff().size(); i++) {
            // Collect local z values
            for (int j = 0; j < PHDIM+1; j++) {
                local_z(j) = z(this->mesh.getConnectivity().col(i)(j));
            }
            
            // Compute the gradient of z
            grad_z.row(i) = apsc::computeGradient<PHDIM>(this->mesh.getGradientCoeff()[i], local_z);

            // Update the Lagrangian multipliers
            lagrangians.row(i) = lagrangianCoeff(i) * lagrangians.row(i) + gradzCoeff(i) * grad_z.row(i);
        }
    }


private:
    double r;
    Eigen::Matrix<double, Eigen::Dynamic, PHDIM> lagrangians;
};


#endif // EIKONAL_EQUATION_HPP