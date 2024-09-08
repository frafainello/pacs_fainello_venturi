#ifndef EIKONAL_EQUATION_HPP
#define EIKONAL_EQUATION_HPP

#include <mpi.h>

template<std::size_t PHDIM, std::size_t INTRINSIC_DIM=PHDIM>
class EikonalEquation {
public:

    using Traits = Eikonal::Eikonal_traits<PHDIM, INTRINSIC_DIM>;
    using Index = typename Traits::Index;
    using Indexes = typename Traits::Indexes;
    using Nodes = typename Traits::Nodes;
    using Elements = typename Traits::Elements;
    using Values = typename Traits::Values;
    
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
          local_w(mesh.getConnectivity().rows()), // PHDIM+1
          localStiffness(mesh.getConnectivity().rows(), mesh.getConnectivity().rows()), // PHDIM+1 x PHDIM+1
          local_rhs(mesh.getConnectivity().rows()), // PHDIM+1
          grad_w(this->mesh.getGradientCoeff().size(), PHDIM), // mesh.getNumElements() x PHDIM
        //   norm_grad_w(gradientCoeff.size(), 1),
          rhs(Values::Zero(stiffnessMatrix.rows())), 
          linearFiniteElement(),
          localNodes(PHDIM, PHDIM+1),
          globalNodeNumbers(PHDIM+1),
          M(1.0 * Eigen::Matrix<double, PHDIM, PHDIM>::Identity())  // Initialize M as 2 * Identity matrix
          {}
    
    virtual ~EikonalEquation() {}

    void reset_attributes() {

        local_w.setZero();
        localStiffness.setZero();
        local_rhs.setZero();

        grad_w.setZero();
        // norm_grad_w.setZero();
        rhs.setZero();

        apsc::LinearFiniteElement<PHDIM> linearFiniteElement;
        
        localNodes.setZero();
        globalNodeNumbers.resize(globalNodeNumbers.size(), 0);

    }

    // Function to compute the anisotropic norm
    double anisotropicNorm(const Eigen::Matrix<double, 1, PHDIM>& u) const {
        // Compute the anisotropic norm using the matrix M
        return std::sqrt((u * M * u.transpose()).value());
    }

    // void computeLocalReactionAndStiffness(int i) {
    //     for (auto k = 0; k < PHDIM+1; ++k) // node numbering
    //         {
    //         globalNodeNumbers[k] = mesh.getConnectivity()(k, i);
    //         for (auto j = 0; j < PHDIM; ++j) {// local node coordinates
    //             localNodes(j, k) = mesh.getNodes()(j, mesh.getConnectivity()(k, i));
    //         }
    //     }
        
    //     // Compute local nodes and update global node numbers
    //     linearFiniteElement.update(localNodes);
    //     linearFiniteElement.updateGlobalNodeNumbers(globalNodeNumbers);
        
    //     // Compute the local stiffness matrix and update the global stiffness matrix
    //     linearFiniteElement.computeLocalStiffness();
    //     linearFiniteElement.computeLocalReaction();
    //     // linearFiniteElement.computeLocalMass();
    // }

    void computeLocalRhs(int i) {
        
        // compute local w

        for (int j = 0; j < mesh.getConnectivity().col(i).size(); j++) {
            local_w(j) = w(mesh.getConnectivity().col(i)(j));
        }
        
        // compute local gradient
        grad_w.row(i) = apsc::computeGradient<PHDIM>(this->mesh.getGradientCoeff()[i], local_w);
        // double norm_grad = grad_w.row(i).norm();
        // norm_grad_w.coeffRef(i) = grad.norm(); // scalar
        
        Values stiffnessTerm = computeStiffnessTerm(i);
        Values reactionTerm = computeReactionTerm(i);

        local_rhs = stiffnessTerm + reactionTerm; 
    }

    void updateGlobalRhs(int i) {
        for (int j = 0; j < mesh.getConnectivity().col(i).size(); ++j) {
            rhs(mesh.getConnectivity().col(i)(j)) += local_rhs(j);
        }
    }

    bool updateSolution(const int rank, const int size) {
        
        reset_attributes();
        int n = this->mesh.getGradientCoeff().size();

        int local_n = n / size;
        int start = rank * local_n;
        int end = (rank == size - 1) ? n : start + local_n; // Last process handles the remainder
        
        double start_time, end_time, elapsed_time; // Variables for timing
        start_time = MPI_Wtime();
        for (int i = start; i < end; i++) {
            computeLocalRhs(i);
            updateGlobalRhs(i);
        }
        
        Values global_rhs(rhs.size());
        MPI_Reduce(rhs.data(), global_rhs.data(), rhs.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        end_time = MPI_Wtime();
        elapsed_time = end_time - start_time;
        if (rank == 0) std::cout << "Time taken for rhs with " << size << " cores is " << elapsed_time << " s" << std::endl;

        Eigen::Matrix<double, Eigen::Dynamic, PHDIM> global_grad_w(grad_w.size() / PHDIM, PHDIM);
        MPI_Reduce(grad_w.data(), global_grad_w.data(), grad_w.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

        bool localConverged = true; // Initialize to true for ranks other than 0

        if (rank == 0) { // Only the master process handles the solution update
            rhs = global_rhs;
            grad_w = global_grad_w;
            // Continue with the rest of the updateSolution as before
            for (int idx : mesh.getBoundaryNodes()) {
                rhs(idx) = 0.0 * 1e40;
            }
            
            z = solver.solve(rhs);
            w += z;
            std::cout << "z : " << z.norm() << std::endl;
            updateLagrangians(z);

            localConverged = (z.norm() < tol);
        }
    
        // Broadcast the updated value of w from rank 0 to all other ranks
        MPI_Bcast(w.data(), w.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
        // std::cout << "Rank:" << rank << " Maximum value of w: " << w.maxCoeff() << std::endl;
        return localConverged;
    }

    // Method to return the maximum and minimum of the norm of the residual z
    std::tuple<double, double, double, double> getMetrics() {
        // This method should be called only by rank 0
        double maxZ = z.maxCoeff();
        double minZ = z.minCoeff();
        
        // Compute norms for each row of grad_w
        double maxGrad = -std::numeric_limits<double>::infinity();
        double minGrad = std::numeric_limits<double>::infinity();

        for (int i = 0; i < grad_w.rows(); ++i) {
            Eigen::Matrix<double, 1, PHDIM> gradientRow = grad_w.row(i);
            double norm = anisotropicNorm(gradientRow);

            if (norm > maxGrad) {
                maxGrad = norm;
            }
            if (norm < minGrad) {
                minGrad = norm;
            }
    }
        
        // Return the result as a pair
        return std::make_tuple(maxZ, minZ, maxGrad, minGrad);
    }

    virtual Values computeStiffnessTerm(int i) = 0;
    virtual Values computeReactionTerm(int i) = 0;
    virtual void updateLagrangians(const Values& z) = 0;

protected:
    MeshData<PHDIM> mesh;
    Values& w;
    Values z;
    const Eigen::SparseMatrix<double>& stiffnessMatrix;
    const Eigen::SparseLU<Eigen::SparseMatrix<double>>& solver;
    double gamma;
    double tol;

    // Additional members that are initialized in the constructor body (and those that are reset at the beginning of each iteration)
    Values local_w;
    Eigen::SparseMatrix<double> localStiffness;
    Values local_rhs;
    Eigen::Matrix<double, Eigen::Dynamic, PHDIM> grad_w;
    // Values norm_grad_w;
    Values rhs;
    apsc::LinearFiniteElement<PHDIM> linearFiniteElement;
    Nodes localNodes;
    Indexes globalNodeNumbers;
    Eigen::Matrix<double, PHDIM, PHDIM> M;
};

template<std::size_t PHDIM, std::size_t INTRINSIC_DIM=PHDIM>
class StandardEikonal : public EikonalEquation<PHDIM> {
public:
    
    using Traits = Eikonal::Eikonal_traits<PHDIM, INTRINSIC_DIM>;
    using Index = typename Traits::Index;
    using Indexes = typename Traits::Indexes;
    using Nodes = typename Traits::Nodes;
    using Elements = typename Traits::Elements;
    using Values = typename Traits::Values;
    
    StandardEikonal(MeshData<PHDIM>& mesh,
                    Values& w,
                    const Eigen::SparseMatrix<double>& stiffnessMatrix,
                    const Eigen::SparseLU<Eigen::SparseMatrix<double>>& solver)
        : EikonalEquation<PHDIM>(mesh, w, stiffnessMatrix, solver) {}

    Values computeStiffnessTerm(int i) override {
        Eigen::Matrix<double, 1, PHDIM> grad = this->grad_w.row(i);
        // double stiffnessCoeff = (1.0 - grad.norm()) / (grad.norm() + this->gamma);
        double stiffnessCoeff = (1.0 - this->anisotropicNorm(grad)) / (this->anisotropicNorm(grad) + this->gamma);
        return stiffnessCoeff * (this->mesh.getLocalStiffnessMatrices()[i] * this->local_w);
    }
    Values computeReactionTerm(int i) override {
        return Values::Zero(PHDIM+1);
    }
    void updateLagrangians(const Values& z) override { return; }
};

template<std::size_t PHDIM, std::size_t INTRINSIC_DIM=PHDIM>
class PenaltyEikonal : public EikonalEquation<PHDIM> {
public:
    using Traits = Eikonal::Eikonal_traits<PHDIM, INTRINSIC_DIM>;
    using Index = typename Traits::Index;
    using Indexes = typename Traits::Indexes;
    using Nodes = typename Traits::Nodes;
    using Elements = typename Traits::Elements;
    using Values = typename Traits::Values;
    
    PenaltyEikonal(MeshData<PHDIM>& mesh,
                    Values& w,
                    const Eigen::SparseMatrix<double>& stiffnessMatrix,
                    const Eigen::SparseLU<Eigen::SparseMatrix<double>>& solver,
                    double r)
        : EikonalEquation<PHDIM>(mesh, w, stiffnessMatrix, solver), r(r) {
            std::cout << "r_penalty: " << r << std::endl;
        }

    Values computeStiffnessTerm(int i) override {
        Eigen::Matrix<double, 1, PHDIM> grad = this->grad_w.row(i);
        // double stiffnessCoeff = (1.0 - grad.norm()) / ((1.0 + this->r) * grad.norm() + this->gamma);
        double stiffnessCoeff = (1.0 - this->anisotropicNorm(grad)) / ((1.0 + this->r) * this->anisotropicNorm(grad) + this->gamma);
        return stiffnessCoeff * (this->mesh.getLocalStiffnessMatrices()[i] * this->local_w);
    }
    Values computeReactionTerm(int i) override {
        return Values::Zero(PHDIM+1);
    }
    void updateLagrangians(const Values& z) override { return; }
private:
    double r;
};

template<std::size_t PHDIM, std::size_t INTRINSIC_DIM=PHDIM>
class LagrangianEikonal : public EikonalEquation<PHDIM> {
public:
    using Traits = Eikonal::Eikonal_traits<PHDIM, INTRINSIC_DIM>;
    using Index = typename Traits::Index;
    using Indexes = typename Traits::Indexes;
    using Nodes = typename Traits::Nodes;
    using Elements = typename Traits::Elements;
    using Values = typename Traits::Values;

    LagrangianEikonal(MeshData<PHDIM>& mesh,
                Values& w,
                const Eigen::SparseMatrix<double>& stiffnessMatrix,
                const Eigen::SparseLU<Eigen::SparseMatrix<double>>& solver,
                double r,
                Eigen::Matrix<double, Eigen::Dynamic, PHDIM>& lagrangians
                )
        : EikonalEquation<PHDIM>(mesh, w, stiffnessMatrix, solver), r(r), lagrangians(lagrangians) {}
        
    Values computeStiffnessTerm(int i) override{
        Eigen::Matrix<double, 1, PHDIM> grad = this->grad_w.row(i);
        Eigen::Matrix<double, 1, PHDIM> localLagrangian = lagrangians.row(i);
        // double q_norm = (grad - localLagrangian/this->r).norm();
        double q_norm = this->anisotropicNorm(grad - localLagrangian/this->r);
        double stiffnessCoeff = (1.0 - q_norm) / ((1.0 + this->r)*q_norm + this->gamma);
        return stiffnessCoeff * (this->mesh.getLocalStiffnessMatrices()[i] * this->local_w);
    }
    Values computeReactionTerm(int i) override {
        Eigen::Matrix<double, 1, PHDIM> grad = this->grad_w.row(i);
        Eigen::Matrix<double, 1, PHDIM> localLagrangian = lagrangians.row(i);
        double q_norm = this->anisotropicNorm(grad - localLagrangian/r);
        double reactionCoeff = (q_norm - 1.0) / ((r * (1.0+r) * q_norm) + this->gamma);
        
        Values uu(PHDIM+1);
        for (int j = 0; j < PHDIM+1; j++) {
            for (int k = 0; k < PHDIM+1; k++) {
                uu(j) += this->mesh.getLocalReactionMatrices()[i][j][k].dot(localLagrangian.transpose()); // (PHDIMx1) x (PHDIMx1)
            }
        }
        return reactionCoeff * uu;
    }

    void updateLagrangians(const Values& z) override {
        // Eigen::Matrix<double, Eigen::Dynamic, 1> q_norm = (this->grad_w - lagrangians/r).rowwise().norm();
        // With Anisotropic norm
        Eigen::Matrix<double, Eigen::Dynamic, 1> q_norm(this->grad_w.rows());

        // Compute the anisotropic norm for each row of the matrix ////
        for (int i = 0; i < this->grad_w.rows(); ++i) {
            Eigen::Matrix<double, 1, PHDIM> diff = this->grad_w.row(i) - lagrangians.row(i) / this->r;
            q_norm(i) = this->anisotropicNorm(diff);
        }
        ///////

        // std::cout << "q_norm: \n" << q_norm << std::endl;
        
        Eigen::Matrix<double, Eigen::Dynamic, 1> c = (1.0 + this->r * q_norm.array()) / ((1 + this->r) * q_norm.array());
        // std::cout << "c: \n" << c << std::endl;
        
        // Eigen::Matrix<double, Eigen::Dynamic, 1> lagrangianCoeff = (q_norm.array() - 1.0) / ((1.0 + this->r)*q_norm.array());
        Eigen::Matrix<double, Eigen::Dynamic, 1> lagrangianCoeff = 1 - c.array();
        // std::cout << "lagrangianCoeff: \n" << lagrangianCoeff << std::endl;
          
        Eigen::Matrix<double, Eigen::Dynamic, 1> gradzCoeff = this->r * c;
        // std::cout << "gradzCoeff: \n" << gradzCoeff << std::endl;
        
        Values local_z = Values::Zero(this->mesh.getConnectivity().rows());
        Eigen::Matrix<double, Eigen::Dynamic, PHDIM> grad_z(this->mesh.getNumElements(), PHDIM);

        for(int i = 0; i < this->mesh.getGradientCoeff().size(); i++) {
            for (int j = 0; j < this->mesh.getConnectivity().col(i).size(); j++) {
                local_z(j) = z(this->mesh.getConnectivity().col(i)(j));
                }
            grad_z.row(i) = apsc::computeGradient<PHDIM>(this->mesh.getGradientCoeff()[i], local_z);

            lagrangians.row(i) = lagrangianCoeff(i)*lagrangians.row(i) + gradzCoeff(i)*grad_z.row(i);
        }
        // for (int i = 0; i < this->mesh.getNumElements(); i++) {
        //     // std::cout << "lagrangians: " << lagrangians.row(i) << std::endl;
        // }
    }

private:
    // Eigen::SparseMatrix<double> reactionMatrix;
    double r;
    Eigen::Matrix<double, Eigen::Dynamic, PHDIM> lagrangians;
};


#endif // EIKONAL_EQUATION_HPP