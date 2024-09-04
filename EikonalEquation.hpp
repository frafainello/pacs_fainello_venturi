#ifndef EIKONAL_EQUATION_HPP
#define EIKONAL_EQUATION_HPP

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
                    const std::vector<Eigen::Matrix<double, PHDIM, PHDIM>>& gradientCoeff,
                    const Eigen::SparseLU<Eigen::SparseMatrix<double>>& solver,
                    double gamma = 1e-3,
                    double tol = 1e-3)
        : mesh(mesh), 
          w(w), 
          stiffnessMatrix(stiffnessMatrix), 
          gradientCoeff(gradientCoeff),
          solver(solver), 
          gamma(gamma), 
          tol(tol),

          local_w(mesh.getConnectivity().rows()), // PHDIM+1
          localStiffness(mesh.getConnectivity().rows(), mesh.getConnectivity().rows()), // PHDIM+1 x PHDIM+1
          local_rhs(mesh.getConnectivity().rows()), // PHDIM+1
          grad_w(gradientCoeff.size(), PHDIM), // mesh.getNumElements() x PHDIM
        //   norm_grad_w(gradientCoeff.size(), 1),
          rhs(Values::Zero(stiffnessMatrix.rows())), 
          linearFiniteElement(),
          localNodes(PHDIM, PHDIM+1),
          globalNodeNumbers(PHDIM+1) 
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

    void computeLocalReactionAndStiffness(int i) {
        for (auto k = 0; k < PHDIM+1; ++k) // node numbering
            {
            globalNodeNumbers[k] = mesh.getConnectivity()(k, i);
            for (auto j = 0; j < PHDIM; ++j) {// local node coordinates
                localNodes(j, k) = mesh.getNodes()(j, mesh.getConnectivity()(k, i));
            }
        }
        
        // Compute local nodes and update global node numbers
        linearFiniteElement.update(localNodes);
        linearFiniteElement.updateGlobalNodeNumbers(globalNodeNumbers);
        
        // Compute the local stiffness matrix and update the global stiffness matrix
        linearFiniteElement.computeLocalStiffness();
        linearFiniteElement.computeLocalReaction();
        // linearFiniteElement.computeLocalMass();
    }

    void computeLocalRhs(int i) {
        
        computeLocalReactionAndStiffness(i);

        // compute local w
        for (int j = 0; j < mesh.getConnectivity().col(i).size(); j++) {
            local_w(j) = w(mesh.getConnectivity().col(i)(j));
        }
        
        // compute local gradient
        grad_w.row(i) = apsc::computeGradient<PHDIM>(gradientCoeff[i], local_w);
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

    bool updateSolution() {
        
        reset_attributes();
        for(int i = 0; i < gradientCoeff.size(); i++) {
            computeLocalRhs(i);
            updateGlobalRhs(i);
        }
        
        // Update for Dirichlet BC. If == 0, set with value*TGV, where TGV=1e40
        for (int idx : mesh.getBoundaryNodes()) {
            rhs(idx) = 0.0 * 1e40;
        }
        
        // Solve the lienar system
        Values z = solver.solve(rhs);

        // Update the solution
        w += z;
        std::cout << "z : " << z.norm() << std::endl;

        updateLagrangians(z);

        return (z.norm() < tol);
    }

    virtual Values computeStiffnessTerm(int i) = 0;
    virtual Values computeReactionTerm(int i) = 0;
    virtual void updateLagrangians(const Values& z) = 0;

protected:
    MeshData<PHDIM> mesh;
    Values& w;
    const Eigen::SparseMatrix<double>& stiffnessMatrix;
    const std::vector<Eigen::Matrix<double, PHDIM, PHDIM>>& gradientCoeff;
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
                    const std::vector<Eigen::Matrix<double, PHDIM, PHDIM>>& gradientCoeff,
                    const Eigen::SparseLU<Eigen::SparseMatrix<double>>& solver)
        : EikonalEquation<PHDIM>(mesh, w, stiffnessMatrix, gradientCoeff, solver) {}

    Values computeStiffnessTerm(int i) override {
        Eigen::Matrix<double, 1, PHDIM> grad = this->grad_w.row(i);
        double norm_grad = grad.norm();
        double stiffnessCoeff = (1.0 - norm_grad) / (norm_grad + this->gamma);
        return stiffnessCoeff * (this->linearFiniteElement.getLocalStiffness() * this->local_w);
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
                    const std::vector<Eigen::Matrix<double, PHDIM, PHDIM>>& gradientCoeff,
                    const Eigen::SparseLU<Eigen::SparseMatrix<double>>& solver,
                    double r)
        : EikonalEquation<PHDIM>(mesh, w, stiffnessMatrix, gradientCoeff, solver), r(r) {}

    Values computeStiffnessTerm(int i) override {
        Eigen::Matrix<double, 1, PHDIM> grad = this->grad_w.row(i);
        double stiffnessCoeff = (1.0 - grad.norm()) / (grad.norm() + this->gamma);
        return stiffnessCoeff * (this->linearFiniteElement.getLocalStiffness() * this->local_w);
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
                const std::vector<Eigen::Matrix<double, PHDIM, PHDIM>>& gradientCoeff,
                const Eigen::SparseLU<Eigen::SparseMatrix<double>>& solver,
                double r,
                Eigen::Matrix<double, Eigen::Dynamic, PHDIM>& lagrangians
                )
        : EikonalEquation<PHDIM>(mesh, w, stiffnessMatrix, gradientCoeff, solver), r(r), lagrangians(lagrangians) {}
        
    Values computeStiffnessTerm(int i) override{
        Eigen::Matrix<double, 1, PHDIM> grad = this->grad_w.row(i);
        Eigen::Matrix<double, 1, PHDIM> localLagrangian = lagrangians.row(i);
        double q_norm = (grad - localLagrangian/this->r).norm();
        double stiffnessCoeff = (1.0 - q_norm) / ((1.0 + this->r)*(q_norm) + this->gamma);
        return stiffnessCoeff * (this->linearFiniteElement.getLocalStiffness() * this->local_w);
    }
    Values computeReactionTerm(int i) override {
        Eigen::Matrix<double, 1, PHDIM> grad = this->grad_w.row(i);
        Eigen::Matrix<double, 1, PHDIM> localLagrangian = lagrangians.row(i);
        double q_norm = (grad - localLagrangian/r).norm();
        double reactionCoeff = (q_norm - 1.0) / ((r * (1.0+r) * q_norm) + this->gamma);
        
        Values uu(PHDIM+1);
        for (int j = 0; j < PHDIM+1; j++) {
            for (int k = 0; k < PHDIM+1; k++) {
                uu(j) += this->linearFiniteElement.getLocalReaction()[j][k].dot(localLagrangian.transpose()); // (PHDIMx1) x (PHDIMx1)
            }
        }
        return reactionCoeff * uu;
    }

    void updateLagrangians(const Values& z) override {
        Eigen::Matrix<double, Eigen::Dynamic, 1> q_norm = (this->grad_w - lagrangians/r).rowwise().norm();
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

        for(int i = 0; i < this->gradientCoeff.size(); i++) {
            for (int j = 0; j < this->mesh.getConnectivity().col(i).size(); j++) {
                local_z(j) = z(this->mesh.getConnectivity().col(i)(j));
                }
            grad_z.row(i) = apsc::computeGradient<PHDIM>(this->gradientCoeff[i], local_z);

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