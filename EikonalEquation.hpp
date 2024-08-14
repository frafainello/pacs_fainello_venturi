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

    using Values = Eigen::Matrix<double, Eigen::Dynamic, 1>;
    using Gradients = typename Eigen::Matrix<double, Eigen::Dynamic, PHDIM>;
    
    EikonalEquation(const MeshData<PHDIM>& mesh, 
                    Values& w, 
                    const Eigen::SparseMatrix<double>& stiffnessMatrix, 
                    const std::vector<Eigen::Matrix<double, PHDIM, PHDIM>>& gradientCoeff,
                    const std::vector<long int>& boundaryIndices,
                    const Eigen::SparseLU<Eigen::SparseMatrix<double>>& solver,
                    double gamma = 1e-3,
                    double tol = 1e-4)
        : mesh(mesh), 
          w(w), 
          stiffnessMatrix(stiffnessMatrix), 
          gradientCoeff(gradientCoeff), 
          boundaryIndices(boundaryIndices), 
          solver(solver), 
          gamma(gamma), 
          tol(tol),

          local_w(mesh.getConnectivity().rows()),
          localStiffness(mesh.getConnectivity().rows(), mesh.getConnectivity().rows()),
          local_rhs(mesh.getConnectivity().rows()),
          grad_w(gradientCoeff.size(), PHDIM),
          norm_grad_w(gradientCoeff.size(), 1),
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
        norm_grad_w.setZero();
        rhs.setZero();

        apsc::LinearFiniteElement<PHDIM> linearFiniteElement;
        
        localNodes.setZero();
        globalNodeNumbers.resize(globalNodeNumbers.size(), 0);

    }

    void computeLocalRhs(int i) {
        
        for (int j = 0; j < mesh.getConnectivity().col(i).size(); j++) {
            local_w(j) = w(mesh.getConnectivity().col(i)(j));
        }
        
        grad_w.row(i) = apsc::computeGradient<PHDIM>(gradientCoeff[i], local_w);

        double norm_grad = grad_w.row(i).norm();
        norm_grad_w.coeffRef(i) = norm_grad; // scalar
        
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

        double coeff = computeCoeff(norm_grad);
        local_rhs = coeff * (linearFiniteElement.getLocalStiffness() * local_w); // contributo locale di ogni elemento al rhs globale
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
        for (int idx : boundaryIndices) {
            rhs(idx) = 0.0 * 1e40;
        }
        
        // Solve the lienar system
        Values z = solver.solve(rhs);

        // Update the solution
        w += z;
        std::cout << "z : " << z.norm() << std::endl;
        return (z.norm() < tol);
    }

    virtual double computeCoeff(double norm_grad, double gamma = 1e-3) = 0;

protected:
    MeshData<PHDIM> mesh;
    Values& w;
    const Eigen::SparseMatrix<double>& stiffnessMatrix;
    const std::vector<Eigen::Matrix<double, PHDIM, PHDIM>>& gradientCoeff;
    const std::vector<long int>& boundaryIndices;
    const Eigen::SparseLU<Eigen::SparseMatrix<double>>& solver;
    double gamma;
    double tol;

    // Additional members that are initialized in the constructor body (and those that are reset at the beginning of each iteration)
    Values local_w;
    Eigen::SparseMatrix<double> localStiffness;
    Values local_rhs;
    Gradients grad_w;
    Values norm_grad_w;
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

    using Values = Eigen::Matrix<double, Eigen::Dynamic, 1>;
    using Gradients = typename Eigen::Matrix<double, Eigen::Dynamic, PHDIM>;

    StandardEikonal(MeshData<PHDIM>& mesh,
                    Values& w,
                    const Eigen::SparseMatrix<double>& stiffnessMatrix,
                    const std::vector<Eigen::Matrix<double, PHDIM, PHDIM>>& gradientCoeff,
                    const std::vector<long int>& boundaryIndices,
                    const Eigen::SparseLU<Eigen::SparseMatrix<double>>& solver)
        : EikonalEquation<PHDIM>(mesh, w, stiffnessMatrix, gradientCoeff, boundaryIndices, solver) {}

    double computeCoeff(double norm_grad, double gamma = 1e-3) override {
        return (1.0 - norm_grad) / (norm_grad + gamma);
    }
};

template<std::size_t PHDIM, std::size_t INTRINSIC_DIM=PHDIM>
class PenaltyEikonal : public EikonalEquation<PHDIM> {
public:
    using Traits = Eikonal::Eikonal_traits<PHDIM, INTRINSIC_DIM>;
    using Index = typename Traits::Index;
    using Indexes = typename Traits::Indexes;
    using Nodes = typename Traits::Nodes;
    using Elements = typename Traits::Elements;

    using Values = Eigen::Matrix<double, Eigen::Dynamic, 1>;
    using Gradients = typename Eigen::Matrix<double, Eigen::Dynamic, PHDIM>;
    
    PenaltyEikonal(MeshData<PHDIM>& mesh,
                    Values& w,
                    const Eigen::SparseMatrix<double>& stiffnessMatrix,
                    const std::vector<Eigen::Matrix<double, PHDIM, PHDIM>>& gradientCoeff,
                    const std::vector<long int>& boundaryIndices,
                    const Eigen::SparseLU<Eigen::SparseMatrix<double>>& solver)
        : EikonalEquation<PHDIM>(mesh, w, stiffnessMatrix, gradientCoeff, boundaryIndices, solver) {}

    double computeCoeff(double norm_grad, double gamma = 1e-3) override {
        return (1.0 - norm_grad) / (norm_grad + gamma);
    }
};

#endif // EIKONAL_EQUATION_HPP