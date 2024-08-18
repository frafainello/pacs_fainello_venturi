#ifndef HH_LINEARFINITEELEMENT_HPP
#define HH_LINEARFINITEELEMENT_HPP
#include "gradient.hpp"
#include <Eigen/Sparse>
#include <iostream>
namespace apsc {
/*! A class for linear finite elements on a simplex
@tparam N The dimension of the space
*/
template <unsigned int N> class LinearFiniteElement {
// class for just one finite element (N=3 for 3D elements)
public:
  using Node = Eigen::Matrix<double, N, 1>;
  using Nodes = Eigen::Matrix<double, N, N + 1>; 
  using Index = long int;
  using Indexes = std::vector<Index>;

  using LocalMatrix = Eigen::Matrix<double, N + 1, N + 1>;
  using GlobalMatrix = Eigen::SparseMatrix<double>;
  using LocalVector = Eigen::Matrix<double, N + 1, 1>;
  using GlobalVector = Eigen::Matrix<double, Eigen::Dynamic, 1>;
  
  /*!
  @brief Constructor
  @details The class has only a default constructor
  You should use the update method to set the nodes
  */
  LinearFiniteElement() = default; 
  // apsc::LinearFiniteElement<3> linearFiniteElement;
  // specify N=3 in the template (parametric class on the parameter N)
  

  // Getters
  auto getNodes() const{
    return nodes_;
  }
  auto getGlobalNodeNumbers() const{
    return globalNodeNumbers_;
  }
  // auto getRefGradients() const{
  //   return refGradients_;
  // }
  auto getLocalStiffness() const{
    return localStiffness_;
  }
  auto getLocalMass() const{
    return localMass_;
  }

  /*!
  @brief Update the nodes
  */
  void update(const Nodes &nodes) {
    nodes_ = nodes;
    // Compute the measure
    measure_ = apsc::computeMeasure<N>(nodes_);
  }

  /*!
  @brief update the global node numbers
  */
  void updateGlobalNodeNumbers(const Indexes &globalNodeNumbers) {
    globalNodeNumbers_ = globalNodeNumbers;
  }

  void computeLocalIntegral() const {
    localIntegral_ = measure() / (N + 1); // measure() * localRefIntegral_;
    // std::cout << "Local integral " << localIntegral_ << std::endl;
  }

  void updateGlobalIntegrals(GlobalVector& globalIntegrals) const {
    for (auto i = 0u; i < N + 1; ++i) {
      globalIntegrals.coeffRef(globalNodeNumbers_[i]) += localIntegral_;
      }
    }

  /*!
  @brief Compute the local stiffness matrix
  @details It computes the local stiffness matrix for a linear finite element
  on an N dimensional simplex given the coordinate of the vertices, the measure
  of the simplex, and a possible constant elemental multiplicative factor.
  @param factor The constant multiplicative factor
  @return The elemental stiffness matrix
  @note It also caches the result so it can be accessed later
  */
  auto computeLocalStiffness(double factor = 1.0) const -> LocalMatrix const & {
    localStiffness_ = apsc::computeLocalStiffness<N>(nodes_, measure_, factor);
    return localStiffness_;
  }

  auto computeLocalReaction(double factor = 1.0) const -> std::vector<std::vector<Eigen::Matrix<double, N, 1>>> const & {
    localReaction_ = apsc::computeLocalReaction<N>(nodes_, measure_, factor);
    return localReaction_;
  }
  
  /*!
  @brief Compute the local mass matrix
  @details It computes the (consistent) mass matrix for a linear finite element
  on an N dimensional simplex given the measure of the simplex, and a possible
  constant elemental multiplicative factor.
  @param factor The constant multiplicative factor
  @return The elemental mass matrix
  @note It also caches the result so it can be accessed later
  */
  auto computeLocalMass(double factor = 1.0) const -> LocalMatrix const & {
    localMass_ = apsc::computeLocalMass<N>(measure_, factor);
    return localMass_;
  }

  auto computeGradientCoeff() const -> Eigen::Matrix<double, N, N> const {
    return apsc::computeGradCoeff<N>(nodes_);
  }

  auto computeGradient(const LocalVector& values) const {
    return apsc::computeGradient<N>(computeGradientCoeff(), values);
  }
  
  /*!
  @brief Get the measure of the simplex
  */
  auto measure() const -> double { return measure_; }
  /*!
  @brief Compute lumped mass matrix
  @details It computes the lumped mass matrix for a linear finite element
  on an N dimensional simplex given the measure of the simplex, and a possible
  constant elemental multiplicative factor.
  @param factor The constant multiplicative factor
  @return The lumped mass matrix
  @note It also caches the result so it can be accessed later
  */
  auto computeLumpedMass(double factor = 1.0) const -> LocalMatrix const & {
    localMass_ = apsc::computeLocalLumpedMass<N>(measure_, factor);
    return localMass_;
  }
  /*! @brief Compute the mass matrix using a source term
  @details It computes the mass matrix for a linear finite element with an
  additional source term given by its nodal values. \f[ M_{ij} = \int_{K}
  \phi_i(x) \phi_j(x) (\sum_{k=1}^{N+1} \phi_k(x) s_k) \, dx \f] where \f$
  \phi_i \f$ are the basis functions of the linear finite element.
  @param source The source term
  @return The mass matrix
  @note It also caches the result so it can be accessed later
  */
  auto computeMassWithSource(const LocalVector &source) const
      -> LocalMatrix const & {
    localMass_ = apsc::computeLocalMassWithSource<N>(measure_, source);
    return localMass_;
  }
  /*!
  @brief Get the chached stiffness matrix
  */
  auto localStiffness() const -> LocalMatrix const & { return localStiffness_; }
  /*!
  @brief Get the chached reaction matrix
  */
  auto localReaction() const -> LocalMatrix const & { return localReaction_; }
  /*!
  @brief Get the chached mass matrix
  */
  auto localMass() const -> LocalMatrix const & { return localMass_; }
  /*!
  // @brief Get the gradients of the element
  // */
  // auto localGradient() const -> Nodes const & { return localGradient_; }
  /*!
  @brief Update global stiffness matrix with the cached local stiffness matrix
  @param globalStiffnessMatrix The global stiffness matrix
  @note It is assumed that the global stiffness matrix is already correctly
  initialized! Being a sparse matrix it means that the non-zero elements are
  already set.
  */
  void updateGlobalStiffnessMatrix(GlobalMatrix &globalStiffnessMatrix) const {
    for (auto i = 0u; i < N + 1; ++i) {
      for (auto j = 0u; j < N + 1; ++j) {
        globalStiffnessMatrix.coeffRef(globalNodeNumbers_[i],
                                       globalNodeNumbers_[j]) +=
            localStiffness_(i, j);
      }
    }
  }

  void updateGlobalReactionMatrix(std::vector<std::vector<Eigen::Matrix<double, N, 1>>> &globalReactionMatrix) const {
    for (auto i = 0u; i < N + 1; ++i) {
      for (auto j = 0u; j < N + 1; ++j) {
        globalReactionMatrix[globalNodeNumbers_[i]][globalNodeNumbers_[j]] += localReaction_[i][j];
      }
    }
  }

  /*!
  @brief Update global mass matrix with the cached local mass matrix
  @param globalMassMatrix The global mass matrix
  @note It is assumed that the global mass matrix is already correctly
  initialized! Being a sparse matrix it means that the non-zero elements are
  already set.
  */
  void updateGlobalMassMatrix(GlobalMatrix &globalMassMatrix) const {
    for (auto i = 0u; i < N + 1; ++i) {
      for (auto j = 0u; j < N + 1; ++j) {
        globalMassMatrix.coeffRef(globalNodeNumbers_[i],
                                  globalNodeNumbers_[j]) += localMass_(i, j);
      }
    }
  }
  /*!
  @brief Impose Dirichlet boundary condition on stiffness/mass matrix
  @param globalMatrix The global matrix.
  @param indices A vector of indices corresponding to boundary nodes.
  @note It is assumed that the global matrix is already correctly initialized and
  that the diagonal elements at the specified indices exist.
  */
void updateMatrixWithDirichletBoundary(GlobalMatrix &globalMatrix, const std::vector<long int> &indices) {
    constexpr double TGV = 1e40;
    for (auto idx : indices) {
        // Check if the diagonal element at position (idx, idx) exists
        if (globalMatrix.coeffRef(idx, idx) != 0.0) {
            // Element exists, update its value
            globalMatrix.coeffRef(idx, idx) = TGV;
        } else {
            // Element does not exist, insert a new one
            globalMatrix.insert(idx, idx) = TGV;
        }
    }
}

private:
  double measure_ = 0.0;
  Nodes nodes_;
  Indexes globalNodeNumbers_;
  mutable double localIntegral_;

  mutable LocalMatrix localStiffness_;
  mutable std::vector<std::vector<Eigen::Matrix<double, N, 1>>> localReaction_;
  mutable LocalMatrix localMass_;
  
  const double localRefIntegral_ = 1.0/apsc::factorial<N+1>();

};
// namespace apsc
};
#endif