#ifndef GRADIENT_HPP
#define GRADIENT_HPP
#include <Eigen/Dense>
#include <iostream>
/// @brief Computing the gradient of the linear interpolant on a simplex
namespace apsc {

/*!
@brief Compute the coefficients of the gradient of the interpolant
@details It computes the N times N matrix that when with the vector containing
the differences in values if provides the gradient for the linear interpolant
on an N dimensional simplex

@tparam N The dimension of the space
@param vertices The vertices of the simplex
*/
// the function computeGradCoeff returns the inverse of the Jacobian matrix D, 
// which is used to compute the gradients of the shape functions in the actual
// element coordinates. This step is essential in finite element methods to 
// transform local (reference) gradients to global (actual) gradients.
template <unsigned int N>
auto computeGradCoeff(const Eigen::Matrix<double, N, N + 1> &vertices)
    -> Eigen::Matrix<double, N, N> {
  using MatrixNd = Eigen::Matrix<double, N, N>;
  // Compute the matrix of differences between vertices
  MatrixNd D;
  for (auto i = 0u; i < N; ++i) {
    D.row(i) = vertices.col(i + 1) - vertices.col(0); // D=B_trasposto
  }
  // Compute the inverse of the matrix of differences
  return D.inverse();
}
/*!
@brief Compute the local stiffness matrix if you already have the gradients of
the shape functions
@details It computes the local stiffness matrix for a linear finite element
on an N dimensional simplex given the gradients of the shape functions obtained
by computeGradCoeff, the measure of the simplex, and a possible constant
elemental multiplicative factor.
@tparam N The dimension of the space
@param gradCoeff The gradients of the shape functions computed with
computeGradCoeff
@param measure The measure of the simplex
@param factor The constant multiplicative factor
@return The elemental stiffness matrix
@note This version is useful if you have already computed the gradients of the
shape functions and you want to avoid recomputing them Note that the matrix
gradCoeff contains the gradients of the shape functions 1 to N in its columns,
excluding the gradient of shape function 0 since the latter is just the negative
sum of the other gradients
*/
template <unsigned int N>
auto computeLocalStiffnessWithGradient(
    const Eigen::Matrix<double, N, N> &gradCoeff, 
    double const &measure,
    double const &factor = 1.0) 
    -> Eigen::Matrix<double, N + 1, N + 1> {
  Eigen::Matrix<double, N, N + 1> gradCoeffExt;
  Eigen::Matrix<double, N + 1, N + 1> localStiffness;

  // Constructs gradCoeffExt, an extended gradient coefficient matrix of size N x (N + 1):
  // The first column is the negative sum of the rows of gradCoeff.
  // The remaining columns are copied directly from gradCoeff.
  gradCoeffExt.template rightCols<N>() = gradCoeff;
  gradCoeffExt.template leftCols<1>() = -gradCoeff.rowwise().sum();

  // compute the local stiffness matrix
  for (auto i = 0u; i < N + 1; ++i) {
    localStiffness(i, i) =
        factor * measure * gradCoeffExt.col(i).dot(gradCoeffExt.col(i));
    for (auto j = 0u; j < i; ++j) { // compute lower triangular part
      auto value =
          factor * measure * gradCoeffExt.col(i).dot(gradCoeffExt.col(j));
      localStiffness(i, j) = value;
      localStiffness(j, i) = value;
    }
  }
  return localStiffness;
}

/*!
@brief Compute the local stiffness matrix for a linear finite element
@details It computes the local stiffness matrix for a linear finite element
on an N dimensional simplex given the coordinate of the vertices, the measure of
the simplex, and a possible constant elemental multiplicative factor.
@tparam N The dimension of the space
@param vertices The vertices of the simplex
@param measure The measure of the simplex
@param factor The constant multiplicative factor
@return The elemental stiffness matrix
*/
template <unsigned int N>
auto computeLocalStiffness(const Eigen::Matrix<double, N, N + 1> &vertices,
                           double measure, double factor = 1.0)
    -> Eigen::Matrix<double, N + 1, N + 1> {
  return computeLocalStiffnessWithGradient<N>(computeGradCoeff<N>(vertices),
                                              measure, factor);
}

/**
@brief Compute the factorial of a number via a constexpr function
@tparam N The number to compute the factorial of
*/
template <unsigned int N> constexpr auto factorial() -> unsigned int {
  if constexpr (N == 0 || N == 1) {
    return 1;
  } else {
    return N * factorial<N - 1>();
  }
}

/*!
@brief Compute the measure of a simplex
@details It computes the measure of an N dimensional simplex given the vertices
@tparam N The dimension of the space
@param vertices The vertices of the simplex
@return The measure of the simplex
*/
template <unsigned int N>
auto computeMeasure(const Eigen::Matrix<double, N, N + 1> &vertices) -> double {
  // Compute the matrix of differences between vertices
  Eigen::Matrix<double, N, N> D;
  for (auto i = 0u; i < N; ++i) {
    D.row(i) = vertices.col(i + 1) - vertices.col(0);
  }
  // Compute the measure
  return std::abs(D.determinant()) / factorial<N>();
}

/*!
@brief Compute the mass matrix for a linear finite element
@details It computes the (consistent) mass matrix for a linear finite element
on an N dimensional simplex given the measure of the simplex, and a possible
constant elemental multiplicative factor.
@tparam N The dimension of the space
@param measure The measure of the simplex
@param factor The constant multiplicative factor
@return The elemental mass matrix
*/
template <unsigned int N>
auto computeLocalMass(double measure, double factor = 1.0)
    -> Eigen::Matrix<double, N + 1, N + 1> {

  constexpr double scale = 1.0 / ((N + 1.0) * (N + 2.0));

  Eigen::Matrix<double, N + 1, N + 1> localMass =
      Eigen::Matrix<double, N + 1, N + 1>::Constant(N + 1, N + 1,
                                                    scale * measure * factor);
  //
  // localMass.setOnes();
  localMass.diagonal().array() *= 2.0;
  // localMass *= scale * measure * factor;
  return localMass;
}

/*!
@brief Compute the local lumpedmass matrix of a linear finite element
@details It computes the local lumped mass matrix for a linear finite element
on an N dimensional simplex given the measure of the simplex, and a possible
constant elemental multiplicative factor.
@tparam N The dimension of the space
@param measure The measure of the simplex
@param factor The constant multiplicative factor
@return The elemental lumped mass matrix
*/
template <unsigned int N>
auto computeLocalLumpedMass(double measure, double factor = 1.0)
    -> Eigen::Matrix<double, N + 1, N + 1> {
  Eigen::Matrix<double, N + 1, N + 1> localLumpedMass =
      Eigen::Matrix<double, N + 1, N + 1>::Zero();
  // localLumpedMass.setZero();
  localLumpedMass.diagonal().setConstant(measure * factor / (N + 1));
  return localLumpedMass;
}
/*
  @brief Compute the  local mass matrix with a linear source term
  @details It computes the local mass matrix with a linear source term for a
  linear finite element on an N dimensional simplex given the measure of the
  simplex, and the vector of source terms In practice it computes the matrix
  given by \f[ M_{ij}=\int_{\Omega} \phi_i \phi_j + \int_{\Omega} \phi_i(x)
  \phi_j(x) \sum_{k=0}^{N} \phi_k(x)\source_k \f] where the \f$\phi_i\f$ are the
  basis functions of a linear finite element on the simplex of dimension N.
  @tparam N The dimension of the space
  @param measure The measure of the simplex
  @param source The vector of source terms
  @return The elemental mass matrix
  @note I am using the formula
  \f[
  \int_{\Omega} \phi_i^m \phi_j^n \phi_k(x)^p = \frac{m!n!p!N!}{(m+n+p+N)!}
  \f]
  */
template <unsigned int N>
auto computeLocalMassWithSource(double measure,
                                const Eigen::Matrix<double, N + 1, 1> &source)
    -> Eigen::Matrix<double, N + 1, N + 1> {
  Eigen::Matrix<double, N + 1, N + 1> localMass;
  localMass.setZero();
  // I have not found a better way to do it, sorry
  constexpr double noEqual =
      factorial<N>() / static_cast<double>(factorial<N + 3u>()); // 1! = 1
  constexpr double twoEqual = 2.0 * noEqual; // 2! = 2
  constexpr double allEqual = 6.0 * noEqual; // 3! = 6 (threeEqual)
  for (auto i = 0u; i < N + 1; ++i) {
    for (auto j = 0u; j <= i; ++j) {
      for (auto k = 0u; k < N + 1; ++k) {
        if (i == j && j == k) {
          localMass(i, j) += measure * allEqual * source(k);
        } else if (i == j || j == k || i == k) {
          localMass(i, j) += measure * twoEqual * source(k);
        } else {
          localMass(i, j) += measure * noEqual * source(k);
        }
      }
    }
  }
  return localMass.template selfadjointView<Eigen::Lower>();
}

template <unsigned int N>
auto computeLocalReactionWithGradient(
    const Eigen::Matrix<double, N, N> &gradCoeff, 
    double const &measure,
    double const &factor = 1.0) 
    -> std::vector<std::vector<Eigen::Matrix<double, N, 1>>> {

  // constexpr double scale = factorial<N>() / (N + 1.0);
  
  Eigen::Matrix<double, N, N + 1> gradCoeffExt;
  std::vector<std::vector<Eigen::Matrix<double, N, 1>>> localReaction(N + 1, std::vector<Eigen::Matrix<double, N, 1>>(N + 1));

  // Constructs gradCoeffExt, an extended gradient coefficient matrix of size N x (N + 1):
  // The first column is the negative sum of the rows of gradCoeff.
  // The remaining columns are copied directly from gradCoeff.
  gradCoeffExt.template rightCols<N>() = gradCoeff;
  gradCoeffExt.template leftCols<1>() = -gradCoeff.rowwise().sum();

  // compute the local reaction matrix
  for (auto i = 0u; i < N + 1; ++i) {
    for (auto j = 0u; j < N + 1; ++j) { // compute lower triangular part
      localReaction[i][j] = factor * measure * gradCoeffExt.col(j); // !!!!!!!!!!!!!
    }
  }
  // for (auto i = 0u; i < N + 1; ++i) {
  //   for (auto j = 0u; j < N + 1; ++j) { // compute lower triangular part
  //     std::cout << "localReaction[" << i << "][" << j << "]: " << localReaction[i][j] << std::endl;
  //     }
  //   }
  
  return localReaction;
}

template <unsigned int N>
auto computeLocalReaction(const Eigen::Matrix<double, N, N + 1> &vertices,
                           double measure, double factor = 1.0)
    -> std::vector<std::vector<Eigen::Matrix<double, N, 1>>> {
      return computeLocalReactionWithGradient<N>(computeGradCoeff<N>(vertices),
                                              measure, factor);
}


/*!
@brief Compute the gradient of the linear interpolant
@details It computes the gradient of the linear interpolant on an N dimensional
simplex given the coefficients of the gradient and the values

@tparam N The dimension of the space
@param Coefficients The coefficients of the gradient
@param values The differences in values
*/
template <unsigned int N>
auto computeGradient(const Eigen::Matrix<double, N, N> &Coefficients,
                     const Eigen::Matrix<double, N + 1, 1> &values)
    -> Eigen::Matrix<double, N, 1> {
  // Compute the differences in values
  Eigen::Matrix<double, N, 1> d_values;
  for (auto i = 0u; i < N; ++i) {
    d_values(i) = values(i + 1) - values(0);
  }
  // Compute the gradient
  return Coefficients * d_values;
}
/*!
@brief Fix the gradient so that its norm is one
@tparam N The dimension of the space
@param grad The gradient
@return The gradient with norm one
*/
template <unsigned int N>
auto fixGrad(Eigen::Matrix<double, N, 1> grad) -> Eigen::Matrix<double, N, 1> {
  auto norm = grad.norm();
  if (norm == 0) {
    return grad;
  }
  return grad / grad.norm();
}
} // namespace apsc
#endif