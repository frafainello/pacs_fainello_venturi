/*
 * Eikonal_traits.hpp
 *
 *  Created on: May 31, 2022
 *      Author: forma
 */
#ifndef EIKONAL_EIKONAL_TRAITS_HPP_
#define EIKONAL_EIKONAL_TRAITS_HPP_
// clang-format off
#include <vector>
//#include<array>
#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include <memory>
//#include <unordered_set>
//#include <unordered_map>
#include <algorithm>
//PHDIM phisical dimension of the domain
//New version to accomodate also manifolds
namespace Eikonal {
  /*!
    * \brief Traits for the Eikonal equation
    * \tparam PHDIM physical dimension of the domain
    * \tparam INTRINSIC_DIM intrinsic dimension of the domain
    * \details The traits are used to define the types used in the Eikonal equation
    * \note The traits are defined for a physical dimension PHDIM and an intrinsic dimension INTRINSIC_DIM
    * The intrinsic dimension is the dimension of the manifold on which the Eikonal equation is solved.
    * The physical dimension is the dimension of the embedding space.
    * The intrinsic dimension is always less or equal to the physical dimension and, by default is equal to the physical dimension.
    * Intrinsic dimension must be >1 
  */
  template<std::size_t PHDIM, std::size_t INTRINSIC_DIM=PHDIM>
  
  struct Eikonal_traits {
      using Index = long int; //!< Index type used to address indexes
      using Indexes = std::vector<Index>; //!< A generic vector of indexes
      using Node = Eigen::Matrix<double, PHDIM, 1, Eigen::ColMajor>;  //!< The array storing the coordinate of a node      
      using Nodes = Eigen::Matrix<double, PHDIM, Eigen::Dynamic, Eigen::ColMajor>;//!< The array storing the coordinates of the nodes
      using Values = Eigen::Matrix<double, Eigen::Dynamic, 1>;//!< The array storing values at the the nodes
      template <class T>
              using allocator=std::allocator<T>;
      using AnisotropyM = Eigen::Matrix<double, PHDIM, PHDIM>;//!< The array storing the anisotropy matrix M 
      using Element = Eigen::Matrix<Index, INTRINSIC_DIM+1,1>;//!< The array storing the indexes of the nodes of an element
      using Elements = Eigen::Matrix<Index, INTRINSIC_DIM+1, Eigen::Dynamic>;//<! The array storing the indexes of the nodes of all elements
      //! The array storing the indexes of the nodes describing an elements of co-dimension 2: edges of the mesh if intrinsic dimension is 3, not used if intrinsic dimension is 2
      using Ridges = Eigen::Matrix<Index, INTRINSIC_DIM-1, Eigen::Dynamic>;
      //! The array storing the indexes of the ridges of the elements of the mesh: 6 if intrinsic dimension is 3, not used if intrinsic dimension is 2
      using Elements_ridges = Eigen::Matrix<Index, 3*(INTRINSIC_DIM - 1), Eigen::Dynamic>;
      //! The array storing the indexes of the nodes describing an element of codimension 1: faces if implicit dimension is 3, edges if intrinsic dimension is 2
      using Faces = Eigen::Matrix<Index, INTRINSIC_DIM, Eigen::Dynamic>; // INTRINSIC_DIM+1
      //! The array storing the indexes of the faces (elements of codimension 1) of the elements of the mesh: 4 per element 
      //! if intrinsic dimension is 3, 3 if intrinsic dimension is 2
      using Elements_faces = Eigen::Matrix<Index, INTRINSIC_DIM+1, Eigen::Dynamic>;
      using Vector = Eigen::Matrix<Index, Eigen::Dynamic,1>; //!< A generic vector of double
      using ListOfVectors = std::vector<Vector>; //!< A list of vectors
      using setData = std::vector<bool>;//!< A vector of boolean used as a mask
      using BooleanMatrix=Eigen::SparseMatrix<long int,Eigen::ColMajor>;//!< A sparse matrix of boolean

  };

}


#endif /* EIKONAL_EIKONAL_TRAITS_HPP_ */