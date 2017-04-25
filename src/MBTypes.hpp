/*
 * MBTypes.hpp
 *
 *  Created on: 07.03.2017
 *      Author: jgeisler
 */

#ifndef MBTYPES_HPP_
#define MBTYPES_HPP_

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <vector>


typedef Eigen::Matrix<double, 3, 1> Vec;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VecX;
typedef Eigen::Matrix<bool, Eigen::Dynamic, 1> bVecX;
typedef Eigen::Matrix<double, 3, 3> Mat;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatX;
typedef Eigen::Transform<double, 3, Eigen::Affine> Trans;

using std::string;

#endif /* MBTYPES_HPP_ */
