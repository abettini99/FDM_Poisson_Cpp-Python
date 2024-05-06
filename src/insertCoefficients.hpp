#pragma once

#include "defines_standard.hpp" 
#include "defines_eigen.hpp" 

#include <vector>

void insertCoefficient(i64 idx, i64 i, i64 j, f64 weight,
                       std::vector<  Eigen::Triplet<f64>  > coefficients,
                       std::vector<  Eigen::Triplet<f64>  > forceVector,
                       std::vector<  EigenDefs::Array1D<f64>  > boundaryConditions);