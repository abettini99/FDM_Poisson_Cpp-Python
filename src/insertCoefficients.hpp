#pragma once

#include "defines_standard.hpp" 
#include "defines_eigen.hpp" 
#include "mesh.hpp"

#include <vector>

void insertCoefficient( i64 idx, i64 ii, i64 ij, f64 weight,
                        std::vector<  Eigen::Triplet<f64>  > &coefficients,
                        EigenDefs::Vector<f64>               &forceVector,
                        Mesh::boundaryStruct                 &boundaryConditions);