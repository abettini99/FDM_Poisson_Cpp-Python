#include "insertCoefficients.hpp"

#include "defines_standard.hpp" 
#include "defines_eigen.hpp" 
#include "mesh.hpp"

void insertCoefficient( i64 idx, i64 i, i64 j, f64 weight,
                        std::vector<  Eigen::Triplet<f64>  > coefficients,
                        EigenDefs::Vector<f64>               forceVector,
                        Mesh::boundaryStruct                 boundaryConditions){

};