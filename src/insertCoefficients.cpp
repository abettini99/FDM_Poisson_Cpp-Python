#include "insertCoefficients.hpp"

#include "defines_standard.hpp" 
#include "defines_eigen.hpp" 
#include "mesh.hpp"

void insertCoefficient( i64 idx, i64 ii, i64 ij, f64 weight,
                        std::vector<  Eigen::Triplet<f64>  > &coefficients,
                        EigenDefs::Vector<f64>               &forceVector,
                        Mesh::boundaryStruct                 &boundaryConditions){
    
    i32 iimax = boundaryConditions.East.size() - 2;  /**< Internal maximum size (x-dir) */
    i32 ijmax = boundaryConditions.South.size() - 2; /**< Internal maximum size (y-dir) */
    const i32 idx1 = ii + ij*iimax;


    if     (ii==-1   ) forceVector[idx] -= weight * boundaryConditions.West[ij+1];   // Remember that ii is offset from i by 1
    else if(ii==iimax) forceVector[idx] -= weight * boundaryConditions.East[ij+1];   // Remember that ij is offset from j by 1
    else if(ij==-1   ) forceVector[idx] -= weight * boundaryConditions.South[ii+1];
    else if(ij==ijmax) forceVector[idx] -= weight * boundaryConditions.North[ii+1];
    else               coefficients.push_back(Eigen::Triplet<f64>(idx,idx1,weight));
};