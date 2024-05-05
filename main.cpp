/************************************************************************************************************************
 * In these two header files, we include:
 *    - Typedefs:   used to create an additional name for another data type, but does not create a new type.
 *    - Type alias: same as typedefs, but allows for templates, e.g. Vector<f64>, Vector<i32>.
 ************************************************************************************************************************/
#include "defines_standard.hpp"
#include "defines_eigen.hpp" // arguably, we do not really need this. 

#include <vector>
#include <iostream>

int main(){

    // Provide parameters
    const u32 imax  = 21;        /**< #gridpoints in x */
    const u32 jmax  = 21;        /**< #gridpoints in y */
    const f64 Lx[2] = {0., 2.};  /**< domain endpoints in x */
    const f64 Ly[2] = {0., 1.};  /**< domain endpoints in y */

    // Calculate parameters
    const u32 n = (imax-2)*(jmax-2);     /**< sparse matrix size component (n,n), boundaries excluded */

    // Declare and initalise 2D grid
    EigenDefs::Array1D<f64> xGrid(imax); xGrid.setLinSpaced(Lx[0], Lx[1]); /**< x-domain gridpoints */
    EigenDefs::Array1D<f64> yGrid(jmax); yGrid.setLinSpaced(Ly[0], Ly[1]); /**< y-domain gridpoints*/

    // Declare problem matrices and vectors to solve: Ax = b
    Eigen::SparseMatrix<f64> A(n, n);
    EigenDefs::Vector<f64>   x(n);
    EigenDefs::Vector<f64>   b(n);
    
    // Fill out sparse matrix using a list of triplets (i,j,value)
    std::vector< Eigen::Triplet<f64> > coefficients;

    std::cout<<xGrid<<"\n";

    return EXIT_SUCCESS;
}
