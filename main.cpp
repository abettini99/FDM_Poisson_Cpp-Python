/************************************************************************************************************************
 * In these two header files, we include:
 *    - Typedefs:   used to create an additional name for another data type, but does not create a new type.
 *    - Type alias: same as typedefs, but allows for templates, e.g. Vector<f64>, Vector<i32>.
 ************************************************************************************************************************/
#include "defines_standard.hpp"
#include "defines_eigen.hpp"     // arguably, we do not really need this.
#include "mesh.hpp"
#include "insertCoefficients.hpp"

#include <vector>
#include <iostream>

/************************************************************************************************************************
 * Solve -div(grad(u)) = f, f = 0
 ************************************************************************************************************************/
int main(){

    // Provide parameters
    const u32 imax  = 11;                /**< #gridpoints in x */
    const u32 jmax  = 11;                /**< #gridpoints in y */
    const f64 Lx[2] = {0., 2.*EIGEN_PI}; /**< domain endpoints in x */
    const f64 Ly[2] = {0., EIGEN_PI};    /**< domain endpoints in y */

    // Calculate parameters
    const u32 n = (imax-2)*(jmax-2);     /**< sparse matrix size component (n,n), boundaries excluded */

    // Declare and initalise 2D gridpoint list
    Mesh::gridStruct grid;
    grid.x.setLinSpaced(imax, Lx[0], Lx[1]);
    grid.y.setLinSpaced(jmax, Ly[0], Ly[1]);

    // Declare and initialise boundary condition list
    Mesh::boundaryStruct boundaries;
    boundaries.North.setZero(imax);
    boundaries.West = Eigen::sin(grid.y);
    boundaries.South.setZero(imax);
    boundaries.East.setZero(jmax);

    // Declare problem matrices and vectors to solve: Au = b
    Eigen::SparseMatrix<f64> A(n, n); /**< Sparse weights matrix */
    EigenDefs::Vector<f64>   u(n);    /**< Solution vector */
    EigenDefs::Vector<f64>   b(n);    /**< Forcing vector */
    b.setZero();
    
    // Fill out sparse matrix using a list of triplets (i,j,value)
    std::vector<  Eigen::Triplet<f64>  > coefficients; /**< List of triplets to fill out sparse matrix with */
    for (i64 i=1; i<imax-1; i++){
        for (i64 j=1; j<jmax-1; j++){

            // Calculate grid spacing necessary from full grid
            f64 dx1 = grid.x[i]   - grid.x[i-1];
            f64 dx2 = grid.x[i+1] - grid.x[i];
            f64 dy1 = grid.y[j]   - grid.y[j-1];
            f64 dy2 = grid.y[j+1] - grid.y[j];

            // We do not consider boundary conditions in the sparse matrix, so we need to consider correct for index
            const i32 ii = i-1;
            const i32 ij = j-1;
            const i32 idx = ii + ij*(imax-2); // iimax = imax - 2
            
            insertCoefficient(idx,ii,  ij-1, -2/( dy1*(dy1+dy2) ),     coefficients, b, boundaries);
            insertCoefficient(idx,ii-1,ij,   -2/( dx1*(dx1+dx2) ),     coefficients, b, boundaries);
            insertCoefficient(idx,ii,  ij,    2/(dx1*dx2)+2/(dy1*dy2), coefficients, b, boundaries);
            insertCoefficient(idx,ii+1,ij,   -2/( dx2*(dx1+dx2) ),     coefficients, b, boundaries);
            insertCoefficient(idx,ii,  ij+1, -2/( dy2*(dy1+dy2) ),     coefficients, b, boundaries);
        }
    }
    A.setFromTriplets(coefficients.begin(), coefficients.end());

    // Solve problem for u (Au = b)
    // see https://eigen.tuxfamily.org/dox/group__TutorialSparse.html
    Eigen::SimplicialCholesky<Eigen::SparseMatrix<f64>> chol(A);  /**< Cholesky factorization of A */
    u = chol.solve(b);                                            // use the factorization to solve for the given right hand side
 
    return EXIT_SUCCESS;
}
