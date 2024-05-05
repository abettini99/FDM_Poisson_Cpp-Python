/************************************************************************************************************************
 * In these two header files, we include:
 *    - Typedefs:   used to create an additional name for another data type, but does not create a new type.
 *    - Type alias: same as typedefs, but allows for templates, e.g. sparseMatrix<f64>, sparseMatrix<i32>.
 ************************************************************************************************************************/
#include "defines_standard.hpp"
#include "defines_eigen.hpp" // arguably, we do not really need this. 

int main(){

    // Provide parameters
    u32 imax = 401; /**< #gridpoints in x */
    u32 jmax = 201; /**< #gridpoints in y */
    f64 Lx[2] = {0., 2.}; /**< domain endpoints in x */
    f64 Ly[2] = {0., 1.}; /**< domain endpoints in y */

    // Calculate parameters
    u32 n, m;       /**< sparse matrix size component (n,m) */
    n = m = imax*jmax;

    return EXIT_SUCCESS;
}
