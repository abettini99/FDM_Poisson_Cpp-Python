#pragma once

#include "CoreIncludes.hpp"

/************************************************************************************************************************ 
 *  @brief Any mesh-related functions/classes are represented in this namespace.
 * 
 *  @details
 *  This namespace serves to identify any-and-all operations related to the mesh. In the FEM-sense, this can include
 *  element spaces, elements, assembly, etc.
 ************************************************************************************************************************/
namespace Mesh{

/**< Simplistic structure representing the mesh's boundary conditions */
struct boundaryStruct{
    EigenDefs::Array1D<f64> North; /**< North boundary values*/
    EigenDefs::Array1D<f64> West;  /**< West boundary values*/
    EigenDefs::Array1D<f64> South; /**< South boundary values*/
    EigenDefs::Array1D<f64> East;  /**< East boundary values*/
};

/**< Simplistic structure representing the mesh's gridpoints. */
struct gridStruct{
    EigenDefs::Array1D<f64> x; /**< x grid points */
    EigenDefs::Array1D<f64> y; /**< y grid points */
};

} // namespace Mesh