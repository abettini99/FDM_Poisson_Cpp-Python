#pragma once

#include "definesStandard.hpp" 
#include "definesEigen.hpp" 

namespace Mesh{

struct boundaryStruct{
    EigenDefs::Array1D<f64> North; /**< North boundary values*/
    EigenDefs::Array1D<f64> West;  /**< West boundary values*/
    EigenDefs::Array1D<f64> South; /**< South boundary values*/
    EigenDefs::Array1D<f64> East;  /**< East boundary values*/
};

struct gridStruct{
    EigenDefs::Array1D<f64> x; /**< x grid points */
    EigenDefs::Array1D<f64> y; /**< y grid points */
};

} // namespace Mesh