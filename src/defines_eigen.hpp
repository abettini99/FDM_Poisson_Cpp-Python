#include "defines_standard.hpp"
#include "Eigen/Sparse"

// ------------ //
// matrix types //
// ------------ // 
/** column-major sparse matrix with elements of type T, e.g. i32, f64, i8, etc... */
template<typename T> using ei_sparseMatrix = Eigen::SparseMatrix<T>;
