/************************************************************************************************************************
 * In these two header files, we include:
 *    - Typedefs:   used to create an additional name for another data type, but does not create a new type.
 *    - Type alias: same as typedefs, but allows for templates, e.g. Vector<f64>, Vector<i32>.
 ************************************************************************************************************************/
#include "CoreIncludes.hpp"
#include "solver/CG.hpp"
#include "solver/BiCGstab_l_.hpp"
#include "mesh/mesh.hpp"
#include "mesh/valueSource.hpp"

#include <vector>
#include <fstream>
#include <iostream>

/************************************************************************************************************************
 * Solve -div(grad(u)) = f, using FDM
 ************************************************************************************************************************/
int main(){

    //## ================== ##//
    //## Provide parameters ##//
    //## ================== ##//
    const u32 imax  = 1001;              /**< #gridpoints in x */
    const u32 jmax  = 1001;              /**< #gridpoints in y */
    const f64 Lx[2] = {0., 1.*EIGEN_PI}; /**< domain endpoints in x */
    const f64 Ly[2] = {0., 1.*EIGEN_PI}; /**< domain endpoints in y */

    //## ==================== ##//
    //## Calculate parameters ##//
    //## ==================== ##//
    const u32 n = (imax-2)*(jmax-2);     /**< sparse matrix size component (n,n), boundaries excluded */

    //## ============= ##//
    //## Problem Setup ##//
    //## ============= ##//
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
    // Internal
    std::vector<  Eigen::Triplet<f64>  > coefficients; /**< List of triplets to fill out sparse matrix with */
    for (u32 j=2; j<jmax-2; j++){
        for (u32 i=2; i<imax-2; i++){
            // Calculate internal matrix internal index relative to grid position
            u32 jj = j-1;
            u32 ii = i-1;
            u32 iimax = imax-2;

            // Calculate grid spacing necessary from full grid
            f64 dx1 = grid.x[i]   - grid.x[i-1];
            f64 dx2 = grid.x[i+1] - grid.x[i];
            f64 dy1 = grid.y[j]   - grid.y[j-1];
            f64 dy2 = grid.y[j+1] - grid.y[j];

            // We consider BCs in the sparse matrix problem, fill general pattern in A matrix.
            const u32 idx = jj*iimax + ii;
            u32 idx1;
            idx1 = (jj-1)*iimax + (ii)  ; coefficients.push_back(  Eigen::Triplet<f64>(idx,idx1, -2./( dy1*(dy1+dy2) ))  );
            idx1 = (jj)  *iimax + (ii-1); coefficients.push_back(  Eigen::Triplet<f64>(idx,idx1, -2./( dx1*(dx1+dx2) ))  );
            idx1 = (jj)  *iimax + (ii)  ; coefficients.push_back(  Eigen::Triplet<f64>(idx,idx1,  2./(dx1*dx2)+2./(dy1*dy2))  );
            idx1 = (jj)  *iimax + (ii+1); coefficients.push_back(  Eigen::Triplet<f64>(idx,idx1, -2./( dx2*(dx1+dx2) ))  );
            idx1 = (jj+1)*iimax + (ii)  ; coefficients.push_back(  Eigen::Triplet<f64>(idx,idx1, -2./( dy2*(dy1+dy2) ))  );
            
            // Fill source term in b vector.
            b[idx] = valueSource(grid.x[i], grid.y[j]);
        }
    }
    // South Boundary
    u32 j=1;
    for (u32 i=1; i<imax-1; i++){
        // Calculate internal matrix internal index relative to grid position
        u32 jj = j-1;
        u32 ii = i-1;
        u32 iimax = imax-2;

        // Calculate grid spacing necessary from full grid
        f64 dy1 = grid.y[j]   - grid.y[j-1];
        f64 dy2 = grid.y[j+1] - grid.y[j];

        const u32 idx = jj*iimax + ii;
        u32 idx1;
        idx1 = (jj-1)*iimax + (ii)  ; b[idx] -= -2./( dy1*(dy1+dy2) ) * boundaries.South[i];
        idx1 = (jj)  *iimax + (ii)  ; coefficients.push_back(  Eigen::Triplet<f64>(idx,idx1,  2./(dy1*dy2)        )  );
        idx1 = (jj+1)*iimax + (ii)  ; coefficients.push_back(  Eigen::Triplet<f64>(idx,idx1, -2./( dy2*(dy1+dy2) ))  );
    
        // b[idx] += valueSource(grid.x[i], grid.y[j]);
    }
    // North Boundary
    j=jmax-2;
    for (u32 i=1; i<imax-1; i++){
        // Calculate internal matrix internal index relative to grid position
        u32 jj = j-1;
        u32 ii = i-1;
        u32 iimax = imax-2;

        // Calculate grid spacing necessary from full grid
        f64 dy1 = grid.y[j]   - grid.y[j-1];
        f64 dy2 = grid.y[j+1] - grid.y[j];

        const u32 idx = jj*iimax + ii;
        u32 idx1;
        idx1 = (jj-1)*iimax + (ii)  ; coefficients.push_back(  Eigen::Triplet<f64>(idx,idx1, -2./( dy1*(dy1+dy2) ))  );
        idx1 = (jj)  *iimax + (ii)  ; coefficients.push_back(  Eigen::Triplet<f64>(idx,idx1,  2./(dy1*dy2)        )  );
        idx1 = (jj+1)*iimax + (ii)  ; b[idx] -= -2./( dy2*(dy1+dy2) ) * boundaries.South[i];
        
        // b[idx] += valueSource(grid.x[i], grid.y[j]);
    }
    // West Boundary
    i64 i=1;
    for (u32 j=1; j<jmax-1; j++){
        // Calculate internal matrix internal index relative to grid position
        u32 jj = j-1;
        u32 ii = i-1;
        u32 iimax = imax-2;

        // Calculate grid spacing necessary from full grid
        f64 dx1 = grid.x[i]   - grid.x[i-1];
        f64 dx2 = grid.x[i+1] - grid.x[i];

        const u32 idx = jj*iimax + ii;
        u32 idx1;
        idx1 = (jj)  *iimax + (ii-1); b[idx] -= -2./( dx1*(dx1+dx2) ) * boundaries.West[j];
        idx1 = (jj)  *iimax + (ii)  ; coefficients.push_back(  Eigen::Triplet<f64>(idx,idx1,  2./(dx1*dx2)        )  );
        idx1 = (jj)  *iimax + (ii+1); coefficients.push_back(  Eigen::Triplet<f64>(idx,idx1, -2./( dx2*(dx1+dx2) ))  );
    
        // b[idx] += valueSource(grid.x[i], grid.y[j]);
    }
    // East Boundary
    i=imax-2;
    for (u32 j=1; j<jmax-1; j++){
        // Calculate internal matrix internal index relative to grid position
        u32 jj = j-1;
        u32 ii = i-1;
        u32 iimax = imax-2;

        // Calculate grid spacing necessary from full grid
        f64 dx1 = grid.x[i]   - grid.x[i-1];
        f64 dx2 = grid.x[i+1] - grid.x[i];

        const u32 idx = jj*iimax + ii;
        u32 idx1;
        idx1 = (jj)  *iimax + (ii-1); coefficients.push_back(  Eigen::Triplet<f64>(idx,idx1, -2./( dx1*(dx1+dx2) ))  );
        idx1 = (jj)  *iimax + (ii)  ; coefficients.push_back(  Eigen::Triplet<f64>(idx,idx1,  2./(dx1*dx2)        )  );
        idx1 = (jj)  *iimax + (ii+1); b[idx] -= -2./( dx2*(dx1+dx2) ) * boundaries.East[j];
        
        // b[idx] += valueSource(grid.x[i], grid.y[j]);
    }
    // Fill out sparse matrix
    A.setFromTriplets(coefficients.begin(), coefficients.end());

    INFO_MSG("Matrix-Vector setup finished");

    //## ================ ##//
    //## Solution Routine ##//
    //## ================ ##//
    // Initialization    
    const u8 l = 8;
    u32 kappa;          /**< iterate number*/
    u32 kappaMax = 5000;    /**< max iterate allowed */
    f64 tol = 1e-15;        /**< acceptable tolerance */
    f64 err = 1/0.0;        /**< residual error */


    EigenDefs::Vector<f64> &hx0 = u;
    EigenDefs::Vector<f64> rk(n);            /**< residual vector */
    std::vector< EigenDefs::Vector<f64> > hu(l+1);
    std::vector< EigenDefs::Vector<f64> > hr(l+1);

    EigenDefs::Vector<f64> tr0(n);
    
    f64 alpha, beta, omega;                              /**< update coefficients */
    f64 rho0, rho1;
    f64 tau[l+1][l+1];  // ignore index 0, hence the +1;
    f64 sigma[l+1];     // ignore index 0;
    f64 gam;
    f64 gamma[l+1];
    f64 gammap[l+1];
    f64 gammapp[l+1];



    rk = b - A*u;
    kappa=-l;
    tr0 = rk;
    for (u8 i=0; i<=l; i++){
        hu[i].setZero(n);
        hr[i].setZero(n);
    }
    hr[0] = b - A*u;
    rho0 = 1.;
    alpha = 0.;
    omega = 1.;

do {
    kappa += l;
    rho0 = -omega*rho0;
    
    //## ---- ##//
    //## BiCG ##//
    //## ---- ##//
    for (u8 j=0; j<=l-1; j++){
        rho1 = hr[j].dot(tr0);
        beta = alpha * rho1/rho0;
        rho0 = rho1;
        for (u8 i=0; i<=j; i++){
            hu[i] = hr[i] - beta*hu[i];
        }
        hu[j+1] = A*hu[j];
        gam = hu[j+1].dot(tr0);
        alpha = rho0/gam;
        for (u8 i=0; i<=j; i++){
            hr[i] = hr[i] - alpha*hu[i+1];
        }
        hr[j+1] = A*hr[j];
        hx0 = hx0 + alpha*hu[0];
    }
        
    //## ------- ##//
    //## mod.G-S ##//
    //## ------- ##//
    sigma[1] = hr[1].dot(hr[1]);
    gammap[1] = 1/sigma[1] * hr[0].dot(hr[1]);
    for (u8 j=2; j<=l; j++){
        for (u8 i=1; i<=j-1; i++){
            tau[i][j] = 1/sigma[i] * hr[j].dot(hr[i]);
            hr[j] = hr[j] - tau[i][j]*hr[i];
        }
        sigma[j] = hr[j].dot(hr[j]);
        gammap[j] = 1/sigma[j] * hr[0].dot(hr[j]); 
    }

    gamma[l] = gammap[l];
    omega = gamma[l];

    for (u8 j=l-1; j>=1; j--){
        f64 sum = 0;
        for (u8 i=j+1; i<=l; i++) sum += tau[j][i]*gamma[i];
        gamma[j] = gammap[j] - sum;
    }
    for (u8 j=1; j<=l-1; j++){
        f64 sum = 0;
        for (u8 i=j+1; i<=l-1; i++) sum += tau[j][i]*gamma[i+1];
        gammapp[j] = gamma[j+1] + sum;
    }

    //## ------ ##//
    //## update ##//
    //## ------ ##//
    hx0 = hx0 + gamma[1]*hr[0];
    hr[0] = hr[0] - gammap[l]*hr[l];
    hu[0] = hu[0] - gamma[l]*hu[l];

    for (u8 j=1; j<=l-1; j++){
        hu[0] = hu[0] - gamma[j]*hu[j];
        hx0 = hx0 + gammapp[j]*hr[j];
        hr[0] = hr[0] - gammap[j]*hr[j];
    }

    rk = hr[0];
    err = std::sqrt( rk.dot(rk)/rk.size() );
    CHECK_FATAL_ITERERROR(kappa, err);
    INFO_MSG("kappa = %-5u err = %1.4e", kappa, err); 
    if (tol > err) break;

    } while (kappa < kappaMax); 




    // u32 kappa;          /**< iterate number*/
    // u32 kappaMax = 5000;    /**< max iterate allowed */
    // f64 tol = 1e-15;        /**< acceptable tolerance */
    // f64 err = 1/0.0;        /**< residual error */

    // EigenDefs::Vector<f64> pk(n), pkm1(n);                   /**< search/conjugate direction vector */
    // EigenDefs::Vector<f64> vk(n), vkm1(n);                   /**< search/conjugate direction vector */
    // EigenDefs::Vector<f64> hr0(n), rk(n), rkm1(n);   /**< residual vector */
    // EigenDefs::Vector<f64> s(n), t(n);   /**< residual vector */
    // f64 alpha, beta, omegak, omegakm1;                              /**< update coefficients */
    // f64 rhok, rhokm1;

    // kappa = 0;
    // rkm1 = b - A*u;
    // hr0 = rkm1;
    // rhokm1 = 1;
    // alpha = 1; 
    // omegakm1 = 1;
    // vkm1.setZero(n);
    // pkm1.setZero(n);
    // kappa++;

    // do {
    
    //     rhok = hr0.dot(rkm1);
    //     beta = (rhok/rhokm1)*(alpha/omegakm1);
    //     pk = rkm1 + beta*(pkm1 - omegakm1*vkm1);
    //     vk = A*pk;
    //     alpha = rhok/ hr0.dot(vk);
    //     s = rkm1 - alpha*vk;
    //     t = A*s;
    //     omegak = t.dot(s)/t.dot(t);
    //     u = u + alpha*pk + omegak*s;
    //     rk = s - omegak*t;

    //     err = std::sqrt( rk.dot(rk)/rk.size() );
    //     CHECK_FATAL_ITERERROR(kappa, err);
    //     INFO_MSG("kappa = %-5u err = %1.4e", kappa, err); 
    //     if (tol > err) break;

    //     kappa++;
    //     omegakm1 = omegak;
    //     rhokm1 = rhok;
    //     rkm1 = rk;
    //     vkm1 = vk;
    //     pkm1 = pk;

    // } while (kappa < kappaMax); 





    // Diagonal Preconditioner
    // - see https://diamhomes.ewi.tudelft.nl/~mvangijzen/PhDCourse_DTU/LES5/TRANSPARANTEN/les5.pdf 
    // - see Section 4.1 https://homepage.tudelft.nl/d2b4e/burgers/lin_notes.pdf
    // Eigen::SparseMatrix<f64> Mm1(n,n);          /**< Inversed Preconditioner M^-1*/
    // Mm1 = A.diagonal().asDiagonal().inverse();
    // f64 Mm1 = 1.;

    // KrylovSolver::CG solver(A);
    // solver.solve(u,b);




    // - see https://eigen.tuxfamily.org/dox/group__TutorialSparse.html
    // - see https://eigen.tuxfamily.org/dox/classEigen_1_1SparseLU.html
    // Eigen::SparseLU<Eigen::SparseMatrix<f64>> solver(A);  /**< LU factorization of A */
    // solver.analyzePattern(A); // Compute the column permutation to minimize the fill-in
    // solver.factorize(A);      // Compute the numerical factorization 
    // u = solver.solve(b);      // use the factorization to solve for the given right hand side

    //## =============== ##//
    //## Export solution ##//
    //## =============== ##//
    std::ofstream dataFile("data.bin", std::ios::out | std::ios::binary | std::ios::trunc); /**< Data output file, not using std::ios::app */ 
    dataFile.write(reinterpret_cast<const char*>(&jmax), sizeof(jmax));
    dataFile.write(reinterpret_cast<const char*>(&imax), sizeof(imax));

    u32 idx = 0;
    for (u32 j=0; j<jmax; j++){
        for (u32 i=0; i<imax; i++){
            // Plotting / postprocessing currently does not need to be done in such high precision
            f32 xPoint = (float) grid.x[i]; /**< f32 x-position */
            f32 yPoint = (float) grid.y[j]; /**< f32 y-position */
            f32 uPoint;                     /**< f32 u-solution */
            if      (i==0     ) {uPoint = (float) boundaries.West[j]; }
            else if (i==imax-1) {uPoint = (float) boundaries.East[j]; }
            else if (j==0     ) {uPoint = (float) boundaries.South[i];}
            else if (j==jmax-1) {uPoint = (float) boundaries.North[i];}
            else                {uPoint = (float) u[idx]; idx++;      }

            // INFO_MSG("idx=%d, uPoint=%f", idx, uPoint)
            dataFile.write(reinterpret_cast<const char*>(&xPoint), sizeof(xPoint));
            dataFile.write(reinterpret_cast<const char*>(&yPoint), sizeof(yPoint));
            dataFile.write(reinterpret_cast<const char*>(&uPoint), sizeof(uPoint));

        }
    }
    dataFile.close();

    INFO_MSG("Solution saved.");

    return EXIT_SUCCESS;
}
