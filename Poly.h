#ifndef POLY_H
#define POLY_H

#include "armadillo-code-15.0.x/armadillo-code-15.0.x/include/armadillo"

/**
 * @brief Utility class for generating Hermite and Laguerre polynomial tables.
 *
 * The class precomputes polynomial values on supplied grids so subsequent
 * lookups are O(1) without redundant recurrence evaluations.
 */
class Poly
{
public:
    int ordreMax;        
    int nombreValeur;    

    arma::mat hermiteTable; 
    arma::cube laguerreTable;

    /**
     * @brief Compute Hermite polynomials up to order n on the provided points.
     * @param n highest order to tabulate (inclusive).
     * @param zVals vector of evaluation points.
     */
    void calcHermite(int n , const arma::vec& zVals);

    /**
     * @brief Compute generalized Laguerre polynomials for multiple (m,n) pairs.
     * @param mCount number of alpha values (rows) to generate.
     * @param nCount number of polynomial orders per alpha.
     * @param zVals vector of evaluation points.
     */
    void calcLaguerre(int mCount , int nCount , const arma::vec& zVals);

    /**
     * @brief Debug helper that prints a dense matrix to stdout.
     */
    void static printMatrix(arma::mat mat);

private :
    arma::mat matFact1; // Ca aussi c'est indep de n donc j'en fais des attributs
    arma::mat matFact2;

    arma::mat calcSliceN(int n , arma::mat sliceNMoins1 , arma::mat sliceNMoins2 , arma::mat mat1 , arma::mat mat2); 

};

#endif
