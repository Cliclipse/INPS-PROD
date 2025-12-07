#ifndef POLY_H
#define POLY_H

#include "armadillo-code-15.0.x/armadillo-code-15.0.x/include/armadillo"

class Poly
{
public:
    int ordreMax;        
    int nombreValeur;    

    arma::mat hermiteTable; 
    arma::cube laguerreTable;
    
    void calcHermite(int n , const arma::vec& zVals);
    void calcLaguerre(int mCount , int nCount , const arma::vec& zVals);
    void static printMatrix(arma::mat mat);
    arma::mat calcSliceN(int n , arma::mat sliceNMoins1 , arma::mat sliceNMoins2); 
    arma::vec hermite(int order) const;
    arma::vec getHermiteRow(int order) const;
    arma::vec laguerre(int mOrder, int nOrder) const;


};

#endif
