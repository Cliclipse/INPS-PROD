#ifndef BASIS_H
#define BASIS_H

#include "armadillo-code-15.0.x/armadillo-code-15.0.x/include/armadillo"

class Basis
{
public:
    Basis(double br, double bz, int N, double Q);

    arma::vec rPart(const arma::vec& r_perps, int m, int n) const;
    arma::vec zPart(const arma::vec& zValues, int n_z) const;

    int mMax;
    arma::ivec nMax;
    arma::imat n_zMax;

private:
    double br_;
    double bz_;
    int N_;
    double Q_;

    static double compute_Nu(unsigned int i, double N, double Q);
    int compute_mMax(double N, double Q);
    arma::ivec compute_nMax();
    arma::imat compute_nzMax(double N, double Q);
    static arma::vec associatedLaguerre(int n, int alpha, const arma::vec& x);
    static arma::vec hermitePolynomial(int order, const arma::vec& x);
};

#endif // BASIS_H
