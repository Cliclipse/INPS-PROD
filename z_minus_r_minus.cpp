#include "z_minus_r_minus.h"
#include "pr akito/Hermit.h"

#include <cmath>
#include <stdexcept>

arma::rowvec ZMinus(const arma::rowvec& zValues, int n_z, double b_z)
{
	if (n_z < 0)
	{
		throw std::invalid_argument("n_z must be non-negative");
	}

	if (b_z <= 0.0)
	{
		throw std::invalid_argument("b_z must be positive");
	}

	double factorial = std::tgamma(static_cast<double>(n_z) + 1.0);
	double coeff = (1.0 / std::sqrt(b_z)) * (1.0 / (std::sqrt(std::pow(2.0, n_z) * std::sqrt(std::numbers::pi) * factorial)));

	arma::rowvec expo_coeff = arma::exp(-arma::square(zValues) / (2.0 * std::pow(b_z, 2.0)));

    arma::rowvec hermiteVals = poly.calcHermite(n_z, zValues / b_z);

	return coeff * expo_coeff % hermiteVals;
}

arma::rowvec RMinus(double r_perp, double b_perp, int m, int n)
{
    if (n < 0)
    {
        throw std::invalid_argument("n must be non-negative");
    }

    double factorial1 = std::tgamma(static_cast<double>(n) + 1.0);
    double factorial2 = std::tgamma(static_cast<double>(n) + static_cast<double>(std::abs(m)) + 1.0);
    double coeff = (1.0 / b_perp * std::sqrt(std::numbers::pi)) * std::sqrt(factorial1 / factorial2) * (std::pow(r_perp / b_perp, std::abs(m)));

    arma::rowvec expo_coeff = arma::exp(-std::pow(r_perp, 2.0) / (2.0 * std::pow(b_perp, 2.0)));

    arma::rowvec laguerreVals = poly.calcLaguerre(m, n, std::pow(r_perp, 2.0) / std::pow(b_perp, 2.0));
}
