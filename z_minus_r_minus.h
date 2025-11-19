#ifndef Z_MINUS_H
#define Z_MINUS_H

#include "armadillo-code-15.0.x/armadillo-code-15.0.x/include/armadillo"

/**
 * @brief Computes the normalized wavefunction for a quantum harmonic oscillator in the z-direction.
 * 
 * @param zValues A row vector containing the z-coordinate values at which to evaluate the wavefunction
 * @param n_z The quantum number (energy level) of the harmonic oscillator state (must be non-negative)
 * @param b_z The characteristic length scale parameter of the harmonic oscillator (must be positive)
 * 
 * @return arma::rowvec A row vector containing the normalized wavefunction values at each z-coordinate
 * 
 * @throws std::invalid_argument if n_z is negative
 * @throws std::invalid_argument if b_z is not positive
 */
arma::rowvec ZMinus(const arma::rowvec& zValues, int n_z, double b_z);

/**
 * @brief Computes the normalized radial wavefunction for a quantum harmonic oscillator in the perpendicular direction.
 * 
 * @param r_perp The radial coordinate value at which to evaluate the wavefunction
 * @param z_perp The axial coordinate value (not used in this function but included for completeness)
 * @param b_perp The characteristic length scale parameter of the harmonic oscillator (must be positive)
 * @param m The azimuthal quantum number (can be negative, zero, or positive)
 * @param n The radial quantum number (must be non-negative)
 * 
 * @return arma::rowvec A row vector containing the normalized radial wavefunction values at the given radial coordinate
 *
 * @throws std::invalid_argument if n is negative
 */
arma::rowvec RMinus(double r_perp, double z_perp, double b_perp, int m, int n);

#endif // Z_MINUS_H
