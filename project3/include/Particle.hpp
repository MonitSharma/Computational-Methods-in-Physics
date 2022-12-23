#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include <armadillo>

class Particle
{

  public:

    double q_, m_;
    arma::vec pos_, vel_;
    bool outofbounds_;


    // Constructor
    Particle(double q_in, double m_in, arma::vec pos_in, arma::vec vel_in);


    /* Only if private parameters:
    // Method that returns the charge
    int charge();

    // Method that returns the mass_
    double mass();

    // Method that returns the position vector of the particle
    arma::vec position();

    // Method that returns the velocity vector of the particle
    arma::vec velocity();
    */

};

#endif
