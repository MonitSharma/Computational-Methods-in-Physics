#include "../include/Particle.hpp"


// Constructor
Particle::Particle(double q_in, double m_in, arma::vec pos_in, arma::vec vel_in)
{
  q_ = q_in;          // charge of particle, [e] elementary charge
  m_ = m_in;          // mass of particle in atomic mass unit u
  pos_.swap(pos_in);  // position vector r
  vel_.swap(vel_in);  // velocity vector v
  outofbounds_ = false; // Is the particle out of bounds?

}
