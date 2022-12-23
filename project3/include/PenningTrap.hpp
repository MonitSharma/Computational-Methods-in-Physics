#ifndef PENNINGTRAP_HPP
#define PENNINGTRAP_HPP

#include <vector>
#include <armadillo>

class PenningTrap
{

  public:
    std::vector<Particle> particles_, particles_old_state_;  // container of particles

    double B0_, V0_, d_, extreme_, E_, sense_;
    bool particle_interactions_, pertrubation;

    // Constructor
    PenningTrap(double B0_in, double V0_in, double d_in);

    // Add a particle to the trap
    void add_particle(Particle p_in);

    // External electric field at point r=(x,y,z)
    arma::vec external_E_field(arma::vec r);

    // External magnetic field at point r=(x,y,z)
    arma::vec external_B_field(arma::vec r);

    // Force on particle_i from particle_j
    arma::vec force_particle(int i, int j);

    // The total force on particle_i from the external fields
    arma::vec total_force_external(int i);

    // The total force on particle_i from the other particles
    arma::vec total_force_particles(int i);

    // The total force on particle_i from both external fields and other particles
    arma::vec total_force(int i);

    // Evolve the system one time step (dt) using Runge-Kutta 4th order
    void evolve_RK4(double dt);

    // Evolve the system one time step (dt) using Forward Euler
    void evolve_Euler_Cromer(double dt);

    // Count how many particles is inside the d-region
    double particles_inside();

};

#endif
