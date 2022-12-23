#include "../include/Particle.hpp"
#include "../include/PenningTrap.hpp"

#include <armadillo>

//q = particle_[0].q_;


// Constructor
PenningTrap::PenningTrap(double B0_in, double V0_in, double d_in) {
  B0_ = B0_in; // definer disse
  V0_ = V0_in;
  E_ = V0_in;
  pertrubation = false;
  particle_interactions_ = true;
  d_ = d_in;
  sense_ = 1.0; // A sensible threshold for how close particles can get along one axis while not close on another axis
}

// Add a particle to the trap
void PenningTrap::add_particle(Particle p_in) {
  particles_.push_back(p_in);
}

// External electric field at point r=(x,y,z)
arma::vec PenningTrap::external_E_field(arma::vec r) {
  double v0d;
  if (pertrubation) {v0d = E_/pow(d_, 2);
  } else {
    v0d = V0_/pow(d_, 2);
  }
  arma::vec E_field = arma::vec(3).fill(0);
  E_field(0) = r(0)*v0d;
  E_field(1) = r(1)*v0d;
  E_field(2) = -2*r(2)*v0d;
  return E_field;

}

// External magnetic field at point r=(x,y,z)
arma::vec PenningTrap::external_B_field(arma::vec r) {
  arma::vec B_field = arma::vec(3).fill(0);
  B_field(2) = B0_;
  return B_field;

}

// Force on particle_i from particle_j, ignoring the magnetic forces
arma::vec PenningTrap::force_particle(int i, int j) {
  double ke = 1.38935333 * pow (10, 5);
  double qj = particles_[j].q_;
  arma::vec ipos = particles_old_state_[i].pos_;
  arma::vec jpos = particles_old_state_[j].pos_;
  arma::vec result = ipos - jpos;

  arma::vec out = ke*qj*result/pow(abs(result), 3);
  double r_temp = sqrt(pow(result(0), 2) + pow(result(1), 2) + pow(result(2), 2));
  for (int i = 0; i < 3; ++i) {
    if (std::isnan(out(i))) {out(i) = 0.0;} // Deals with division by zero resulting in nan and other numerical errors
    // Set force to 0 if below threshold for how close particles can get along a plane for sensible results but not if they are actually close
    if (abs(result(i)) < sense_ && r_temp > 2*sense_) {out(i) = 0.0;}
  }

  return out;
}

// The total force on particle_i from the external fields
arma::vec PenningTrap::total_force_external(int i) {
  // Lorentz force
  double q = particles_[i].q_;
  arma::vec v = particles_[i].vel_;
  arma::vec E = external_E_field(particles_[i].pos_);
  arma::vec B = external_B_field(particles_[i].pos_);
  return q*E + q*cross(v,B);

}

// The total force on particle_i from the other particles
arma::vec PenningTrap::total_force_particles(int i) {
  arma::vec total_force_internal = arma::vec(3).fill(0.);
  for(int j=0 ; j < particles_.size(); j++){

     if (i!= j) {
      total_force_internal += force_particle(i, j);
    }
  }
  return total_force_internal;


}

// The total force on particle_i from both external fields and other particles
arma::vec PenningTrap::total_force(int i) {
  arma::vec Fexternal = total_force_external(i);
  arma::vec Finternal;
  if (particle_interactions_) {
    Finternal = total_force_particles(i);
    return Fexternal+Finternal;
  } else {
    return Fexternal;
  }
}

// Evolve the system one time step (dt) using Runge-Kutta 4th order
void PenningTrap::evolve_RK4(double dt) {
  particles_old_state_ = particles_;
  for (int i = 0; i < particles_.size(); i++){
    arma::vec r = particles_[i].pos_;
    arma::vec v = particles_[i].vel_;
    arma::vec a = arma::vec(3).fill(0.);
    double m = particles_[i].m_;

    if (sqrt(pow(particles_[i].pos_(0), 2) + pow(particles_[i].pos_(1), 2) +pow(particles_[i].pos_(2), 2)) > d_) {
      particles_[i].outofbounds_ = true; // Particle is now considered out of bounds.
    } else {
      a = total_force(i)/m;
    }
    //std::cout << "hmm " << i << std::endl;
    //std::cout << a[1] << std::endl;
    // 1
    arma::vec k1r = dt * v;
    arma::vec k1v = dt * a;

    // 2
    particles_[i].pos_ = r + 0.5*k1r;
    particles_[i].vel_ = v + 0.5*k1v;

    if (!particles_[i].outofbounds_) {a = total_force(i)/m;}

    arma::vec k2r = dt * particles_[i].vel_;
    arma::vec k2v = dt * a;

    // 3
    particles_[i].pos_ = r + 0.5*k2r;
    particles_[i].vel_ = v + 0.5*k2v;
    if (!particles_[i].outofbounds_) {a = total_force(i)/m;}

    arma::vec k3r = dt * particles_[i].vel_;
    arma::vec k3v = dt * a;

    // 4
    particles_[i].pos_ = r + k3r;
    particles_[i].vel_ = v + k3v;
    if (!particles_[i].outofbounds_) {a = total_force(i)/m;}

    arma::vec k4r = dt * particles_[i].vel_;
    arma::vec k4v = dt * a;

    // 5
    particles_[i].pos_ = r + (1./6) * (k1r + 2 * k2r + 2 * k3r + k4r);
    particles_[i].vel_ = v + (1./6) * (k1v + 2 * k2v + 2 * k3v + k4v);
  }
}

// Evolve the system one time step (dt) using Euler-Cromer
void PenningTrap::evolve_Euler_Cromer(double dt) {
  particles_old_state_ = particles_;
  for (int i = 0; i < particles_.size(); i++){
    arma::vec r = particles_[i].pos_;
    arma::vec v = particles_[i].vel_;
    double m = particles_[i].m_;
    arma::vec a = arma::vec(3).fill(0.);

    if (sqrt(pow(particles_[i].pos_(0), 2) + pow(particles_[i].pos_(1), 2) + pow(particles_[i].pos_(2), 2)) > d_) {
      particles_[i].outofbounds_ = true;
    } else {
      a = total_force(i)/m;
    }

    v = v + (dt * a);
    r = r + dt*v;

    particles_[i].vel_ = v;
    particles_[i].pos_ = r;
  }
}
  // Count how many particles is inside the d-region
double PenningTrap::particles_inside() {
  double count = 0;
  for (int i = 0; i < particles_.size(); i++) {
    if (!particles_[i].outofbounds_) {
      count += 1;
    }
   }
  return count;
 }
