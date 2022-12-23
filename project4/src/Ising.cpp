#include "../include/Ising.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <assert.h>
#include <vector>
#include <cmath>
#include <random>
#include <string>
#include <iomanip>
#include <armadillo>

using namespace std;

Ising::Ising(int lattice_side_length, double T, int seed, int ordered_spin, int burn_in) {
    L_ = lattice_side_length;
    N_ = L_*L_;
    T_ = T;
    mag_per_spin_ = 0;
    epsilon_ = 0;
    tot_cycles_ = 0;
    burn_in_cycles_ = burn_in;
    burned_ = 0;
    totalenergy_ = 0;
    sample_ = 0;
    accumulatedtotalenergy_ = 0;
    accumulatedtotalmagnetization_ = 0;
    magnetisation_ = 0;
    E2_ = 0; M2_ = 0;

    //Initalise randomness with Mersenne Twister 19937 random number generator
    generator_.seed(seed);
    proposal_pdf_ = normal_distribution<double>(0.0, 1.0);
    uniform_real_ = uniform_real_distribution<double>(0.0, 1.0);
    lattice_uniform_distribution_ = uniform_int_distribution<int>(0, L_-1);
    up_or_down_spin_ = uniform_int_distribution<int>(0, 1);

    boltzmann_factors_ = calc_boltzmann_factors(T);

    if (ordered_spin != 0){
      generate_ordered_lattice(ordered_spin);
    } else {
      generate_unordered_lattice();
    }

    calc_tot_magnetisation_of_state();
    calc_energy_of_lattice_state();
}

void Ising::generate_ordered_lattice(int spin) {
    vector<vector<int>> lattice(L_, vector<int>(L_, spin));
    s_ = lattice;
}

void Ising::generate_unordered_lattice() {
     vector<vector<int>> lattice(L_, vector<int>(L_, 1));
     for (int i=0; i<L_; i++){
        for (int j=0; j<L_; j++){
            lattice[i][j] = lattice[i][j] - 2 * up_or_down_spin_(generator_); // Generates a 1 or a -1
        }
    }
    s_ = lattice;
}

void Ising::run_metropolis_MCMC(){
  int randRow, randCol, index, deltaE;
  //eps_ = 0;
  eps_cycle_.clear(); // reset each cycle
  //eps_ = (totalenergy_) /N_;
  for (int c = 0; c < N_; c++){ // one MC cycle; attempt N spin flips
    // flip random spin
    randRow = lattice_uniform_distribution_(generator_);
    randCol = lattice_uniform_distribution_(generator_);

    // examining surrounding spins to figure out index in boltzmann_factor vector
    // for computing the probability ratio. Computes the difference in energy after flipping a spin
    deltaE = s_[randRow][randCol] * 2 // The spin flip happens here (no minus).
           *(s_[randRow][(randCol - 1 + L_) % L_] // Neighbour to the left
           + s_[randRow][(randCol + 1) % L_] // Neighbour to the right
           + s_[(randRow + 1) % L_][randCol] // Neighbour above
           + s_[(randRow - 1 + L_) % L_][randCol]); // Neighbour below

    index = deltaE/4 + 2; // [4] = 8/4+2, [3] = 4/4+2, [2] = 0/4+2, [1] = -4/4+2, [0] = -8/4+2
    // Acceptance ratio
    double probability_ratio = boltzmann_factors_[index]; // w_i/w_j = exp(-beta*deltaE)
    double r = uniform_real_(generator_);

    if (r <= probability_ratio){ //abs(totalenergy_ + deltaE) < abs(totalenergy_)
      // Accept spin configuration candidate
      // Always accept for energy reducing flips
      // Set new state of system:
      s_[randRow][randCol] *= -1;
      totalenergy_ += deltaE;
      eps_cycle_.push_back(totalenergy_/N_);
      magnetisation_ += 2 * s_[randRow][randCol]; // Equation 13.7 in lectures2015 M_(i+1) = M_i + 2*s_(i+1) (= +/- 2 )
    }
  }

  // Adding the values from each cycle once sytem has reached equilibrium,
  // so it can be used to find exp values.
  burned_ += 1;
  if (burn_in_cycles_ < burned_) {
    tot_cycles_ += 1;
    epsilon_ += 1.0*totalenergy_/N_;
    mag_per_spin_ += 1.0*abs(magnetisation_)/ N_;
    accumulatedtotalenergy_ += totalenergy_; //accumulatedtotalenergy_ er sum(E_i) over alle cycles i
    accumulatedtotalmagnetization_ += abs(magnetisation_);
    M2_ += pow(magnetisation_, 2);
    E2_ += pow(totalenergy_, 2);
  }
}

void Ising::burn_in_lattice() {
  for (int i = 0; i<burn_in_cycles_; i++) {
    run_metropolis_MCMC();
  }
}

double Ising::mean(double value, int n_cycles){
  return value / n_cycles;
}

double Ising::expval_epsilon(int n_cycles){
  return mean(epsilon_, n_cycles);
}

double Ising::expval_mag_per_spin(int n_cycles){
  return mean(mag_per_spin_, n_cycles);
}

double Ising::heat_capacity(int n_cycles){
  return (1./N_)*(1./(T_*T_))*((E2_/tot_cycles_) - pow((mean(accumulatedtotalenergy_, tot_cycles_)), 2));
}

double Ising::susceptibility(int n_cycles){
  return (1./N_)*(1./T_)*(mean(M2_, n_cycles) - pow(mean(accumulatedtotalmagnetization_, n_cycles), 2));
}

void Ising::calc_energy_of_lattice_state() {
  double energy = 0;
  for (int i=0; i<L_; i++){
    for (int j=0; j<L_; j++){
      energy +=  - s_[i][j] * (s_[(i+1)%L_][j] + s_[i][(j+1)%L_]);
    }
  }
  totalenergy_ = energy;
}

void Ising::calc_tot_magnetisation_of_state(){
  int magnetisation = 0;
  for(int i=0 ; i<L_ ; i++){ //the first row will be the Lth row
    for(int j=0 ; j<L_ ; j++){ //the first column will be the Lth column
      magnetisation +=s_[i][j];
    }
  }
  magnetisation_ = magnetisation;
}

vector<double> Ising::calc_boltzmann_factors(double T){
  double beta = 1. / (T);
  vector<double> boltzmann_values;
  boltzmann_values.push_back(exp(-beta*(-8))); // 0 +1 spins
  boltzmann_values.push_back(exp(-beta*(-4))); // 1 +1 spins
  boltzmann_values.push_back(exp(-beta*(0))); // = 1? // 2 +1 spins
  boltzmann_values.push_back(exp(-beta*(4))); // 3 +1 spins
  boltzmann_values.push_back(exp(-beta*(8))); // 4 +1 spins
  return boltzmann_values;
}

void Ising::write_parameters_to_file(ofstream& ofile) {
  int width = 16;
  ofile << setw(width) << tot_cycles_;
  ofile << setw(width) << totalenergy_;
  ofile << setw(width) << magnetisation_;
  ofile << setw(width) << expval_epsilon(tot_cycles_-burn_in_cycles_);
  ofile << setw(width) << expval_mag_per_spin(tot_cycles_-burn_in_cycles_);
  ofile << setw(width) << heat_capacity(tot_cycles_-burn_in_cycles_);
  ofile << setw(width) << susceptibility(tot_cycles_-burn_in_cycles_);
  ofile << endl;
  sample_+=1;
}

void Ising::write_eps_to_file(ofstream& ofile) {
  int width = 16;
  for (int i = 0; i < eps_cycle_.size(); i++){
    ofile << setw(width) << eps_cycle_[i] << endl;
  }
  sample_+=1;
}

arma::vec Ising::sample_average_over_sampled_values(int samples) {
  arma::vec results = arma::vec(4).fill(0);

  results(0) = expval_epsilon(tot_cycles_);
  results(1) = expval_mag_per_spin(tot_cycles_);
  results(2) = heat_capacity(tot_cycles_);
  results(3) = susceptibility(tot_cycles_);
  for (int i = 0; i<samples-1; i++) {
    run_metropolis_MCMC();
    results(0) += expval_epsilon(tot_cycles_);
    results(1) += expval_mag_per_spin(tot_cycles_);
    results(2) += heat_capacity(tot_cycles_);
    results(3) += susceptibility(tot_cycles_);
  }
  return results;
}

void Ising::print(){
  cout << "Energy of lattice: " << totalenergy_ << "\n";
  cout << "Magnetisation of lattice: " << magnetisation_ << "\n";
    for (int i=0; i<L_; i++){
        for (int j=0; j<L_; j++){
          if(s_[i][j] > 0) {
            cout << "+";
          }
            cout << s_[i][j] << " ";
        }
        cout << "\n";
    }
}
