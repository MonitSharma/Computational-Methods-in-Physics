#ifndef ISING_HPP
#define ISING_HPP
#include <random>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <armadillo>

using namespace std;

class Ising {
    public:
        int L_, N_, seed, magnetisation_,sample_, tot_cycles_;
        double T_, totalenergy_, accumulatedtotalenergy_, accumulatedtotalmagnetization_,  epsilon_, mag_per_spin_, burn_in_cycles_;
        double E2_, M2_, burned_;
        vector<double> eps_cycle_;
        mt19937 generator_;
        vector<double> boltzmann_factors_;
        vector<vector<int>> s_;
        normal_distribution<double> proposal_pdf_;
        uniform_int_distribution<int> lattice_uniform_distribution_, up_or_down_spin_;
        uniform_real_distribution<double> uniform_real_;

      // A Class that models a 2D model of a ferromagnetic lattice using the Ising model
      Ising(int lattice_side_length, double T, int seed, int ordered_spin, int burn_in);

      //Generates lattices in ordered and random states
      void generate_unordered_lattice();
      void generate_ordered_lattice(int spin);

      // The meat of the program, performes a Monte Carlo cycle
      void run_metropolis_MCMC();

      // Burns in the lattice to equilibrium state
      void burn_in_lattice();

      // Calculates the current energy of the system and updates the totalenergy_ variable
      void calc_energy_of_lattice_state();

      // Calculates the current magnetisation of the system and updates the magnetisation_ variable
      void calc_tot_magnetisation_of_state();

      // Pre-calculates the possible results of the boltzmann factors
      vector<double> calc_boltzmann_factors(double T);

      // Returns mean of input
      double mean(double value, int n_cycles);

      // Calculates and returns the expectation value epsilon (of energy per spin)
      double expval_epsilon(int n_cycles);

      // Calculates and returns the expectation value of magnetisation per spin
      double expval_mag_per_spin(int n_cycles);

      // Calculates and returns the current heat capacity of the system
      double heat_capacity(int n_cycles);

      // Calculates and returns the current magnetic suseptibility of the system.
      double susceptibility(int n_cycles);

      // Writes parameters of the system state to file
      void write_parameters_to_file(ofstream& ofile);

      // Samples and returns some parameters of the current system state
      arma::vec sample_average_over_sampled_values(int samples);

      // Prints the lattice, E and M to terminal
      void print();

      // Write eps per sample to file
      void write_eps_to_file(ofstream& ofile);

};
#endif
