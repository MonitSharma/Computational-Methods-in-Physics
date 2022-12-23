#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <assert.h>
#include <vector>
#include <cmath>
#include <iomanip>
#include <omp.h>
#include <chrono>
#include <armadillo>
#include "./include/Ising.hpp"

using namespace std;
using namespace arma;
// Performs simulations based on parameter inputs
void simulator(int n_cycles, int lattice_side_length, double T, int seed, int ordered_spin, string filen, int burn_in, const int sample_rate = 100);
void problem4();
void problem5(int cycles);
void problem6(int cycles);
void epsilon_per_sample(int n_cycles,int lattice_side_length, double T, int seed, int ordered_spin, string output_file_name, int burn_in);
void problem7_8();
void analytical_2x2(double T);
void get_phase_transition_averages(double T_start, double T_end, int steps, int lattice_side_length, int seed, int ordered_spin, int burn_in, const int avg = 3);
mat phase_transitions_parallel(double T_start, double T_end, int steps, int lattice_side_length, int seed, int ordered_spin, int burn_in);
mat phase_transitions_serial(double T_start, double T_end, int steps, int lattice_side_length, int seed, int ordered_spin, int burn_in);

template <typename T> string to_string_with_precision(const T a_value, const int n = 2) {
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return out.str();
}

int main(int argc, char const *argv[]) {
  int T, L, n_cycles, ordered_spin, seed;
  string output_file_name;
   if (argc != 6)
    {
      // Get the name of the executable file
      std::string executable_name = argv[0];

      std::cerr << "Running simulations for Project 4 specified in main.cpp. To run a specific simulation use 5 parameters like so:" << std::endl;
      std::cerr << executable_name << " <temperature (float)>"
      << " <lattice side size (integer)>" << " <MCMC cycles (integer)>"
      << " <unordered lattice: use 0, ordered lattice: use -1 or 1>"
      << " <output_file_name> " << std::endl;
      problem4();
      problem5(10000);
      problem6(1e6);
      problem7_8();
      return 0; // quit program

    } else if (argc == 6) {
      T = atof(argv[1]);
      L = atoi(argv[2]);
      n_cycles = atoi(argv[3]);
      int ordered_spin = atoi(argv[4]); // 0 = unordered, ordered: -1 or 1
      output_file_name = argv[5];
      seed = 2134;
      simulator(n_cycles, L, T, seed, ordered_spin, output_file_name, 10000);
    }
  return 0;
}

void simulator(int n_cycles, int lattice_side_length, double T, int seed, int ordered_spin, string output_file_name, int burn_in, int sample_rate) {
   // ------ Output-file --------
  string filename = "datafiles/" + output_file_name;
  ofstream ofile;
  ofile.open(filename);
  // Some width and precision parameters we will use to format the output
  int width = 16;
  int prec  = 8;
  ofile << setw(width) << "Sample#";
  ofile << setw(width) << "E";
  ofile << setw(width) << "M";
  ofile << setw(width) << "<eps>";
  ofile << setw(width) << "<m>";
  ofile << setw(width) << "C_V";
  ofile << setw(width) << "Suscept.";
  ofile << endl;
  // -----------------------------
  // Run the sim
  Ising ising(lattice_side_length, T, seed, ordered_spin, burn_in);
  ising.burn_in_lattice();

  // Run MCMC cycles:
  for (int i = 0; i < ((n_cycles/sample_rate)+1); i++) {
    ising.write_parameters_to_file(ofile);
    for (int j = 0; j < sample_rate; j++) {
      ising.run_metropolis_MCMC();
    }
  }
  ofile.close();
}

void problem4() {
  // Do all the things we need for Problem 4 here
  int cycles = 100000;
  double temp = 1.0;
  simulator(cycles, 2, temp, 1337, 0, "task4.txt", 10000); //unordered
  analytical_2x2(temp);
}

void problem5(int cycles) {
  int L = 20;
  double T_1 = 1.0;
  double T_2 = 2.4;
  int seed = 3572;
  // --- Change 1e4 in textfilename if cycles change

  // T = 1.0 :
  simulator(cycles, L, T_1, seed, 1, "ncyc_1e4_L_20_T_1.0_ordered.txt", 0); // for -1 also?
  simulator(cycles, L, T_1, seed, 0, "ncyc_1e4_L_20_T_1.0_unordered.txt", 0);

  // T = 2.4
  simulator(cycles, L, T_2, seed, 1, "ncyc_1e4_L_20_T_2.4_ordered.txt", 0); // for -1 also?
  simulator(cycles, L, T_2, seed, 0, "ncyc_1e4_L_20_T_2.4_unordered.txt", 0);
}


void epsilon_per_sample(int n_cycles, int lattice_side_length, double T, int seed, int ordered_spin, string output_file_name, int burn_in) {
   // ------ Output-file --------
  string filename = "datafiles/" + output_file_name;
  ofstream ofile;
  ofile.open(filename);
  // Some width and precision parameters we will use to format the output
  int width = 16;
  int prec  = 8;
  ofile << setw(width) << "eps" << endl;
  // -----------------------------
  // Run the sim
  Ising ising(lattice_side_length, T, seed, ordered_spin, burn_in);
  ising.burn_in_lattice();

  // Run MCMC cycles:
  for (int i = 0; i < n_cycles; i++) {
    ising.write_eps_to_file(ofile);
    ising.run_metropolis_MCMC();
  }
  ofile.close();
}


void problem6(int cycles) {
  // Need its own simulator to get epsilon per sample
  int L = 20;
  double T_1 = 1.0;
  double T_2 = 2.4;
  int seed = 8276;

  // T = 1.0 :
  epsilon_per_sample(cycles, L, T_1, seed, 0, "histogram_T_1.0_unordered.txt", 10000);

  // T = 2.4
  epsilon_per_sample(cycles, L, T_2, seed, 0, "histogram_T_2.4_unordered.txt", 10000);
}

void problem7_8() {
  // Problem 7: Speedup
  double sum_average=0;
  int avg_iterations=1;
  int resolution = 10; //+1 for endpoints
  for (int i=0; i<avg_iterations; i++){
    double time_parallel = 0;
    double time_serial = 0;
    auto start = chrono::high_resolution_clock::now();
    phase_transitions_parallel(1.0, 4.0, 10, 40, 1337+i, 0, 1000);
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> p = end - start;
    time_parallel = p.count();

    start = chrono::high_resolution_clock::now();
    phase_transitions_serial(1.0, 4.0, 10, 40, 1337+i, 0, 1000);
    end = chrono::high_resolution_clock::now();
    chrono::duration<double> s = end - start;
    time_serial = s.count();
 
    sum_average += time_serial/time_parallel;
    cout << "Mesuring average speedup, at step:" << i+1 <<"/"<<avg_iterations<<endl;
  }
  cout << "Average speedup over "<<avg_iterations<<" runs, was found to be: "<< sum_average/avg_iterations<<endl;

  // Problem 8+9: Critical T
  get_phase_transition_averages(2.0, 2.6, resolution, 40, 41337, 0, 5, 10000);
  get_phase_transition_averages(2.0, 2.6, resolution, 60, 41337, 0, 5, 10000);
  get_phase_transition_averages(2.0, 2.6, resolution, 80, 41337, 0, 5, 10000);
  get_phase_transition_averages(2.0, 2.6, resolution, 100, 41337, 0, 5, 10000);

}


void analytical_2x2(double T){
  double beta = 1./T;
  double Z = 2*exp(beta*8) + 2*exp(-beta*8) + 12;
  double exp_val_epsilon = (4./Z) * (exp(-beta*8) - exp(beta*8));
  double exp_val_abs_mag = (2./Z) * (exp(beta*8) + 2);
  double eps2 = (8./Z) * (exp(-beta*8) + exp(beta*8));
  double m2 = (1./Z) * (2*exp(beta*8) + 2);
  double heat_capacity = (32./(pow(T,2)*Z))*(exp(-beta*8)+exp(beta*8) - ((2./Z)*(exp(-beta*16)+exp(beta*16)-2)));
  double susceptibility = (8./(T*Z)) * (exp(8*beta) + 1 - ((2./Z)*(exp(16*beta)+ 4*exp(8*beta)+4)));

  // Write to file
  ofstream ofile;
  // To have 2 decimals in output-filename
  std::ostringstream temp;
  temp << std::fixed << std::setprecision(1) << T;

  ofile.open("./datafiles/analytical_2x2_T=" + temp.str() +".txt");
  int width = 18;
  int prec  = 8;

  ofile << setw(2)<< "T" << setw(width) << "<eps>"
  << setw(width) << "<m>" << setw(width) << "C_v"<< setw(width) << "chi" << endl;

  ofile << setprecision(2) << scientific << T;
  ofile << setw(width) << setprecision(prec) << scientific << exp_val_epsilon;
  ofile << setw(width) << setprecision(prec) << scientific << exp_val_abs_mag;
  ofile << setw(width) << setprecision(prec) << scientific << heat_capacity;
  ofile << setw(width) << setprecision(prec) << scientific << susceptibility;
  ofile << endl;
  ofile.close();
}

void get_phase_transition_averages(double T_start, double T_end, int steps, int lattice_side_length, int seed, int ordered_spin, int avg, int burn_in) {
  double h = (T_end - T_start) / (steps); //h=2.6-20
  mat averages1 = mat((steps/2)+1, 5, fill::zeros);
  mat averages2 = mat((steps/2)+1, 5, fill::zeros);
  for (int i = 0; i<avg; i++) {
    averages1 += phase_transitions_parallel( T_start,  T_end,  steps/2,  lattice_side_length,  seed+i,  ordered_spin, burn_in);
    averages2 += phase_transitions_parallel( T_start+h,  T_end+h,  steps/2,  lattice_side_length,  seed+i,  ordered_spin, burn_in);
  }

  averages1 = averages1/avg;
  averages2 = averages2/avg;
  string filename = "datafiles/phase_transitions_parallel_T("+to_string_with_precision(T_start)+"-"+to_string_with_precision(T_end)+")_L("+to_string(lattice_side_length)+")_steps("+to_string(steps)+").txt";
  ofstream ofile;
  ofile.open(filename);
  int width = 16;
  ofile << setw(width) << "T";
  ofile << setw(width) << "<eps>";
  ofile << setw(width) << "<m>";
  ofile << setw(width) << "C_V";
  ofile << setw(width) << "Sucept.";
  ofile << endl;

  for(int row = 0; row < (steps/2)+1; row++) {
    for(int column = 0; column<5; column++) {
      ofile << setw(width) << averages1(row, column);
    }
    ofile << endl;
    for(int column = 0; column<5; column++) {
      ofile << setw(width) << averages2(row, column);
    }
    ofile << endl;
  }
  ofile.close();
}

mat phase_transitions_parallel(double T_start, double T_end, int steps, int lattice_side_length, int seed, int ordered_spin, int burn_in){
  double h = (T_end - T_start) / steps;
  int n_cycles = 100000;
  int sample_rate = 100;
  int samples = 3; // Number of samples to collect per sampling
  int spin = ordered_spin;

  //matrix to store results
  mat results = mat(steps+1, 5, fill::zeros);

  #pragma omp parallel for
  for (int i = 0; i < steps+1; i++){
    int thread_id = omp_get_thread_num();
    double T = T_start + i * h;
    Ising ising(lattice_side_length, T, seed+thread_id, spin, burn_in);
    // Burn in system
    for (int j = 0; j < burn_in+1000; j++) {
      ising.run_metropolis_MCMC();
    }

    // Collect samples
    int count = 0;
    vec params = arma::vec(4).fill(0);
    for (int j = 0; j < ((n_cycles/sample_rate)+1); j++) {
      params += ising.sample_average_over_sampled_values(samples);
      count += samples; // the params variable receives n=samples samples every loop
      for (int k = 0; k < sample_rate; k++) {
        ising.run_metropolis_MCMC();
      }
    }
    params = params/count;
    results(i,0) = T;
    results(i,1) = params(0);
    results(i,2) = params(1);
    results(i,3) = params(2);
    results(i,4) = params(3);

  }
  return results;
}

mat phase_transitions_serial(double T_start, double T_end, int steps, int lattice_side_length, int seed, int ordered_spin, int burn_in) {
  double h = (T_end - T_start) / steps;
  int n_cycles = 100000;
  int burn_int = burn_in;
  int sample_rate = 100;
  int samples = 10; // Number of samples to collect per sampling
  int spin = ordered_spin;
  //matrix to store results
  mat results = mat(steps+1, 5, fill::zeros);

  for (int i = 0; i < steps+1; i++){
    double T = T_start + i * h;
    if (T<1.5) {spin = 1;} else {spin = 0;}
    Ising ising(lattice_side_length, T, seed, spin, burn_in);
    // Burn in system
    for (int j = 0; j < burn_in+1000; j++) {
      ising.run_metropolis_MCMC();
    }

    // Collect samples
    int count = 0;
    vec params = arma::vec(4).fill(0);
    for (int j = 0; j < ((n_cycles/sample_rate)+1); j++) {
      params += ising.sample_average_over_sampled_values(samples);
      count += samples; // the params variable receives n=samples samples every loop

      for (int k = 0; k < sample_rate; k++) {
        ising.run_metropolis_MCMC();
      }
    }
    params = params/count;
    results(i,0) = T;
    results(i,1) = params(0);
    results(i,2) = params(1);
    results(i,3) = params(2);
    results(i,4) = params(3);
  }
  return results;
}
