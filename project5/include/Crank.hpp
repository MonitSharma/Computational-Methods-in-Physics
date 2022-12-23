#ifndef CRANK_HPP
#define CRANK_HPP
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <chrono>

#include <armadillo>

#include <complex>
#include <cmath>

using namespace std::complex_literals; // to use imaginary number i | DEMANDS c++14!
using namespace std;
using namespace arma;

// Class that numerically solves the time evolution of a Gaussian wavepacket in a box in a 0-3 slit potential barrier setup
class Crank {
    private:
    sp_cx_mat A_, B_;
    cx_mat U_, U_empty;
    mat V_;
    cx_vec u_;
    double deltat_, h_;
    int M_, t_, t_steps_; 
    complex<double> r_;

    public:
    double v_0_;

    // Class that numerically solves the time evolution of a Gaussian wavepacket in a box in a 0-3 slit potential barrier setup
    Crank(double h, double deltat, double T, double x_c, double y_c, double sigma_x, double sigma_y, double p_x, double p_y, double v_0, int slits);

        // Performs the simulation
        cx_cube run_simulation();

        // Overload: Performs the simulation and only returning the state at the last time step (for Problem 9).
        cx_mat run_simulation(int last_slice);

            // Numerically solves the next state of the system using the Crank-Nicolson method. 
            cx_vec time_step(cx_vec u);

        // Creates a normalised Gaussian wavepacket in a box.
        cx_mat make_insert_wavepacket(int M, double h, double x_c, double y_c, double sigma_x, double sigma_y, double p_x, double p_y);

        // Constructs a column vector u based on the matrix U
        cx_vec construct_u_vec(cx_mat U, bool normalise);

        // Constructs the A and B matrices used by the time_step/Crank-Nicolsonmethod and stores them as class variables.
        void make_matrices(int M, double h, double deltat, mat V, complex<double> r);

            // Helper method for constructing the A and B matrices used by the time_step/Crank-Nicolsonmethod and stores them as class variables.
            sp_cx_mat make_matrix(complex<double> r, cx_vec d);

        // Translates matrix (i,j) index to column (k) index which have values from 1 to M-1.
        int get_k_index(int i, int j, int M);

        // Converts a column vector u to matrix form.
        cx_mat col_to_mat(cx_vec u);

        // Makes the box with a potential barrier V that traps the particle with no barriers within. Parameters fetched from class variables.
        mat make_potential_box();

        // Makes the box with a potential barrier V divided in two by a barrier with a single slit opening. Parameters fetched from class variables.
        mat make_potential_single_slit();

        // Makes the box with a potential barrier V divided in two by a barrier with double slit openings. Parameters fetched from class variables.
        mat make_potential_double_slit();

        // Makes the box with a potential barrier V divided in two by a barrier with triple slit openings. Parameters fetched from class variables.
        mat make_potential_triple_slit();  
};
#endif
