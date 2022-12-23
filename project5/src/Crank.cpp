#include "../include/Crank.hpp"
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

/// <summary>Class that numerically solves the time evolution of a Gaussian wavepacket in a box in a 0-3 slit potential barrier setup</summary>
/// <param name="h">The step size along the x and y axis.</param>
/// <param name="deltat">The temporal step.</param>
/// <param name="T">The time duration.</param>
/// <param name="x_c">The starting x-position of the center of the wavepacket.</param>
/// <param name="y_c">The starting y-position of the center of the wavepacket.</param>  
/// <param name="sigma_x">The x-dimension/size of the wavepacket.</param>
/// <param name="sigma_y">The y-dimension/size of the wavepacket.</param>
/// <param name="p_x">The starting momentum of the wavepacket in the x-direction.</param>
/// <param name="p_y">The starting momentum of the wavepacket in the y-direction.</param> 
/// <param name="v_0">The value of the potential barrier.</param> 
/// <param name="slits">The number of slits in the simulation</param>  
Crank::Crank(double h, double deltat, double T, double x_c, double y_c, double sigma_x, double sigma_y, double p_x, double p_y, double v_0, int slits=0) {
  int M = 1/h+1; //To avvoid using M as a paramater
  t_steps_ = round(T/deltat)+1;
  M_ = M;
  h_ = h;
  v_0_ = v_0;
  deltat_ = deltat;
  U_empty = cx_mat(M_, M_).fill(0); // Makes a blank canvas to be reused by the col_to_mat function
  r_ = 1i*deltat/(2*pow(h,2)); //definition of r

  if(slits==0){
    V_ = make_potential_box();
  }
  if(slits==1){
    V_ = make_potential_single_slit();
  }
  if(slits==2){
    V_ = make_potential_double_slit();
  }
  if(slits==3){
    V_ = make_potential_triple_slit(); 
  }
  //makes matrices A and B and stores them as class variables A_ and B_
  make_matrices(M_, h, deltat, V_, r_);
  U_ = make_insert_wavepacket(M_, h, x_c, y_c, sigma_x, sigma_y, p_x, p_y); //Tryin to manually correct syntax error where dims are swapped.
}

/// <summary>Performs the simulation</summary>
/// <param name="inputfile">Filename of file that specifies what simulations to run.</param>  
/// <returns>The state of the system at each time step in the simulation in the form of an armadillo complex cube, cx_cube.</returns>  
cx_cube Crank::run_simulation() {
  cx_vec u = construct_u_vec(U_, true); // calculates the initial u column vector
  cx_cube results = cx_cube(M_, M_, t_steps_);
  results.slice(0) = U_; // Add initial state to results cube
  cx_vec u_next;
  for (int i = 1; i<t_steps_; i++) {
    u_next = time_step(u);
    results.slice(i) = col_to_mat(u_next);
    u = u_next;
  }
  U_ = results.slice(t_steps_-1);
  return results;
}

/// <summary>Overload: Performs the simulation and only returning the state at the last time step (for Problem 9).</summary>
/// <param name="inputfile">Filename of file that specifies what simulations to run.</param>  
/// <returns>The state of the system at the final time step of the simulation in the form of an armadillo complex matrix, cx_mat.</returns>  
cx_mat Crank::run_simulation(int last_slice) {
  cx_vec u = construct_u_vec(U_, true); // calculates the initial u column vector
  cx_vec u_next;
  for (int i = 1; i<t_steps_; i++) {
    u_next = time_step(u);
    u = u_next;
  }
  cx_mat results = col_to_mat(u_next);
  return results;
}

/// <summary>Numerically solves the next state of the system using the Crank-Nicolson method.</summary>
/// <param name="u">Current state of the system, u.</param> 
/// <returns>The next state of the system as an Armadillo complex vector, cx_vec.</returns>  
cx_vec Crank::time_step(cx_vec u) {
  cx_vec b = B_ * u; // Solves eq 26 RHS with u(n)
  return spsolve(A_, b); // Solves eq 26 for u(n+1). Here the vector u(n) (and similarly u(n+1) is a column vector that contains the u^n_ij values for all the internal points of the xy grid at time step n
 }

/// <summary>Creates a normalised Gaussian wavepacket in a box.</summary>
/// <param name="M">The size of one size of the matrix.</param> 
/// <param name="h">The step size along the x and y axis.</param>
/// <param name="x_c">The starting x-position of the center of the wavepacket.</param>
/// <param name="y_c">The starting y-position of the center of the wavepacket.</param>  
/// <param name="sigma_x">The x-dimension/size of the wavepacket.</param>
/// <param name="sigma_y">The y-dimension/size of the wavepacket.</param>
/// <param name="p_x">The starting momentum of the wavepacket in the x-direction.</param>
/// <param name="p_y">The starting momentum of the wavepacket in the y-direction.</param>     
/// <returns>The initial state of the particle-in-a-box system as an Armadillo complex matrix, cx_mat.</returns> 
cx_mat Crank::make_insert_wavepacket(int M, double h, double x_c, double y_c, double sigma_x, double sigma_y, double p_x, double p_y){
  cx_mat U = cx_mat(M, M).fill(0); 
  double psum = 0;
  // Inserts the wavepacket and calculates normalisation factor
  int x_start = 1;
  int x_end = M-2;
  int y_start = 1;
  int y_end = M-2;
  for(int i = x_start; i<x_end; i++){
    double x = i*h;
    for(int j = y_start; j<y_end; j++){
      double y = j*h;
      complex <double> c = exp(-(pow(x - x_c ,2) / (2 * pow(sigma_x,2))) - (pow(y - y_c, 2) / (2*pow(sigma_y, 2))) + 1i * p_x * (x - x_c) + 1i * p_y * (y - y_c));
      U(i,j) = c;
      psum += real(conj(c)*c);
    }
  }
  // Normalises U to 1
  double psum2 = 0;
  for(int i = x_start; i< x_end; i++){
    for(int j = y_start; j< y_end; j++){
      complex <double> c = cx_double(1/sqrt(psum)) * cx_double(U(i,j));
      U(i,j) = c;
      psum2 += real(conj(c)*c);
    }
  }
  return U;
}

/// <summary>Constructs a column vector u based on the matrix U</summary>
/// <param name="U">Matrix to be converted.</param> 
/// <returns>The Vector as an Armadillo complex vector, cx_vec.</returns>  
cx_vec Crank::construct_u_vec(cx_mat U, bool normalise){
  int M = sqrt(U.size()); //size() gives M², assumes U to be quadratic
  cx_vec u = cx_vec(pow(M-2,2));
  for(int i = 1; i< M-1; i++){
    for(int j = 1; j< M-1; j++){
      u(get_k_index(i, j, M)) = U(i, j);
    }
  }
  return u;
}

/// <summary>Constructs the A and B matrices used by the time_step/Crank-Nicolson method and stores them as class variables.</summary>
/// <param name="M">The size of one size of the matrix.</param> 
/// <param name="h">The step size along the x and y axis.</param>
/// <param name="deltat">The time step to be used in the simulation.</param> 
/// <param name="V">The box discribed as a potential barrier</param>
/// <param name="r">A constant proportional to deltat and h, the spatial step size</param>
/// <returns>void</returns>  
void Crank::make_matrices(int M, double h, double deltat, mat V, complex<double> r){
  int mat_size = pow(M-2,2);
  cx_vec a = cx_vec(mat_size);
  cx_vec b = cx_vec(mat_size);
  for (int i = 1; i < M-1; i++){ // Excluding boundaries in V (infinity) - is that ok?
    for (int j = 1; j < M-1; j++){
      int k = get_k_index(i,j,M);
      a(k) = (1.0 + 4.0*r + 1.0i*(deltat/2*cx_double(V(i,j))));
      b(k) = (1.0 - 4.0*r - 1.0i*(deltat/2*cx_double(V(i,j))));
    }
  }
  A_ = make_matrix(-r, a);
  B_ = make_matrix(r,  b);
}

/// <summary>Helper method for constructing the A and B matrices used by the time_step/Crank-Nicolson method and stores them as class variables.</summary>
/// <param name="r">A constant proportional to deltat and h, the spatial step size</param>
/// <param name="d">The vector to be put on the diagonal as part of the Crank-Nicolson method </param>
/// <returns>The Matrix A or B as an Armadillo sparse complex matrix, sp_cx_mat</returns> 
sp_cx_mat Crank::make_matrix(complex<double> r, cx_vec d){
  int S = d.size(); // (M-2)^2
  int s = sqrt(S); // (M-2)
  sp_cx_mat M = sp_cx_mat(S, S);

  // Making diag
  for (int i = 0; i < S; i+=s){ // last index = S-s
    sp_cx_mat D(s, s);
    D.diag() = d.subvec(i, i+s-1);// diagonal
    D.diag(-1) = cx_vec(s-1).fill(r); // subdiagonal
    D.diag(1) = cx_vec(s-1).fill(r); // superdiagonal
    //submat(first_row, first_col, last_row, last_col)
    M.submat(i, i, i+s-1, i+s-1) = D; //ex: (0,0,2,2), (3,3,5,5)
  }

  // Making non-diag, non-corners
  for (int i = s; i < S; i+=s){ // last index= S-s
    sp_cx_mat ND(s,s);
    ND.diag() = cx_vec(s).fill(r); // fill diagonal with r value
    M.submat(i-s, i, i-1, i+s-1) = ND; //ex: i=s=3: (0,3,2,5) i=2s=6: (3,6,5,8)
    M.submat(i, i-s, i+s-1, i-1) = ND; //ex: i=s=3: (3,0,5,2) i=2s=6: (6,3,8,5)
  }
  return M;
}

/// <summary>Translates matrix (i,j) index to column (k) index which have values from 1 to M-1</summary>
/// <param name="i">The i index of the matrix.</param>
/// <param name="j">The j index of the matrix.</param> 
/// <param name="M">The size of one size of the matrix.</param>   
/// <returns>The column index k.</returns>  
int Crank::get_k_index(int i, int j, int M){
  return ((j - 1) * (M - 2)) + (i - 1);
}

/// <summary>Converts a column vector u to matrix form</summary>
/// <param name="u">Column vector to be converted.</param> 
/// <returns>The Matrix as an Armadillo complex matrix, cx_mat.</returns>  
cx_mat Crank::col_to_mat(cx_vec u) {
  U_empty.submat(1, 1, M_-2, M_-2) = reshape(u, M_-2, M_-2);
  return U_empty;
}

/// <summary>Makes the box with a potential barrier V that traps the particle with no barriers within. Parameters fetched from class variables.</summary>
/// <returns>A matrix describing the box V as an Armadillo matrix, mat.</returns>  
mat Crank::make_potential_box(){
  mat V = mat(M_,M_); //
  V.col(0) = vec(M_).fill(v_0_);
  V.col(M_-1) = vec(M_).fill(v_0_);
  V.row(0) = rowvec(M_).fill(v_0_);
  V.row(M_-1) = rowvec(M_).fill(v_0_);
  V.submat(1, 1, M_-2, M_-2) = mat(M_-2,M_-2).fill(0); // filling inner matrix
  return V;
}

/// <summary>Makes the box with a potential barrier V divided in two by a barrier with a single slit opening. Parameters fetched from class variables</summary>
/// <returns>A matrix describing the box V as an Armadillo matrix, mat.</returns>  
mat Crank::make_potential_single_slit(){
  mat V = make_potential_box(); //Creates a box of size M_* M_
  int center_index = (M_)*0.5; //200 * 0.5 = 100
  int x_thickness = 0.02/h_;// indices i in x direction: (0.02/0.005) + 1 = 5
  int x_start = center_index - x_thickness/2;
  int x_end = center_index + x_thickness/2;
  int aperture = (0.05/h_) ;// (0.05/0.005) + 1 = 11
  int start = center_index - aperture/2;
  int end = start + aperture + 1;
  for (int i = x_start; i<x_end+1; i++) {
    V.row(i).fill(v_0_);
    for(int j = start; j < end+1; j++) {
      V(i,j) = 0;
    }
  }
  return V;
}

/// <summary>Makes the box with a potential barrier V divided in two by a barrier with double slit openings. Parameters fetched from class variables</summary>
/// <returns>A matrix describing the box V as an Armadillo matrix, mat.</returns>  
mat Crank::make_potential_double_slit(){
  mat V = make_potential_box(); //Creates a box of size M_* M_
  int center_index = (M_)*0.5; //200 * 0.5 = 100
  int x_thickness = 0.02/h_;// indices i in x direction: (0.02/0.005) + 1 = 5
  int x_start = center_index - x_thickness/2;
  int x_end = center_index + x_thickness/2;
  int aperture = (0.05/h_) + 1 ;// (0.05/0.005) + 1 = 11
  int center_wall_length = (0.05/h_);// (0.05/0.005) = 10. From j=95 to j=105
  int start = center_index - (center_wall_length)/2 - aperture; // j=84 for lower aperture and j=106 for upper aperture
  int end = start + aperture; // j=94 for lower aperture and j=116 for upper aperture

  for (int i = x_start; i < x_end+1; i++) {
    V.row(i).fill(v_0_);
    for(int j = start; j < end+1; j++) {
      V(i,j) = 0;
      V(i,j+center_wall_length+1+aperture) = 0;
    }
  }
  return V;
}

/// <summary>Makes the box with a potential barrier V divided in two by a barrier with triple slit openings. Parameters fetched from class variables</summary>
/// <returns>A matrix describing the box V as an Armadillo matrix, mat.</returns>  
mat Crank::make_potential_triple_slit(){
  mat V = make_potential_box(); //Creates a box of size M_* M_
  int center_index = (M_)*0.5; //200 * 0.5 = 100
  int x_thickness = 0.02/h_;// indices i in x direction: (0.02/0.005) + 1 = 5
  int x_start = center_index - x_thickness/2;
  int x_end = center_index + x_thickness/2;
  int aperture = (0.05/h_) +1 ;// (0.05/0.005) + 1 = 11
  int wall_length = (0.05/h_) + 1;// (0.05/0.005) + 1 = 11
  int start = center_index - wall_length - aperture - ((aperture-1)/2);
  int end = start + aperture;
  int unit_separation = wall_length + aperture;

  for (int i = x_start; i < x_end+1; i++) {
    V.row(i).fill(v_0_);
    for(int j = start; j < end+1; j++) {
      V(i,j) = 0;
      V(i,j + unit_separation + 1) = 0;
      V(i,j + (unit_separation * 2) + 1) = 0;
    }
  }
  return V;
}
