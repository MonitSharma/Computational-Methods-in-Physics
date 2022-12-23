#include <string>
#include <iostream>
#include <armadillo>
#include <string> //includes string library
#include <math.h>
#include <iomanip>
#include <fstream>
#include <sstream>

#define pi 3.14159265359

using namespace std;

arma::mat create_tridiagonal(const arma::vec& a, const arma::vec& d, const arma::vec& e);
arma::mat create_tridiagonal(int n, double a, double d, double e);
arma::mat create_symmetric_tridiagonal(int n, double a, double d);

arma::vec analytical_eigenvalues(arma::mat A);
arma::mat analytical_eigenvectors(arma::mat A);

double find_max_value(arma::mat A, int& k, int& l);
void problem_4b();

void jacobi_rotate(arma::mat& A, arma::mat& R, int k, int l); // fra code snippets (why ref A?)
void jacobi_eigensolver(arma::mat& A, double& eps, arma::vec& eigenvalues, arma::mat& eigenvectors,
                        const int& maxiter, int& iterations, bool& converged);
void jacobi_scaling(arma::mat& A, int& N, double& eps, arma::vec& eigenvalues, arma::mat& eigenvectors,
                        int& maxiter, int& iterations, bool& converged);

void file_for_plot(int n);

int main(int argc, char const *argv[]) {

  // ----- Problem 3 -----
  int N = 6;         //size of matrix NxN
  int n = N+1;       //steps in matrix
  double h = 1./n;
  double a = -1./(h*h);     //super and sub diagonal elements
  double d = 2./(h*h);      //diagonal elements
  arma::mat A = create_symmetric_tridiagonal(N, a, d);

  // Finding eigenvalues and eigenvectors using Armadillo's arma::eig_sym:
  arma::vec eigval;
  arma::mat eigvec;
  eig_sym(eigval, eigvec, A);

  cout << endl << "------------Solution to task 3------------"<< endl;
  cout << "Eigenvalues with eig_sym:\n" << eigval << endl; // printing out
  cout << "Eigenvector with eig_sym:\n" << eigvec << endl;

  // Finding the analytical exact eigenvalues and eigenvectors
  arma::vec eigval_analytic = analytical_eigenvalues(A);
  arma::mat eigvec_analytic = analytical_eigenvectors(A);

  cout << "Analytical eigenvalues:\n" << eigval_analytic << endl;
  cout << "Analytical eigenvectors:\n" << eigvec_analytic << endl;

  cout << endl << "------------Solution to task 3(end)------------"<< endl;
  // -------- Problem 3 (end) --------

  cout << endl << "------------Solution to task 4 b)------------"<< endl;
  // Problem 4B
  problem_4b();

  cout << endl << "------------Solution to task 4 b)(end)------------"<< endl;

  // Problem 5B using the matrix A from problem 3
  double eps = 1.0e-8; // tolerance
  arma::vec eigenvalues;
  arma::mat eigenvectors;
  int maxiter = N * N * N;
  int iterations;
  bool converged = 0;
  jacobi_eigensolver(A, eps, eigenvalues, eigenvectors, maxiter, iterations, converged);


  cout << endl << "------------Solution to task 5 b)------------"<< endl;

  cout << "Eigenvalues with jacobi:\n" << eigenvalues << endl; // printing out
  cout << "Eigenvector with jacobi:\n" << eigenvectors << endl;


  cout << endl << "------------Solution to task 5 b) (end)------------"<< endl;

  // Problem 6 - comment out to avoid taking up too much time
  cout << endl << "------------Running task 6------------>"<< endl;
  jacobi_scaling(A, N, eps, eigenvalues, eigenvectors, maxiter, iterations, converged);


  // Problem 7
  cout << endl << "------------Running task 7------------>"<< endl;
  file_for_plot(10);
  file_for_plot(100);

  return 0;
}


//-----------Problem 3-------------

// ----------Creating matrix -------

// Create tridiagonal matrix from vectors.
// - lower diagonal: vector a, lenght N-1
// - main diagonal:  vector d, lenght N
// - upper diagonal: vector e, lenght N-1
arma::mat create_tridiagonal(const arma::vec& a, const arma::vec& d, const arma::vec& e)
{
  int N = d.size();
  // Start from identity matrix
  arma::mat A = arma::mat(N, N, arma::fill::eye); // A(row, column)

  // Fill first row (row index 0)
  A(0,0) = d(0);
  A(0,1) = e(0);

  // Loop that fills rows 2 to N-1 (row indices 1 to N-2)
  for (int r = 1; r <= N-2; r++){
    A(r, r-1) = a(r-1);
    A(r, r) = d(r); // diagonal element
    A(r, r+1) = e(r);
  }

  // Fill last row (row index N-1)
  A(N-1, N-2) = a(N-2);
  A(N-1, N-1) = d(N-1);

  return A;
}


// Create a tridiagonal matrix tridiag(a,d,e) of size N*N
// from scalar input a, d and e
arma::mat create_tridiagonal(int N, double a, double d, double e)
{
  arma::vec a_vec = arma::vec(N-1, arma::fill::ones) * a;
  arma::vec d_vec = arma::vec(N, arma::fill::ones) * d;
  arma::vec e_vec = arma::vec(N-1, arma::fill::ones) * e;

  // Call the vector version of this function and return the result
  return create_tridiagonal(a_vec, d_vec, e_vec);
}

// Create a symmetric tridiagonal matrix tridiag(a,d,a) of size N*N
// from scalar input a and d.
arma::mat create_symmetric_tridiagonal(int N, double a, double d)
{
  // Call create_tridiagonal and return the result
  return create_tridiagonal(N, a, d, a);
}
// -------Creating matrix (end)------------

arma::mat analytical_eigenvectors(arma::mat A){
  int N = arma::size(A,0);
  double d = A(0,0);
  double a = A(0,1);

  arma::mat v(N,N);

  for (int i = 1; i <= N; i++){
    for (int j = 1; j <= N; j++){
      v(j-1,i-1) = sin((i*j*pi)/(N+1)); // j = row, i = column (in accordance to definition)
    }
  }
  return arma::normalise(v); // return scaled eigenvectors

}
arma::vec analytical_eigenvalues(arma::mat A){
  // A is tridiagonal (a,d,a)
  int N = arma::size(A,0);
  double d = A(0,0);
  double a = A(0,1);

  arma::vec lambda(N);

  for (int i = 1; i <= N; i++){
    lambda(i-1) = d + 2*a*cos((i*pi)/(N+1));
  }
  return lambda;
}
//-----------Problem 3 (end)-------------

// ----------Problem 4A -------------
double find_max_value(arma::mat A, int& k, int& l){

  double max_value = 0;
  int N = arma::size(A,0); //i is dimension N of matrix A

  //for loop runs through the non-diagonal matrix elements under the diagonal.
  for (int j=0; j<=N-1; j++){

          for (int i=1+j; i<=N-1; i++){
            if(abs(A(i,j))>abs(max_value)){
              max_value= A(i,j);
              k = j, l=i; //column k and row l
            }
      }
  }
  return max_value;
}
// ---------------Problem 4A (end)----------


//---------------Problem 4B-------------
void problem_4b(){

  //defines and fills the matix shown in Problem 4b)
  arma::mat B_4 = arma::mat(4, 4).fill(0.);
  for (int i = 0; i < 3; i++){  // row
      for (int j = 0; j < 3; j++){ // // column
        if (i == j){
          B_4(i,j) = 1;
        }
      }
    }
    B_4(0,3) = 0.5; B_4(1,2) = -0.7; B_4(2,1) = -0.7; B_4(3,0) = 0.5;


  //prints the matrix to terminal and calls the function from Problem 3
  //to find en print the max value of non-diagonal matrix element in the
  //sub triangular matrix.

  int k; int l;
  cout<< B_4 <<endl;
  //returns max value and assigns k as the column index and l as the row index
  cout <<"max value: "<< find_max_value(B_4,k,l) <<" row: "<<l<<" column: "<< k << endl;
}
//---------------Problem 4B(end)-------------

//---------------Problem 5A------------------

void jacobi_rotate(arma::mat& A, arma::mat& R, int k, int l){ // SJEKK INDEXER A(row, column)
    //Computing tan (t), cos (c) and sin (s)
    int N = arma::size(A,0);
    double s, c;
    if ( A(k,l) != 0.0){
        double t, tau;
        tau = (A(l,l)-A(k,k))/(2*A(k,l));
        if ( tau > 0){
            t = 1.0/(tau + sqrt(1.0 + tau*tau));
        }
        else {
            t = -1.0/( -tau + sqrt(1.0 + tau*tau));
        }
        c = 1.0/(sqrt(1+t*t));
        s = c*t;
    }
    else {
        c = 1.0;
        s = 0.0;
    }

    //Transform current A matrix
    double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
    a_kk = A(k,k);
    a_ll = A(l,l);
    A(k,k) = c*c*a_kk - 2.0*c*s*A(k,l) + s*s*a_ll;
    A(l,l) = s*s*a_kk + 2.0*c*s*A(k,l) + c*c*a_ll;
    A(k,l) = 0.0;
    A(l,k) = 0.0;
    for ( int i = 0; i < N; i++ ) {
        if ( i != k && i != l ) {
            a_ik = A(i,k);
            a_il = A(i,l);
            A(i,k) = c*a_ik - s*a_il;
            A(k,i) = A(i,k);
            A(i,l) = c*a_il + s*a_ik;
            A(l,i) = A(i,l);
        }
        //Compute new eigenvectors
        r_ik = R(i,k);
        r_il = R(i,l);
        R(i,k) = c*r_ik - s*r_il;
        R(i,l) = c*r_il + s*r_ik;
    }
    return;
}

// Jacobi method eigensolver:
// - Runs jacobi_rotate until max off-diagonal element < eps
// - Writes the eigenvalues as entries in the vector "eigenvalues"
// - Writes the eigenvectors as columns in the matrix "eigenvectors"
//   (The returned eigenvalues and eigenvectors are sorted using arma::sort_index)
// - Stops if it the number of iterations reaches "maxiter"
// - Writes the number of iterations to the integer "iterations"
// - Sets the bool reference "converged" to true if convergence was reached before hitting maxiter

void jacobi_eigensolver(arma::mat& A, double& eps, arma::vec& eigenvalues, arma::mat& eigenvectors,
                        const int& maxiter, int& iterations, bool& converged)
                        // chose to not define A as a constant to let this method change it
{
  int N = arma::size(A,0);
  iterations = 0;

  int k;
  int l;
  double max_value = find_max_value(A, k, l); //( A, &k, &l);
  arma::mat R = arma::mat(N, N, arma::fill::eye);

  while (fabs(max_value) > eps && (double) iterations < maxiter ) {
      max_value = find_max_value(A, k, l); //(A, &k, &l);
      jacobi_rotate(A, R, k, l);
      iterations++;
  }

  arma::vec diagonals = A.diag(); // eigenvalues are diagonal elements of rotated matrix A
  arma::uvec indices = sort_index(diagonals, "ascending");

  // Sorting and filling eigenvalues and eigenvectors
  eigenvalues = arma::vec(N);
  eigenvectors = arma::mat(N,N);
  for (int i = 0; i < N; i++){
    eigenvalues(i) = diagonals(indices(i));
    eigenvectors.col(i) = R.col(indices(i));
  }

//converged set to 1 means the jacobi rotation converged
 if(iterations+1 != maxiter){
  converged = 1;
 }
}
// ---------Problem 5A(end)-----------


//-----------Problem 6-------------


void jacobi_scaling(arma::mat& A, int& N, double& eps, arma::vec& eigenvalues, arma::mat& eigenvectors,
                        int& maxiter, int& iterations, bool& converged){

  ofstream myfile;

  myfile.open("./datafiles/task6_dataset.txt");
  myfile << "N=    , Rotations=    "<<endl;

  for (int N = 3; N < 110; N++){
    int n = N+1;       //steps in matrix
    double h = 1./n;
    double a = -1./(h*h);     //super and sub diagonal elements
    double d = 2./(h*h);      //diagonal elements


    A = create_symmetric_tridiagonal(N,a,d); //creates an NxN tridaiag symmetric matrix
    maxiter = (int) N * (int) N * (int) N;   //redefines max number of itterations
    jacobi_eigensolver(A, eps, eigenvalues, eigenvectors, maxiter, iterations, converged); //runs jacobi eigensolver
    myfile <<N<<","<< iterations << endl;   //writes N and number of itterations to a txt file
  }

  myfile.close();
}

//-----------Problem 6(end)-------------

// -----------Problem 7------------

void file_for_plot(int n){
  // n = steps

  // create symmetric tridiagonal matrix with:
  int N = n-1;         //size of matrix NxN
  double h = 1./n;
  double a = -1./(h*h);     //super and sub diagonal elements
  double d = 2./(h*h);      //diagonal elements
  arma::mat A = create_symmetric_tridiagonal(N, a, d);

  // Finding analytical eigenvectors:
  arma::vec eigval_analytic = analytical_eigenvalues(A);
  arma::mat eigvec_analytic = analytical_eigenvectors(A);

  // Solve matrix equation Ax = b
  double eps = 1.0e-8; // tolerance
  arma::vec eigenvalues;
  arma::mat eigenvectors;
  int maxiter = N * N * N;
  int iterations;
  bool converged = 0;
  jacobi_eigensolver(A, eps, eigenvalues, eigenvectors, maxiter, iterations, converged);

  arma::vec xhat = arma::vec(n+1);
  xhat(0) = 0; // boundary point
  for (int i = 1; i < n; i++){
    xhat(i) = xhat(0) + i*h;
  }
  xhat(n) = 1; // boundary point

  arma::mat vhat = arma::mat(N, 3);
  arma::mat v_analytic = arma::mat(N, 3);

  // Make sure analytical eigenvectors are sorted in ascending order:
  arma::uvec indices = sort_index(eigval_analytic, "ascending"); // make sure

  // only interested in eigenvectors corresponding to 3 smallest eigenvalues
  for (int i = 0; i < 3; i++){
    vhat.col(i) = eigenvectors.col(i);
    v_analytic.col(i) = eigvec_analytic.col(indices(i));
  }

  // Boundary conditions v(0) = v(n) = 0
  arma::rowvec boundary = arma::rowvec(3, arma::fill::zeros);
  vhat.insert_rows(0, boundary);
  vhat.insert_rows(N+1, boundary); // N+1 = n
  v_analytic.insert_rows(0, boundary);
  v_analytic.insert_rows(N+1, boundary); // N+1 = n


  ofstream ofile;
  std::ostringstream filename;
  filename << "./datafiles/output" << n << ".txt";
  ofile.open(filename.str());
  int width = 20;
  int prec = 12;
  for (int i = 0; i <= n; i++){
    ofile << setw(width) << setprecision(prec) << scientific << xhat(i)
          << setw(width) << setprecision(prec) << scientific << vhat(i, 0)
          << setw(width) << setprecision(prec) << scientific << vhat(i, 1)
          << setw(width) << setprecision(prec) << scientific << vhat(i, 2)
          << setw(width) << setprecision(prec) << scientific << v_analytic(i, 0)
          << setw(width) << setprecision(prec) << scientific << v_analytic(i, 1)
          << setw(width) << setprecision(prec) << scientific << v_analytic(i, 2) << endl;
  }
  ofile.close();

}
//------- Problem 7(end)------------
