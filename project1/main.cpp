#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <armadillo>
#include <time.h>
using namespace std;

double f(double x);
arma::vec general_algorithm(arma::vec a, arma::vec b, arma::vec c, arma::vec g, int n);
arma::vec special_algorithm(arma::vec g, int n);


int main(int argc, const char * argv[]) {

    if (argc != 2)
    {
      // Get the name of the executable file
      std::string executable_name = argv[0];

      std::cerr << "Error: Wrong number of input arguments." << std::endl;
      std::cerr << "Usage: " << executable_name << " <integer number of steps>" << std::endl;
      return 1;
    }

    int N = atoi(argv[1]); // number of steps

    int m = N+1; // size of exact solution
    arma::vec u = arma::vec(m);
    arma::vec x = arma::vec(m);
    double x_min = 0.0;
    double x_max = 1.0;
    double h = (x_max - x_min) / N;
    // Boundary conditions:
    double u_0 = 0.;  // u(0) = 0
    double u_1 = 0.;  // u(1) = 0

    int width = 30; //12
    int prec = 10; // 4

    //opening file
    ofstream ofile;
    std::ostringstream filename;
    filename << "./datafiles/exact_data" << N << ".txt";
    ofile.open(filename.str());

    //setting up the x-array and the solutions to the function u, and printing it to file
    for (int i=0 ; i <= m-1 ; i++){
        x(i) = h*i;
        u(i) = 1 - (1 - exp(-10)) * x(i) - exp(-10 * x(i));
        ofile << setw(width) << setprecision(prec) << scientific << x(i)
              << setw(width) << setprecision(prec) << scientific << u(i) << endl;
    }
    ofile.close(); //close file


    // Problem 7
    int n = N-1; // number of unknowns

    arma::vec a = arma::vec(n).fill(-1.);
    arma::vec b = arma::vec(n).fill(2.);
    arma::vec c = arma::vec(n).fill(-1.);

    // Making sure a and c only have n-1 values, but indices corresponding to row in A
    a(0) = 0.;
    c(n-1) = 0.;


    // Defining the g vector:
    arma::vec g = arma::vec(n);
    /*
    x is a vector of length m = n + 2 including the boundary elements, while g
    has length n where the boundary values are "ignored".
    Therefore, element g(i) will correspond to element x(i+1) for i = 0,...,n-1
    */
    g(0) = h*h*f(x(1)) + u_0;
    g(n-1) = h*h*f(x(n)) + u_1;

    for (int i = 1; i <= n-2; i++){
      g(i) = h*h*f(x(i+1));
    }

    // General algorithm:
    arma::vec v = general_algorithm(a,b,c,g,n);

    //opening file
    ofstream ofile2;
    std::ostringstream filename2;
    filename2 << "./datafiles/approx_general" << N << ".txt";
    ofile2.open(filename2.str());

    //setting up the x-array and the solutions to the function u, and printing it to file
    for (int i=0 ; i <= n-1 ; i++){
        ofile2 << setw(width) << setprecision(prec) << scientific << x(i+1)
              << setw(width) << setprecision(prec) << scientific << v(i) << endl;
    }
    ofile2.close();    //close file
     double largestvalue = 0;
    // Problem 8:
    ofstream ofile3;
    std::ostringstream filename3;
    filename3 << "./datafiles/errors" << N << ".txt";
    ofile3.open(filename3.str());
    for (int i=0 ; i <= n-1 ; i++){
      double abs_err = std::abs(u(i+1)-v(i));
      double rel_err = std::abs(abs_err / u(i+1));
      ofile3 << setw(width) << setprecision(prec) << scientific << x(i+1)//log10(x(i+1))
            << setw(width) << setprecision(prec) << scientific << log10(abs_err)
            << setw(width) << setprecision(prec) << scientific << log10(rel_err) << endl;
        if(largestvalue < rel_err ){
            largestvalue = rel_err;
        }
    }
    ofile3.close();

    cout << "The max relative error for N = " << N<< " is: " << setprecision(prec) << scientific << largestvalue << endl;


    // Problem 9:
    // A is tridiagonal matrix. Solve Av^ = g where v^ denotes v using special algorithm.
    arma::vec vhat = special_algorithm(g,n);

    ofstream ofile4;
    std::ostringstream filename4;
    filename4 << "./datafiles/approx_special" << N << ".txt";
    ofile3.open(filename4.str());

    //setting up the x-array and the solutions to the function u, and printing it to file
    for (int i=0 ; i <= n-1 ; i++){
        ofile4 << setw(width) << setprecision(prec) << scientific << x(i+1)
              << setw(width) << setprecision(prec) << scientific << vhat(i) << endl;
    }
    ofile4.close(); //close file

    cout << "\n" << endl;
    return 0;
}

arma::vec general_algorithm(arma::vec a, arma::vec b, arma::vec c, arma::vec g, int n){
    clock_t t1 = clock();
    // Helpful new variables
    arma::vec btilde = arma::vec(n);
    arma::vec gtilde = arma::vec(n);
    btilde(0) = b(0);
    gtilde(0) = g(0);

    arma::vec v = arma::vec(n); // solution vector
    double tmp; // variable to reduce FLOPs

    for (int i = 1; i <= n - 1; i++){
      tmp = a(i) / btilde(i-1);
      btilde(i) = b(i) - tmp * c(i-1) ;
      gtilde(i) = g(i) - tmp * gtilde(i-1);
    }

    v(n-1) = gtilde(n-1) / btilde(n-1); // Last element can now be found directly

    for (int j = n-2; j >= 0; j--){ // n-1 elements
      v(j) = (gtilde(j) - (c(j) * v(j+1))) / btilde(j);
    }
    clock_t t2 = clock();
    double duration_seconds = ((double) (t2-t1))/CLOCKS_PER_SEC;
    cout << "Timing for general algorithm: " << duration_seconds<< '\n';
    return v;

}

arma::vec special_algorithm(arma::vec g, int n){
    clock_t t3 = clock();
    // Helpful new variables
    arma::vec btilde = arma::vec(n);
    arma::vec gtilde = arma::vec(n);
    btilde(0) = 2.;
    gtilde(0) = g(0);

    arma::vec v = arma::vec(n); // solution vector
    double tmp; // variable to reduce FLOPs

    for (int i = 1; i <= n - 1; i++){
      tmp = 1 / btilde(i-1);
      btilde(i) = 2 - tmp;
      gtilde(i) = g(i) + tmp * gtilde(i-1);
    }

    v(n-1) = gtilde(n-1) / btilde(n-1); // Last element can now be found directly

    for (int j = n-2; j >= 0; j--){ // n-1 elements
      v(j) = (gtilde(j) +v(j+1)) / btilde(j);
    }
    clock_t t4 = clock();
    double duration_seconds2 = ((double) (t4-t3))/CLOCKS_PER_SEC;
    cout <<  "Timing for special algorithm: " << duration_seconds2 << '\n';
    return v;

}

double f(double x){
  return 100*exp(-10*x);
}
