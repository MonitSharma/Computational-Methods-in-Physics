#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <chrono>
#include "./include/Crank.hpp"
#include <armadillo>
#include <omp.h>
#include <complex>
#include <array>

using namespace std::complex_literals; // to use imaginary number i | DEMANDS c++14!
using namespace std;
using namespace arma;

// Performs the simulations specified in the input file
void simulate(string inputfile);

// Converts to string with specified precision
template <typename T> string to_string_with_precision(const T a_value, const int n = 2) {
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return out.str();
}

int main(int argc, char const *argv[]) {
  if (argc != 2) {
    // Get the name of the executable file
    std::string executable_name = argv[0];

    std::cerr << "Error: Wrong number of input arguments." << std::endl;
    std::cerr << "Usage: " << executable_name << " <inputfile.txt>" << std::endl;

    // Exit program with non-zero return code to indicate a problem
    return 1;   
  }
  
  // Runs the simulations specified in the input file
  simulate(argv[1]); 
  return 0;
}


/// <summary>Performs the simulations specified in the input file</summary>
/// <param name="inputfile">Filename of file that specifies what simulations to run.</param>  
/// <returns>void</returns>  
void simulate(string inputfile) {
  fstream myfile;

  myfile.open(inputfile);
  if (myfile.is_open()){
    // This checks that the file was opened OK
    string line;
    std::getline(myfile, line); // skip the first line
    const size_t input_vals = 14;
    std::vector<std::array<double, input_vals>> vars;

    int line_counter = 0;
    while (std::getline(myfile, line)) {
      double prob, h, deltat, T, xc, sx, px, yc, sy, py, v0, reim, last_slice, slits;
      std::stringstream mysstream(line);
      mysstream >> prob >> h >> deltat >> T >> xc >> sx >> px >> yc >> sy >> py >> v0 >> slits >> reim >> last_slice;
      vars.push_back({{prob, h, deltat, T, xc, sx, px, yc, sy, py, v0, slits, reim, last_slice}});

      line_counter +=1;
    }
    int width = 10;
    cout << std::setw(width) << "Problem" << std::setw(width) << "h" << std::setw(width) << "deltat" << std::setw(width) << "T" << std::setw(width) << "x_c" << std::setw(width)
    << "sigma_x" << std::setw(width) << "p_x" << std::setw(width) << "y_c" << std::setw(width) << "sigma_y" << std::setw(width) << "p_y" << std::setw(width) << "v_0"
    << std::setw(width) << "slits"  << std::setw(width) << "ReIm" << std::setw(width) << "Last_slice" << endl;

    for (int i = 0; i<vars.size(); i++) {
      for (int j = 0; j < (int)input_vals; j++) {
        cout << std::setw(width) << vars[i][j] ;
      }
      cout << endl;
    }
    cout << "Number of simulations to run is " << line_counter << ". Running in parallel." << endl;

    #pragma omp parallel for
    for (int i = 0; i <line_counter; i++) {
      int thread_id = omp_get_thread_num();
      double prob, h, deltat, T, xc, sx, px, yc, sy, py, v0, reim, last_slice, slits;
      prob=vars[i][0]; h=vars[i][1]; deltat=vars[i][2]; T=vars[i][3]; xc=vars[i][4]; sx=vars[i][5]; px=vars[i][6]; yc=vars[i][7]; sy=vars[i][8];
      py=vars[i][9]; v0=vars[i][10]; slits=vars[i][11]; reim=vars[i][12]; last_slice=vars[i][13];

      #pragma omp critical
      cout << "\nRunning problem " << to_string_with_precision(prob) << " on thread " << thread_id << " with parameters: \nh=" << h << ", delta t=" << deltat << ", T=" <<  T << ", x_c=" << xc << ", s_x=" << sx << ", p_x=" << px << ", y_c=" << yc
      << ", s_y=" << sy << ", p_y=" << py << ", v_0=" << v0 << ", slits=" << (int)slits << ", ReIm=" << reim << ", last slice=" << last_slice << "." << endl;

      Crank crank(h, deltat, T, xc, yc, sx, sy, px, py, v0, slits);
      cx_cube results_cube;
      cx_mat last_slice_mat;

      if (last_slice == 1) {
        last_slice_mat = crank.run_simulation(1);
      } else {
        results_cube = crank.run_simulation();
      }

      // For problem 8-2. Results cube before being converted to probabilities. Raw Re and Im.
      if (reim == 1) {
        results_cube.save("datafiles/Problem_"+to_string_with_precision(prob)+"_ReIm_outputCube_slits_" + to_string((int)slits) + ".dat");
      }

      if (last_slice == 1) {
        last_slice_mat.save("datafiles/Problem_"+to_string_with_precision(prob)+"_outputMat_slits_" + to_string((int)slits) + ".dat");
      } else {
        results_cube.save("datafiles/Problem_"+to_string_with_precision(prob)+"_outputCube_slits_" + to_string((int)slits) + ".dat");
      }

      //Output box shapes once
      if(i == 0) {
        crank.v_0_ = 1e10;
        mat box = crank.make_potential_box();
        mat box_single_slit = crank.make_potential_single_slit();
        mat box_double_slit = crank.make_potential_double_slit();
        mat box_triple_slit = crank.make_potential_triple_slit();
        box.save("datafiles/box.dat");
        box_single_slit.save("datafiles/box_single_slit.dat");
        box_double_slit.save("datafiles/box_double_slit.dat");
        box_triple_slit.save("datafiles/box_triple_slit.dat");
      }
    }
  }
  else {
    cout << "Unable to open the file " << inputfile << endl;
  }
  myfile.close();
}
