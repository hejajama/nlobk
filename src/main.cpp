/*
 * nloBK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2013
 */


#include "ic.hpp"
#include "mv.hpp"
#include "dipole.hpp"
#include "solver.hpp"
#include <csignal>
#include <iostream>

using std::cout;
using std::cerr;
using std::endl;
std::string version = "0.01-dev";

// We need global variables so that the signal handler works
std::string output="output.dat";


void SaveData();
void SigIntHandler(int param);

int main(int argc, char* argv[])
{

    //gsl_set_error_handler(&ErrHandler);
    std::signal(SIGINT, SigIntHandler);

    MV ic;
    cout << "Initial condition is " << ic.GetString() << endl;

    Dipole dipole(&ic);

    BKSolver solver(&dipole);
    solver.Solve(1);

    return 0;
}

void SaveData()
{
    cout << "Saving data in file " << output << endl;

}

// User pressed ctrl+c or killed the program, save data
void SigIntHandler(int param)
{
    cerr << endl << "# Received SigInt signal, trying to save data..." << endl;


    exit(1);
}
