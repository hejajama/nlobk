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
#include <iomanip>
#include <gsl/gsl_errno.h>

using std::cout;
using std::cerr;
using std::endl;
std::string version = "0.01-dev";

// We need global variables so that the signal handler works
std::string output="output.dat";


void SaveData();
void SigIntHandler(int param);
void ErrHandler(const char * reason,
                        const char * file,
                        int line,
                        int gsl_errno);

int main(int argc, char* argv[])
{

    gsl_set_error_handler(&ErrHandler);
    //std::signal(SIGINT, SigIntHandler);

    MV ic;
	ic.SetQsqr(0.2);
    cout << "#Initial condition is " << ic.GetString() << endl;

    Dipole dipole(&ic);


    cout <<"# r grid size: " << dipole.RPoints() << " minr " << dipole.MinR() << " maxr " << dipole.MaxR() << endl;

    //cout << "N(r=0.001)=" << dipole.N(0.001) <<", N(r=0.1)=" << dipole.N(0.1) <<", N(r=10)=" << dipole.N(10) << endl;

    BKSolver solver(&dipole);
    solver.Solve(2);
    cout << "BK solved!" << endl;

	std::string output=std::string(argv[1]);

    cout << "Saving to file " << output << endl;
    dipole.Save(output);

    return 0;
}



// User pressed ctrl+c or killed the program, save data
void SigIntHandler(int param)
{
    //cerr << endl << "# Received SigInt signal, trying to save data..." << endl;


    exit(1);
}



int errors;
void ErrHandler(const char * reason,
                        const char * file,
                        int line,
                        int gsl_errno)
{
    
    // Errors related to convergence of integrals are handled when
    // gsl_integration functions are called, don't do anything with them here
    // 14 = failed to reach tolerance
    // 18 = roundoff error prevents tolerance from being achieved
    // 11 = maximum number of subdivisions reached
    if (gsl_errno == 14 or gsl_errno == 18 or gsl_errno == 11)
        return;

    // 15: underflows
    if (gsl_errno == 15 ) return;

    errors++;
    std::cerr << file << ":"<< line <<": Error " << errors << ": " <<reason
            << " (code " << gsl_errno << ")." << std::endl;
}
