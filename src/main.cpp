/*
 * nloBK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2013
 */


#include "ic.hpp"
#include "mv.hpp"
#include "ic_datafile.hpp"
#include "dipole.hpp"
#include "solver.hpp"
#include <csignal>
#include <iostream>
#include <iomanip>
#include <gsl/gsl_errno.h>
#include <ctime>
#include <unistd.h>

using namespace std;

string version = "0.01-dev";


std::string NLOBK_CONFIG_STRING();
void SaveData();
void SigIntHandler(int param);
void ErrHandler(const char * reason,
                        const char * file,
                        int line,
                        int gsl_errno);

int main(int argc, char* argv[])
{
    time_t start = time(0);
    string today = ctime(&start);
    
    char *hostname = new char[500];
    gethostname(hostname, 500);
    
    cout <<"#"<<endl<<"# NLOBK solver " << version  << " running on " << hostname << ", start at " << today << "#" << endl;
    delete[] hostname;
    
    gsl_set_error_handler(&ErrHandler);
    //std::signal(SIGINT, SigIntHandler);

    MV ic;
	ic.SetQsqr(0.2);
    //IC_datafile ic;
    //ic.LoadFile("../amplitudelib_v2/amplitudelib2/nlobk_analyysit/loglog/nlo_smallest_y_2");
    
    cout << "# " << NLOBK_CONFIG_STRING() << endl;
    cout << "#Initial condition is " << ic.GetString() << endl;

    Dipole dipole(&ic);


    cout <<"# r grid size: " << dipole.RPoints() << " minr " << dipole.MinR() << " maxr " << dipole.MaxR() << endl;

    //cout << "N(r=0.001)=" << dipole.N(0.001) <<", N(r=0.1)=" << dipole.N(0.1) <<", N(r=10)=" << dipole.N(10) << endl;

    BKSolver solver(&dipole);
    solver.Solve(10);
    cout << "BK solved!" << endl;

	std::string output=std::string(argv[1]);

    cout << "Saving to file " << output << endl;
    dipole.Save(output);

    
    time_t end = time(0);
    int diff = end-start;
    cout << "Solution took " << diff/60.0/60.0 << " hours" << endl;

    return 0;
}



// User pressed ctrl+c or killed the program, save data
void SigIntHandler(int param)
{
    //cerr << endl << "# Received SigInt signal, trying to save data..." << endl;


    exit(1);
}

std::string NLOBK_CONFIG_STRING()
{
    std::stringstream ss;
    ss << "Integration method: ";
    if (MONTECARLO)
        ss <<"MonteCarlo, points=" << MCINTPOINTS;
    else
        ss <<"Direct integration";
    ss<< ". LO Kernel RC: ";
    if (RC_LO == FIXED_LO)
        ss << " fixed as=" << FIXED_AS;
    else if (RC_LO == SMALLEST_LO)
        ss << " smallest dipole";
    else if (RC_LO == BALITSKY_LO)
        ss << " Balitsky";
    else
        ss << " NO STRING IMPLEMENTED!";

    ss<< ". NLO Kernel RC: ";
    if (RC_NLO == FIXED_NLO)
        ss << " fixed as=" << FIXED_AS;
    else if (RC_NLO == SMALLEST_NLO)
        ss << " smallest dipole";
    else
        ss << " NO STRING IMPLEMENTED!";

    ss <<". Nc=" << NC << ", Nf=" << NF;

    if (DOUBLELOG_LO_KERNEL) ss << ". Double log term in LO kernel included";
    else ss << ". Double log term in LO kernel NOT included";

    if (CONFORMAL_DIPOLE) ss << ". Solving for CONFORMAL dipole";
    else ss << ". Solving for standard NON-CONFORMAL dipole";

    return ss.str();
}
