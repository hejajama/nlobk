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
    string output = "output.dat";
    
    char *hostname = new char[500];
    gethostname(hostname, 500);
    
    cout <<"#"<<endl<<"# NLOBK solver " << version  << " running on " << hostname << ", start at " << today << "#" << endl;
    delete[] hostname;
    
    gsl_set_error_handler(&ErrHandler);
    //std::signal(SIGINT, SigIntHandler);

    InitialCondition* ic = NULL;
    for (int i=1; i<argc; i++)
    {
        if (string(argv[i])=="-ic")
        {
            if (string(argv[i+1])=="FILE")
            {
                ic = new IC_datafile();
                ((IC_datafile*)ic)->LoadFile(argv[i+2]);
            }
        }
        else if (string(argv[i])=="-output")
            output = argv[i+1];
        else if (string(argv[i]).substr(0,1)=="-") 
        {
            cerr << "Unrecoginzed parameter " << argv[i] << endl;
            return -1;
        }
    }
    if (ic == NULL)
    {
        ic = new MV();
        ((MV*)ic)->SetQsqr(0.2);
    }

    //MV ic;
	//ic.SetQsqr(0.2);

    
    cout << "# " << NLOBK_CONFIG_STRING() << endl;
    cout << "#Initial condition is " << ic->GetString() << endl;

    Dipole dipole(ic);


    cout <<"# r grid size: " << dipole.RPoints() << " minr " << dipole.MinR() << " maxr " << dipole.MaxR() << endl;

    //cout << "N(r=0.001)=" << dipole.N(0.001) <<", N(r=0.1)=" << dipole.N(0.1) <<", N(r=10)=" << dipole.N(10) << endl;

    BKSolver solver(&dipole);
    solver.SetTmpOutput(output);
    solver.Solve(2);
    cout << "BK solved!" << endl;

	

    cout << "Saving to file " << output << endl;
    dipole.Save(output);

    
    time_t end = time(0);
    int diff = end-start;
    cout << "Solution took " << diff/60.0/60.0 << " hours" << endl;

    delete ic;
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
    if (INTMETHOD_NLO == MISER)
        ss <<"MonteCarlo Miser, points=" << MCINTPOINTS;
    else if (INTMETHOD_NLO == VEGAS)
        ss <<"MonteCarlo Vegas, points=" << MCINTPOINTS;
    else if (INTMETHOD_NLO == MULTIPLE)
        ss << "Multiple integrals (no montecarlo)";
    else
        ss <<"UNKNOWN!";
    ss<< ". LO Kernel RC: ";
    if (RC_LO == FIXED_LO or EQUATION==CONFORMAL_N4)
        ss << " fixed as=" << FIXED_AS;
    else if (RC_LO == SMALLEST_LO)
        ss << " smallest dipole";
    else if (RC_LO == BALITSKY_LO)
        ss << " Balitsky";
    else
        ss << " NO STRING IMPLEMENTED!";

    ss<< ". NLO Kernel RC: ";
    if (RC_NLO == FIXED_NLO or EQUATION==CONFORMAL_N4)
        ss << " fixed as=" << FIXED_AS;
    else if (RC_NLO == SMALLEST_NLO)
        ss << " smallest dipole";
    else
        ss << " NO STRING IMPLEMENTED!";

    ss <<". Nc=" << NC << ", Nf=" << NF;

    if (EQUATION == QCD)
    {
        if (DOUBLELOG_LO_KERNEL) ss << ". QCD, Double log term in LO kernel included";
        else ss << ". QCD, Double log term in LO kernel NOT included";
    }
    else if (EQUATION == CONFORMAL_QCD) ss << ". Solving for CONFORMAL dipole";
    else if (EQUATION == CONFORMAL_N4) ss << ". Solving in N=4 for CONFORMAL dipole";
    else ss << ". UNKNOWN EQUATION!!";

    if (LO_BK)
        ss <<" Solving Leading Order BK";

    return ss.str();
}
