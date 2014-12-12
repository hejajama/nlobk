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
#include <tools/tools.hpp>

using namespace std;

string version = "0.01-dev";

using namespace config;


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
    double maxy=10;
    
    char *hostname = new char[500];
    gethostname(hostname, 500);
    
    cout <<"#"<<endl<<"# NLOBK solver " << version  << " running on " << hostname << ", start at " << today << "#" << endl;
    delete[] hostname;

    if (string(argv[1])=="-help")
    {
        cout <<"-ic FILE filename: set initial condition (default: MV)" << endl;
        cout << "-maxy yval" << endl;
        cout <<"-output filename_to_save_data" << endl;
        cout << "-eq qcd,confqcd,n4: set equation to solve" << endl;
        cout << "-rc smallest,parent,parent_beta,fixed: set running coupling" << endl;
        cout << "-lo: solve LO BK" << endl;
        cout << "-only_nlo: keep only nlo terms" << endl;
        cout << "-nf nf: set number of quark flavors" << endl;
        cout << "-nolimit: do not force N>=0" << endl;
        cout << "-ic [mv,mve,mvgamma] set initial condition" << endl;
        cout << "-dndy: print dn/dy at initial condition and exit" << endl;
        cout << "-alphas_scaling C^2: set C^2 [setting mv/mve/mvgamma ic sets this also]" << endl;

        return 0;
    }
    
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
            else
            {
                ic = new MV();
                if (string(argv[i+1])=="mv")
                {
                    ((MV*)ic)->SetQsqr(0.104);
                    config::ALPHAS_SCALING = 14.5;
                }
                else if (string(argv[i+1])=="mv02")
                {
                    ((MV*)ic)->SetQsqr(0.2);
                    config::ALPHAS_SCALING = 1;
                }
                else if (string(argv[i+1])=="mv1")
                {
                    ((MV*)ic)->SetQsqr(1);
                    config::ALPHAS_SCALING = std::exp(-2.0*0.57721);
                }
                else if (string(argv[i+1])=="mv10")
                {
                    ((MV*)ic)->SetQsqr(10);
                    config::ALPHAS_SCALING = std::exp(-2.0*0.57721);
                }
                else if (string(argv[i+1])=="mve")
                {
                    ((MV*)ic)->SetQsqr(0.06);
                    ((MV*)ic)->SetE(18.9);
                    config::ALPHAS_SCALING = 7.2;
                }
                else if (string(argv[i+1])=="mvgamma")
                {
                    ((MV*)ic)->SetQsqr(0.165);
                    ((MV*)ic)->SetAnomalousDimension(1.135);
                    config::ALPHAS_SCALING = 6.35;
                }
                else if (string(argv[i+1])=="mvgamma_08")
                {
                    ((MV*)ic)->SetQsqr(1);
                    ((MV*)ic)->SetAnomalousDimension(0.8);
                    config::ALPHAS_SCALING = std::exp(-2.0*0.57721);
                }
                else if (string(argv[i+1])=="mvgamma_08_qsqr_100")
                {   
                    ((MV*)ic)->SetQsqr(100);
                    ((MV*)ic)->SetAnomalousDimension(0.8);
                    config::ALPHAS_SCALING = std::exp(-2.0*0.57721);
                } 
                else if (string(argv[i+1])=="mvgamma_09")
                {
                    ((MV*)ic)->SetQsqr(1);
                    ((MV*)ic)->SetAnomalousDimension(0.9);
                    config::ALPHAS_SCALING = std::exp(-2.0*0.57721);
                }
                else if (string(argv[i+1])=="mvgamma_095")
                {
                    ((MV*)ic)->SetQsqr(1);
                    ((MV*)ic)->SetAnomalousDimension(0.95);
                    config::ALPHAS_SCALING = std::exp(-2.0*0.57721);
                }
                else if (string(argv[i+1])=="mvgamma_1")
                {
                    ((MV*)ic)->SetQsqr(1);
                    ((MV*)ic)->SetAnomalousDimension(1.0);
                    config::ALPHAS_SCALING = std::exp(-2.0*0.57721);
                }
                else
                {
                    cerr << "Uknown initial condition " << argv[i+1] << endl;
                    return -1;
                }
            }
        }
        else if (string(argv[i])=="-output")
            output = argv[i+1];
        else if (string(argv[i])=="-maxy")
            maxy = StrToReal(argv[i+1]);
        else if (string(argv[i])=="-eq")
        {
            if (string(argv[i+1])=="qcd")
                config::EQUATION = config::QCD;
            else if (string(argv[i+1])=="confqcd")
                config::EQUATION = config::CONFORMAL_QCD;
            else if (string(argv[i+1])=="n4")
                config::EQUATION = config::CONFORMAL_N4;
            else
            {
                cerr << "Unknown equation " << argv[i+1] << "! " << LINEINFO;
                return -1;
            }
        }
        else if (string(argv[i])=="-rc")
        {
            if (string(argv[i+1])=="parent")
            {
                config::RC_LO = config::PARENT_LO;
                config::RC_NLO = config::PARENT_NLO;
            }
            else if (string(argv[i+1])=="parent_beta")
            {
                config::RC_LO = config::PARENT_BETA_LO;
                config::RC_NLO = config::PARENT_NLO;
            }
            else if (string(argv[i+1])=="smallest")
            {
                config::RC_LO = config::SMALLEST_LO;
                config::RC_NLO = config::SMALLEST_NLO;
            }
            else if (string(argv[i+1])=="fixed")
            {
                config::RC_LO = config::FIXED_LO;
                config::RC_NLO = config::FIXED_NLO;
            }
            else if (string(argv[i+1])=="balitsky")
            {
                config::RC_LO = config::BALITSKY_LO;
                config::RC_NLO = config::PARENT_NLO;
            }
            else
            {
                cerr << "Unknown RC " << argv[i+1] << " at " << LINEINFO << endl;
                return -1;
            }
        }
        else if (string(argv[i])=="-lo")
            config::LO_BK = true;
        else if (string(argv[i])=="-only_nlo")
            config::ONLY_NLO = true;
        else if (string(argv[i])=="-nlo")
            config::LO_BK = false;
        else if (string(argv[i])=="-nf")
        {
            config::NF = StrToInt(argv[i+1]);
            if (config::NF != 0 and config::NF != 3)
            {
                cerr << "Invalid Nf=" << config::NF << " " << LINEINFO << endl;
                return -1;
            }
        }

            
        
        else if (string(argv[i])=="-nolimit")
            config::FORCE_POSITIVE_N = false;

        else if (string(argv[i])=="-dndy")
            config::DNDY = true;

        else if (string(argv[i])=="-alphas_scaling")
            config::ALPHAS_SCALING = StrToReal(argv[i+1]);
    
        else if (string(argv[i]).substr(0,1)=="-") 
        {
            cerr << "Unrecoginzed parameter " << argv[i] << endl;
            return -1;
        }

    }


    /* Check that configs are not contraversal
     */

    if (ic == NULL)
    {
        cerr << "Initial condition was not set!" << endl;
        return -1;
    }

    if (config::ONLY_NLO == true and config::LO_BK == true)
    {
        cerr << "Asked to solve LO bk with only NLO terms!" << endl;
        exit(1);
    }
    
    cout << "# " << NLOBK_CONFIG_STRING() << endl;
    cout << "#Initial condition is " << ic->GetString() << endl;

    Dipole dipole(ic);


    cout <<"# r grid size: " << dipole.RPoints() << " minr " << dipole.MinR() << " maxr " << dipole.MaxR() << endl;

    //cout << "N(r=0.001)=" << dipole.N(0.001) <<", N(r=0.1)=" << dipole.N(0.1) <<", N(r=10)=" << dipole.N(10) << endl;

    if (config::DNDY)
        cout << "# r   dN/dy [K1]   dN/dy [K2]  N" << endl;

    BKSolver solver(&dipole);
    solver.SetTmpOutput(output);
    solver.Solve(maxy);
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
        else if (RC_LO == PARENT_LO)
            ss << " Parent dipole";
        else if (RC_LO == PARENT_BETA_LO)
            ss << " Parent dipole, explicit beta";
        else
            ss << " NO STRING IMPLEMENTED!";

        ss<< ". NLO Kernel RC: ";
        if (RC_NLO == FIXED_NLO or EQUATION==CONFORMAL_N4)
            ss << " fixed as=" << FIXED_AS;
        else if (RC_NLO == SMALLEST_NLO)
            ss << " smallest dipole";
        else if (RC_NLO  == PARENT_NLO)
            ss << " Parent dipole";
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


        if (FORCE_POSITIVE_N)
            ss << ". Amplitude is limited to [0,1].";
        else
            ss << ". Amplitude is not limited!";

        ss << " Alphas scaling C^2=" << config::ALPHAS_SCALING ;
        ss << endl;
        BKSolver sol;
        ss << "# Alphas(r=1 GeV^-1) = " << sol.Alphas(1) << endl;
        ss << "# Order: ";
        if (LO_BK)
            ss <<"LO";
        else
            ss << "NLO";
        if (config::ONLY_NLO) ss << ", keeping only NLO terms";
    
    return ss.str();
}
