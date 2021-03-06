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
#include <sstream>

using namespace std;

string version = "0.01-dev";

using namespace config;



std::stringstream cmd;
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
    
    cout <<"#"<<endl<<"# NLOBK solver " << version  << " running on " << hostname << " build on " << __DATE__ << " " << __TIME__ << endl;
    cout <<"# Start at " << today << "#" << endl;
    delete[] hostname;

    cout << "# Command: ";    
    for (int i=0; i<argc; i++)
        cmd << argv[i] << " ";
    cout << cmd.str() << endl; cout << "#" << endl;
    
    double alphas_scaling = 1;

    if (string(argv[1])=="-help")
    {
        cout <<"-ic FILE filename: set initial condition (default: MV)" << endl;
        cout << "-ic PARAM qsqr anomalous_dimension ln(e_c)" << endl; 
        cout << "-maxy yval" << endl;
        cout <<"-output filename_to_save_data" << endl;
        cout << "-eq qcd,confqcd,n4: set equation to solve" << endl;
        cout << "-rc smallest,parent,parent_beta,fixed: set running coupling" << endl;
        cout << "-lo: solve LO BK" << endl;
        cout << "-only_nlo: keep only nlo terms" << endl;
        cout << "-nodlog: do not include double log term" << endl;
        cout << "-onlydlog: only include double log term" << endl;
        cout << "-onlylnr: only include lnr^2 NLO terms" << endl;
        cout << "-nolnr: do not include ln r^2 terms" << endl;
        cout << "-nf nf: set number of quark flavors" << endl;
        cout << "-nolimit: do not force N>=0" << endl;
        cout << "-ic [mv,mve,mvgamma] set initial condition" << endl;
        cout << "-dndy: print dn/dy at initial condition and exit" << endl;
        cout << "-alphas_scaling C^2: set C^2 [setting mv/mve/mvgamma ic sets this also]" << endl;
        cout << "-ln_alphas_scaling ln C^2: set ln C^2" << endl;
        cout << "-resum_dlog: resum double log when solving non-conformal dipole" << endl;
        cout << "-resum_slog: resum single log" << endl;
        cout << "-Ksub value: K_sub for single log resummation" << endl;
        cout << "-only_subtraction: calculate only effect from subtraction" << endl;
        cout << "-no_k2: do not include K_2 and K_f" << endl;
        cout << "-ONLY_RESUM_DLOG: only calculate the effect of resummation" << endl;
        cout << "-only_k1fin: include only k1fin contribution from k1"; 
        cout << endl;
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
            else if (string(argv[i+1])=="PARAM")
            {
                ic = new MV();
                ((MV*)ic)->SetQsqr(StrToReal(argv[i+2]));
                ((MV*)ic)->SetAnomalousDimension(StrToReal(argv[i+3]));
				((MV*)ic)->SetE(exp(StrToReal(argv[i+4])));
                alphas_scaling = 1; // This is set by -alphas_scaling, hopefully
            }
            else
            {
                ic = new MV();
                if (string(argv[i+1])=="mv")
                {
                    ((MV*)ic)->SetQsqr(0.2);
                    alphas_scaling=0.315237;
                }
                else if (string(argv[i+1])=="mv02")
                {
                    ((MV*)ic)->SetQsqr(0.2);
                    alphas_scaling = 1;
                }
                else if (string(argv[i+1])=="mv1")
                {
                    ((MV*)ic)->SetQsqr(1);
                    alphas_scaling = std::exp(-2.0*0.57721);
                }
                else if (string(argv[i+1])=="mv10")
                {
                    ((MV*)ic)->SetQsqr(10);
                    alphas_scaling = std::exp(-2.0*0.57721);
                }
                else if (string(argv[i+1])=="mve")
                {
                    ((MV*)ic)->SetQsqr(0.06);
                    ((MV*)ic)->SetE(18.9);
                    alphas_scaling = 7.2;
                }
                else if (string(argv[i+1])=="mvgamma")
                {
                    ((MV*)ic)->SetQsqr(0.165);
                    ((MV*)ic)->SetAnomalousDimension(1.135);
                    alphas_scaling = 6.35;
                }
                else if (string(argv[i+1])=="mvgamma_08")
                {
                    ((MV*)ic)->SetQsqr(1);
                    ((MV*)ic)->SetAnomalousDimension(0.8);
                    alphas_scaling = std::exp(-2.0*0.57721);
                }
                else if (string(argv[i+1])=="mvgamma_08_qsqr_100")
                {   
                    ((MV*)ic)->SetQsqr(100);
                    ((MV*)ic)->SetAnomalousDimension(0.8);
                    alphas_scaling = std::exp(-2.0*0.57721);
                }
                else if (string(argv[i+1])=="mvgamma_09")
                {
                    ((MV*)ic)->SetQsqr(1);
                    ((MV*)ic)->SetAnomalousDimension(0.9);
                    alphas_scaling = std::exp(-2.0*0.57721);
                }
                else if (string(argv[i+1])=="mvgamma_095")
                {
                    ((MV*)ic)->SetQsqr(1);
                    ((MV*)ic)->SetAnomalousDimension(0.95);
                    alphas_scaling = std::exp(-2.0*0.57721);
                }
                else if (string(argv[i+1])=="mvgamma_1")
                {
                    ((MV*)ic)->SetQsqr(1);
                    ((MV*)ic)->SetAnomalousDimension(1.0);
                    alphas_scaling = std::exp(-2.0*0.57721);
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
                config::RESUM_RC = config::RESUM_RC_SMALLEST;
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
                config::RESUM_RC = config::RESUM_RC_PARENT;
            }
			else if (string(argv[i+1])=="guillaume")
			{
				config::RC_LO = config::GUILLAUME_LO;
				config::RC_NLO = config::PARENT_NLO;
				config::RESUM_RC = config::RESUM_RC_PARENT;
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
        else if (string(argv[i])=="-nodlog")
        {
            config::DOUBLELOG_LO_KERNEL = false;
            config::ONLY_NLO = true;
        }
        else if (string(argv[i])=="-onlydlog")
        {
            config::ONLY_DOUBLELOG = true;
            config::ONLY_NLO = true;
        }
        else if (string(argv[i])=="-nlo")
            config::LO_BK = false;
        else if (string(argv[i])=="-nf")
        {
            config::NF = StrToInt(argv[i+1]);
            if (config::NF != 0 and config::NF != 3 and config::NF != 5)
            {
                cerr << "Invalid Nf=" << config::NF << " " << LINEINFO << endl;
                return -1;
            }
        }            
        
        else if (string(argv[i])=="-nolimit")
            config::FORCE_POSITIVE_N = false;

        else if (string(argv[i])=="-dndy")
            config::DNDY = true;

        else if (string(argv[i])=="-onlylnr")
        {
            config::ONLY_LNR = true;
            config::ONLY_NLO = true;
        }
        else if (string(argv[i])=="-nolnr")
        {
            config::NO_LNR = true;
            config::ONLY_NLO = true;
        }
        else if (string(argv[i])=="-alphas_scaling")
            alphas_scaling = StrToReal(argv[i+1]);
		else if (string(argv[i])=="-ln_alphas_scaling")
			alphas_scaling = exp(StrToReal(argv[i+1]));
        else if (string(argv[i])=="-ln_alphas_scaling")
            alphas_scaling = std::exp(StrToReal(argv[i+1]));
        
        else if (string(argv[i])=="-resum_dlog")
            config::RESUM_DLOG = true;

        else if (string(argv[i])=="-resum_slog")
            config::RESUM_SINGLE_LOG = true;
        else if (string(argv[i])=="-Ksub")
			config::KSUB = StrToReal(argv[i+1]);
		
        else if (string(argv[i])=="-no_k2")
            config::NO_K2 = true;
        
        else if (string(argv[i])=="-only_subtraction")
			config::ONLY_SUBTRACTION = true;

        else if (string(argv[i])=="-ONLY_RESUM_DLOG")
        {
            config::ONLY_RESUM_DLOG = true;
            config::RESUM_DLOG = true;
            config::NO_K2 = true;
        }
        else if (string(argv[i])=="-only_k1fin")
            config::ONLY_K1FIN = true;
        
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

    if (config::EQUATION != config::QCD and config::RESUM_DLOG == true)
    {
        cerr << "Asked to resum dlog and solve something else than qcd!" << endl;
        exit(1);
    }

    if ( (config::RC_LO == FIXED_LO and config::RC_NLO != FIXED_NLO)
        or (config::RC_NLO == FIXED_NLO and config::RC_LO != FIXED_LO) )
    {
            cerr << "Must have fixed coupling in both LO and NLO!" << endl;
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
    solver.SetAlphasScaling(alphas_scaling);
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

