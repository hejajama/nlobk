/*
 * nloBK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2013-2015  
 */

#include "nlobk_config.hpp"
#include <sstream>
#include <string>

using namespace std;

using namespace config;

// Default configs
namespace config
{
     double NC=3;
     double NF=3;

     double LAMBDAQCD = 0.241;
     double LAMBDAQCD2 = LAMBDAQCD*LAMBDAQCD;

     int RINTPOINTS=85;
     int THETAINTPOINTS = 85;
     double INTACCURACY=0.001;
     double MCINTACCURACY = 0.2;
     double MAXR = 30;          // Quite small, only for testing
     double MINR=1e-6;
     unsigned int RPOINTS =100; // Number of points in r grid

     size_t MCINTPOINTS = 1e5;


     Equation EQUATION = QCD;  

    double DE_SOLVER_STEP = 0.2; 
    double DE_SOLVER_ABSERR = 1e-6;
    double DE_SOLVER_RELERR = 1e-4;

     double FIXED_AS = 0.2;


     RunningCouplingLO RC_LO = BALITSKY_LO;
     RunningCouplingNLO RC_NLO = SMALLEST_NLO;
     SINGLELOG_RESUM_RC RESUM_RC = RESUM_RC_SMALLEST; // Resummation running coupling

     INTEGRATION_METHOD INTMETHOD_NLO = VEGAS;

     bool FORCE_POSITIVE_N = true;
     bool SYMMETRIZE_Z_Z2_INTEGRATION = true;

     bool DNDY=false;

     // Kernel inclusion default: full NLO with double+single log resummation
    ORDER Order = NLO_RESUM_DLOG_SLOG;
    
     
     double KSUB = 1.0;
     
         
    bool KINEMATICAL_CONSTRAINT = false;
    
    bool EULER_METHOD = false;
}


std::string NLOBK_CONFIG_STRING()
{
    std::stringstream ss;
    
    
    ss << "MC integration method: ";
    if (INTMETHOD_NLO == MISER)
    ss <<"MonteCarlo Miser, points=" << MCINTPOINTS;
    else if (INTMETHOD_NLO == VEGAS)
    ss <<"MonteCarlo Vegas, points=" << MCINTPOINTS;
    else if (INTMETHOD_NLO == MULTIPLE)
    ss << "Multiple integrals (no montecarlo)";
    else
    ss <<"UNKNOWN!";
    ss << ". K1 integration accuracy " << INTACCURACY ;
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
        // Double-log behavior is controlled by kernel inclusion enum now.
        ss << ". QCD";
    }
    else if (EQUATION == CONFORMAL_QCD) ss << ". Solving for CONFORMAL dipole";
    else if (EQUATION == CONFORMAL_N4) ss << ". Solving in N=4 for CONFORMAL dipole";
    else ss << ". UNKNOWN EQUATION!!";
    
    
    if (FORCE_POSITIVE_N)
    ss << ". Amplitude is limited to [0,1].";
    else
    ss << ". Amplitude is not limited!";
    
    ss << endl;
    //BKSolver sol;
    //ss << "# Alphas(r=1 GeV^-1) = " << sol.Alphas(1) << endl;
    ss << "# Order: ";
    bool is_lo = (config::Order == config::LO || config::Order == config::LO_RESUM_DLOG || config::Order == config::LO_RESUM_DLOG_SLOG);
    if (is_lo)
    ss <<"LO";
    else
    ss << "NLO";
    // Print resummation info based on kernel inclusion
    bool resum_dlog = (config::Order == config::LO_RESUM_DLOG || config::Order == config::LO_RESUM_DLOG_SLOG
                       || config::Order == config::NLO_RESUM_DLOG || config::Order == config::NLO_RESUM_DLOG_SLOG);
    bool resum_slog = (config::Order == config::LO_RESUM_DLOG_SLOG || config::Order == config::NLO_RESUM_DLOG_SLOG);
    if (resum_dlog)
    {
        ss << endl << "# Resumming double log";
    }
    if (resum_slog)
    {
        ss << endl << "# Resumming single log, K_sub=" << config::KSUB;
        if (config::RESUM_RC == RESUM_RC_PARENT) ss << " resum rc: parent";
        else if (config::RESUM_RC == RESUM_RC_SMALLEST) ss << " resum rc: smallest";
        else if (config::RESUM_RC == RESUM_RC_BALITSKY) ss << " resum rc: balitsky";
        ss << endl;
    }
    
      
    ss << "# Kernel : ";
    switch (config::Order)
    {
        case LO: ss << "LO only"; break;
        case LO_RESUM_DLOG: ss << "LO with double-log resummation"; break;
        case LO_RESUM_DLOG_SLOG: ss << "LO with double+single log resummation"; break;
        case NLO: ss << "NLO (K2/Kf included), no resummation"; break;
        case NLO_RESUM_DLOG: ss << "NLO with double-log resummation"; break;
        case NLO_RESUM_DLOG_SLOG: ss << "NLO with double+single log resummation"; break;
    }
    ss << endl;
    
    if (config::KINEMATICAL_CONSTRAINT)
        ss << endl << "# Kinematical constraint included" << endl;
    return ss.str();
}

