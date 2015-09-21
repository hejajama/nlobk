/*
 * nloBK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2013-2014
 */

// Configs

#ifndef _NLOBK_CONFIG_HPP
#define _NLOBK_CONFIG_HPP

#include <string>
#include <sstream>


#define LINEINFO __FILE__ << ":" << __LINE__
inline double SQR(double x) { return x*x; }

namespace config
{
    extern double NC;
    extern double NF;

    

    extern double LAMBDAQCD;
    extern double LAMBDAQCD2;

    extern int RINTPOINTS;
    extern int THETAINTPOINTS;
    extern double INTACCURACY;
    extern double MCINTACCURACY;
    extern double MAXR;
    extern double MINR;
    extern unsigned int RPOINTS;

    extern size_t MCINTPOINTS;

    // Select what to solve
    enum Equation
    {
        QCD,            // QCD NLO BK
        CONFORMAL_QCD,  // Confromal NLO BK in QCD
        CONFORMAL_N4    // Conformal NLO BK in N=4 SYM
    };
    extern Equation EQUATION;  

    extern double DE_SOLVER_STEP;

    // Alpha_s in LO part
    enum RunningCouplingLO
    {
        FIXED_LO,
        PARENT_LO,      // Parent dipole where all beta terms are included
        PARENT_BETA_LO, // Parent dipole, only renormalization scale term is included in as, the second beta term is explicit in the expression
        SMALLEST_LO,
        BALITSKY_LO
    };
    enum RunningCouplingNLO
    {
        FIXED_NLO,
        PARENT_NLO,
        SMALLEST_NLO
    };
    extern double FIXED_AS;


    extern RunningCouplingLO RC_LO;
    extern RunningCouplingNLO RC_NLO;

    extern bool DOUBLELOG_LO_KERNEL; // include double log term from the LO kernel
    extern bool ONLY_DOUBLELOG;     // only include double log term

    enum INTEGRATION_METHOD
    {
        VEGAS,
        MISER,
        MULTIPLE            // No monte carlo
    };
    extern INTEGRATION_METHOD INTMETHOD_NLO;

    extern bool LO_BK ;    // solve only LO BK

    extern bool ONLY_NLO;   // do not keep as^1 terms

    extern bool FORCE_POSITIVE_N;   // Force N(r)>=0

    extern double ALPHAS_SCALING;   // C^2

    extern bool DNDY;   // Print only dn/dy and exit

    extern bool ONLY_LNR;   // Keep only ln r^2 terms from nlo conformal dipole
    extern bool NO_LNR;     // Do not include ln r^2

    extern bool RESUM_DLOG; // Resum double log
    extern bool RESUM_SINGLE_LOG;
    
    enum SINGLELOG_RESUM_RC
    {
		RESUM_RC_BALITSKY,
		RESUM_RC_PARENT,
        RESUM_RC_SMALLEST
	};
	
	extern SINGLELOG_RESUM_RC RESUM_RC;

    extern bool NO_K2;      // Do not include K_2

    extern bool ONLY_RESUM_DLOG; // Calculate contribution from resummation
    
    extern double KSUB;	// Constant factor in the sigle log resummation log
    
    extern bool ONLY_SUBTRACTION;	// Calculate only effect from subtraction

    
}

#endif
