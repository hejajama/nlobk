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
        PARENT_LO,      // Parent dipole 
        PARENT_BETA_LO, // Parent dipole, only renormalization scale term is included in as, the second beta term is explicit in the expression
        SMALLEST_LO,
        BALITSKY_LO,
		FRAC_LO, // fastest apparent convergence in 1507.03651
		GUILLAUME_LO // 1708.06557
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

    extern bool FORCE_POSITIVE_N;   // Force N(r)>=0


    extern bool DNDY;   // Print only dn/dy and exit

    extern bool ONLY_LNR;   // Keep only ln r^2 terms from nlo conformal dipole
    extern bool NO_LNR;     // Do not include ln r^2

    // Resummation options are encoded in the kernel inclusion enum below.
    
    enum SINGLELOG_RESUM_RC
    {
		RESUM_RC_BALITSKY,
		RESUM_RC_PARENT,
        RESUM_RC_SMALLEST,
        RESUM_RC_FIXED
	};
	
	extern SINGLELOG_RESUM_RC RESUM_RC;

    // Kernel inclusion mode: choose which parts of the kernel to include
    enum ORDER
    {
        LO,                         // LO only, no resummation
        LO_RESUM_DLOG,              // LO with double-log resummation
        LO_RESUM_DLOG_SLOG,         // LO with double-log and single-log resummation
        NLO,                        // Full NLO (including K2/Kf), no resummation
        NLO_RESUM_DLOG,             // Full NLO with double-log resummation
        NLO_RESUM_DLOG_SLOG         // Full NLO with double-log and single-log resummation
    };
    extern ORDER Order;

    // Helper function to check if kernel is LO-based (i.e. no explict K2 and Kf NLO parts)
    inline bool IsLOKernel(ORDER kernel) {
        return (kernel == LO || 
                kernel == LO_RESUM_DLOG || 
                kernel == LO_RESUM_DLOG_SLOG);
    }
    
    extern double KSUB;	// Constant factor in the sigle log resummation log
    
    extern bool ONLY_K1FIN; // Only K1fin contribution from K_1
    
    extern bool KINEMATICAL_CONSTRAINT; // Solve nonlocal kinematically constrained BK (LO part)
    
    extern bool EULER_METHOD;    // Use Euler method instead of Runge Kutta, must be true if KINEMATICA_CONSTRAINT is used

    
}
std::string NLOBK_CONFIG_STRING();

#endif
