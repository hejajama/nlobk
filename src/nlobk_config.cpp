/*
 * nloBK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2013-2015  
 */

#include "nlobk_config.hpp"

// Default configs
namespace config
{
     double NC=3;
     double NF=3;

     double LAMBDAQCD = 0.241;
     double LAMBDAQCD2 = LAMBDAQCD*LAMBDAQCD;

     int RINTPOINTS=55;
     int THETAINTPOINTS = 55;
     double INTACCURACY=0.005;
     double MCINTACCURACY = 0.2;
     double MAXR = 50;
     double MINR=1e-6;
     unsigned int RPOINTS = 170;

     size_t MCINTPOINTS = 1e8;


     Equation EQUATION = QCD;  

     double DE_SOLVER_STEP = 0.05;


     double FIXED_AS = 0.2;


     RunningCouplingLO RC_LO = BALITSKY_LO;
     RunningCouplingNLO RC_NLO = PARENT_NLO;

     bool DOUBLELOG_LO_KERNEL = true; // include double log term from the LO kernel
     bool ONLY_DOUBLELOG = false;

     INTEGRATION_METHOD INTMETHOD_NLO = MISER;

     bool LO_BK = false;    // solve only LO BK

     bool FORCE_POSITIVE_N = true;

     double ALPHAS_SCALING = 1.0;

     bool DNDY=false;

     bool ONLY_NLO = false;

     bool ONLY_LNR = false;
     bool NO_LNR = false;

     bool RESUM_DLOG = false;
     bool RESUM_SINGLE_LOG = false;

     bool NO_K2 = false;

     bool ONLY_RESUM_DLOG = false;
     
     bool ONLY_SUBTRACTION = false;
     
     double KSUB = 1.0;
     
     SINGLELOG_RESUM_RC RESUM_RC = RESUM_RC_PARENT;

     bool ONLY_K1FIN = false;
}
