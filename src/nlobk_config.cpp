/*
 * nloBK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2013-2014
 */

#include "nlobk_config.hpp"

// Default configs
namespace config
{
     double NC=3;
     double NF=3;

     double LAMBDAQCD = 0.241;
     double LAMBDAQCD2 = LAMBDAQCD*LAMBDAQCD;

     int RINTPOINTS=25;
     int THETAINTPOINTS = 25;
     double INTACCURACY=0.05;
     double MCINTACCURACY = 0.2;
     double MAXR = 50;
     double MINR=1e-6;
     unsigned int RPOINTS = 150;

     size_t MCINTPOINTS = 2e6;


     Equation EQUATION = QCD;  

     double DE_SOLVER_STEP = 0.1;


     double FIXED_AS = 0.1;


     RunningCouplingLO RC_LO = SMALLEST_LO;
     RunningCouplingNLO RC_NLO = SMALLEST_NLO;

     bool DOUBLELOG_LO_KERNEL = true; // include double log term from the LO kernel

     INTEGRATION_METHOD INTMETHOD_NLO = MISER;

     bool LO_BK = false;    // solve only LO BK

     bool FORCE_POSITIVE_N = true;

     double ALPHAS_SCALING = 1.0;

     bool DNDY=false;

     bool ONLY_NLO = false;
}
