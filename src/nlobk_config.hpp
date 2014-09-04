/*
 * nloBK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2013-2014
 */

// Configs

#ifndef _NLOBK_CONFIG_HPP
#define _NLOBK_CONFIG_HPP

#include <string>
#include <sstream>

const double NC=3;
const double NF=0;

inline double SQR(double x) { return x*x; }
#define LINEINFO __FILE__ << ":" << __LINE__

const double LAMBDAQCD = 0.241;
const double LAMBDAQCD2 = LAMBDAQCD*LAMBDAQCD;

const int RINTPOINTS=20;
const int THETAINTPOINTS = 20;
const double INTACCURACY=0.1;
const double MCINTACCURACY = 0.2;
const double MAXR = 40;
const double MINR=1e-5;
const unsigned int RPOINTS = 100;

const size_t MCINTPOINTS = 5e5;

// Select what to solve
enum Equation
{
    QCD,            // QCD NLO BK
    CONFORMAL_QCD,  // Confromal NLO BK in QCD
    CONFORMAL_N4    // Conformal NLO BK in N=4 SYM
};
const Equation EQUATION = CONFORMAL_N4;  

const double DE_SOLVER_STEP = 0.5;

// Alpha_s in LO part
enum RunningCouplingLO
{
    FIXED_LO,
    PARENT_LO,
    SMALLEST_LO,
    BALITSKY_LO
};
enum RunningCouplingNLO
{
    FIXED_NLO,
    PARENT_NLO,
    SMALLEST_NLO
};
const double FIXED_AS = 0.05;


const RunningCouplingLO RC_LO = SMALLEST_LO;
const RunningCouplingNLO RC_NLO = SMALLEST_NLO;

const bool DOUBLELOG_LO_KERNEL = true; // include double log term from the LO kernel

enum INTEGRATION_METHOD
{
    VEGAS,
    MISER,
    MULTIPLE            // No monte carlo
};
const INTEGRATION_METHOD INTMETHOD_NLO = VEGAS;

    

#endif
