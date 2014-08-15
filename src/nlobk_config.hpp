/*
 * nloBK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2013
 */

// Configs

#ifndef _NLOBK_CONFIG_HPP
#define _NLOBK_CONFIG_HPP

const double NC=3;
const double NF=0;

inline double SQR(double x) { return x*x; }
#define LINEINFO __FILE__ << ":" << __LINE__

const double LAMBDAQCD = 0.241;
const double LAMBDAQCD2 = LAMBDAQCD*LAMBDAQCD;

const int RINTPOINTS=20;
const int THETAINTPOINTS = 20;
const double INTACCURACY=0.1;
const double MAXR = 100;

const bool MONTECARLO = true;
const size_t MCINTPOINTS = 1e6;

const double DE_SOLVER_STEP = 0.5;

// Alpha_s in LO part
enum RunningCouplingLO
{
    FIXED_LO,
    PARENT_LO,
    BALITSKY_LO
};
enum RunningCouplingNLO
{
    FIXED_NLO,
    PARENT_NLO
};
const double FIXED_AS = 0.05;


const RunningCouplingLO RC_LO = FIXED_LO;
const RunningCouplingNLO RC_NLO = FIXED_NLO;

const bool DOUBLELOG_LO_KERNEL = true; // include double log term from the LO kernel

#endif
