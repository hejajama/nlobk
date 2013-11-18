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

#endif
