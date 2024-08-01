/*
 * AmplitudeLib
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011-2014
 */

/** @file */

#ifndef _CONFIG_HPP
#define _CONFIG_HPP 


#include <iostream>
using std::cout;
using std::cerr;
using std::cin;
using std::endl;
#include <string>
using std::string;
#include <cmath>

typedef unsigned int uint;


// Constants in Amplitude namespace, avoid overlapping with other programs
// using this class

namespace Amplitude
{

    const double ALPHA_e = 1.0/137.035999679; 
    //const double e = sqrt(4.0*M_PI*ALPHA_e);

    const double FMGEV = 5.068;





    // Inline functions

    #define LINEINFO __FILE__ << ":" << __LINE__

    inline double SQR(const double x) { return x*x; }
}


#endif
