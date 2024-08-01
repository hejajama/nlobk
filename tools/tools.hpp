/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#ifndef _TOOLS_H
#define _TOOLS_H

#include <string>
#include <vector>
#include "config.hpp"

double StrToReal(std::string str);
int StrToInt(std::string str);
// GSL error handler
void ErrHandler(const char * reason,
                        const char * file,
                        int line,
                        int gsl_errno);

/*
double Alpha_s(double Qsqr, double scaling=1.0);
double Alpha_s_r(double rsqr, double scaling=1.0);
double Alphabar_s(double Qsqr, double scaling=1.0); // \alpha_s N_C / Pi
*/
int FindIndex(double val, std::vector<double> &array);

std::string Alpha_s_str();

// Running coupling
const double MAXALPHA = 0.7;


#endif

