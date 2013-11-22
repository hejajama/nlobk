/*
 * nloBK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2013
 */

#ifndef _NLOBK_SOLVER_H
#define _NLOBK_SOLVER_H

#include "dipole.hpp"

/* General solver class for the BK equation
 */

class BKSolver
{
    public:
        BKSolver(Dipole* d);    // Constructor takes the dipole amplitude class
        int Solve(double maxy);	// Solve up to maxy

        double Kernel_lo(double r, double v, double theta);
        double Kernel_nlo_1(double r, double v, double theta_v, double w, double theta_w);
        double Kernel_nlo_2(double r, double v, double theta_v, double w, double theta_w);
        double Kernel_nlo_3(double r, double v, double theta_v, double w, double theta_w);

        double RapidityDerivative_lo(double r);
        double RapidityDerivative_nlo(double r);

        Dipole* GetDipole();
    private:
        Dipole* dipole;

};






#endif
