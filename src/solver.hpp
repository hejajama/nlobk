/*
 * nloBK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2013
 */

#ifndef _NLOBK_SOLVER_H
#define _NLOBK_SOLVER_H

#include "dipole.hpp"
#include <tools/interpolation.hpp>

/* General solver class for the BK equation
 */

class BKSolver
{
    public:
        BKSolver(Dipole* d);    // Constructor takes the dipole amplitude class
        int Solve(double maxy);	// Solve up to maxy

        double Kernel_lo(double r, double v, double theta);
        double Kernel_nlo_1(double r, double X, double Y, double X2, double Y2, double z_m_z2);
        double Kernel_nlo_2(double r, double X, double Y, double X2, double Y2, double z_m_z2);
        double Kernel_nlo_3(double r, double X, double Y, double X2, double Y2, double z_m_z2);

        // kernel for evolution of N
        double Kernel_nlo(double r, double X, double Y, double X2, double Y2, double z_m_z2);

        double RapidityDerivative_lo(double r, Interpolator* dipole_interp);
        double RapidityDerivative_nlo(double r, Interpolator* dipole_interp, Interpolator* dipole_interp_s);

        Dipole* GetDipole();
        
        double Alphas(double r);
    private:
        Dipole* dipole;

};






#endif
