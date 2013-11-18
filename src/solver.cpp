/*
 * nloBK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2013
 */

#include "solver.hpp"
#include "dipole.hpp"

#include "nlobk_config.hpp"

#include <cmath>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h> // odeiv2 Requires GSL 1.15
#include <gsl/gsl_sf_bessel.h>

// Integration constants
const double MAXR=99;
const double MINR=1e-7;
const double eps = 1e-20;

BKSolver::BKSolver(Dipole* d)
{
    dipole=d;
}

struct DEHelper{
    BKSolver* solver;
};

int Evolve(double y, const double amplitude[], double dydt[], void *params);
int BKSolver::Solve(double maxy)
{
    /*
     * In order to be able to use Runge Kutta we create one large array
     * which is basically just vector that we evolve
     * Array size is Dipole->RPoints()
     */

    int vecsize = dipole->RPoints();
    double *ampvec = new double [vecsize];
    dipole->InitializeInterpolation(0); // Initialize interpolation at y=yvals[0]=0
    for (unsigned int rind=0; rind<vecsize; rind++)
    {
        ampvec[rind] = dipole->N( dipole->RVal(rind));
    }

    double y=0; double step = 0.1;  // We have always solved up to y
    int yind=0;

    // Intialize GSL
    DEHelper help; help.solver=this;
    gsl_odeiv_system sys = {Evolve, NULL, vecsize, &help};
        
    const gsl_odeiv_step_type * T = gsl_odeiv_step_rkf45; 
    gsl_odeiv_step * s    = gsl_odeiv_step_alloc (T, vecsize);
    gsl_odeiv_control * c = gsl_odeiv_control_y_new (0.0, 0.01);    //abserr relerr
    gsl_odeiv_evolve * e  = gsl_odeiv_evolve_alloc (vecsize);
    double h = step;  // Initial ODE solver step size
    
    do
    {
        double  nexty = y+step;
        while (y<nexty)
        {
            int status = gsl_odeiv_evolve_apply(e, c, s, &sys,
                &y, nexty, &h, ampvec);
            if (status != GSL_SUCCESS) {
                cerr << "Error in gsl_odeiv_evolve_apply at " << LINEINFO
                << ": " << gsl_strerror(status) << " (" << status << ")"
                << " y=" << y << ", h=" << h << endl;
            }
            cout << "Evolved up to " << y << "/" << nexty << ", h=" << h << endl;
        }

        

        yind = dipole->AddRapidity(y, ampvec);
        // Change Dipole interpolator to the new rapidity
        dipole->InitializeInterpolation(yind);
        cout << "Rapidity index is " << yind << " corresponding to rapidity " << dipole->YVal(yind) << endl;
     
        
    } while (y <= maxy);

    gsl_odeiv_evolve_free (e);
    gsl_odeiv_control_free (c);
    gsl_odeiv_step_free (s);
    delete[] ampvec;
}

int Evolve(double y, const double amplitude[], double dydt[], void *params)
{
    DEHelper* par = reinterpret_cast<DEHelper*>(params);
    Dipole* dipole = par->solver->GetDipole();
    for (unsigned int i=0; i< dipole->RPoints(); i++)
    {
        dydt[i]=0;
    }

    return GSL_SUCCESS;
}

/*
 * LO BK Kernel
 * Evaluates the kernel at given dipole size
 *
 * Parameters:
 *  Parent dipole size r
 *  Daughter dipole size v
 *  Daughter dipole angle theta [0,2\pi]
 */

double BKSolver::Kernel_lo(double r, double v, double theta)
{
    // Y = y-z = v
    double Y=v;
    // X = x-z = r + v
    double X = std::sqrt( r*r + v*v + 2.0*r*v*std::cos(theta) );
    
    
    double result=0;
    result = r*r / (X*X * Y*Y);

    // Handle divergences
    if (X<eps or Y<eps)
        return 0;
    

    return result;
}

/*
 * NLO Kernel, part 1
 *
 * Parent dipole size v
 * 
 * Daughter dipole sizes and angles: (v,theta_v) and (w, theta_w)
 */
double BKSolver::Kernel_nlo_1(double r, double v, double theta_v, double w, double theta_w)
{
    // Y = y-z = v
    double Y=v;
    // X = x-z = r + v
    double X = std::sqrt( r*r + v*v + 2.0*r*v*std::cos(theta_v) );
    // Y' = y-z' = w
    double Y2=w;
    // X' = x-z' = r + w
    double X2 = std::sqrt( r*r + w*w + 2.0*r*w*std::cos(theta_w) );
    // z - z' = w - v
    double z_m_z2 = std::sqrt( w*w + v*v - 2.0*v*w*std::cos(theta_w - theta_v) );

    double kernel=0;

    ///TODO: Check divergences/limits
    kernel = 2.0* ( SQR(X*Y2) + SQR(X2*Y) - 4.0*SQR(z_m_z2 * r)) / ( std::pow(z_m_z2,4.0) * (X*X * Y2*Y2 - X2*X2 * Y*Y) );

    kernel += std::pow(r,4.0) / ( SQR(X*Y2) - SQR(X2*Y)) * (1.0/(X*X * Y2*Y2) + 1.0/(X2*X2 * Y*Y) );

    kernel += SQR(r/z_m_z2) * (1.0/SQR(X*Y2) - 1.0/SQR(X2*Y) );

    kernel *= std::log( SQR(X*Y2/(X2*Y) ) );

    kernel += -4.0 / std::pow(z_m_z2, 4.0);


    if (isinf(kernel) or isnan(kernel))          ///TODO: Kill me
        return 0;


    return kernel;


}

/*
 * NLO Kernel, part 2
 *
 */
 
double BKSolver::Kernel_nlo_2(double r, double v, double theta_v, double w, double theta_w)
{
     // Y = y-z = v
    double Y=v;
    // X = x-z = r + v
    double X = std::sqrt( r*r + v*v + 2.0*r*v*std::cos(theta_v) );
    // Y' = y-z' = w
    double Y2=w;
    // X' = x-z' = r + w
    double X2 = std::sqrt( r*r + w*w + 2.0*r*w*std::cos(theta_w) );
    // z - z' = w - v
    double z_m_z2 = std::sqrt( w*w + v*v - 2.0*v*w*std::cos(theta_w - theta_v) );

    double kernel=0;

    kernel = SQR(r / z_m_z2) * ( 1.0/SQR(X*Y2) + 1.0/SQR(Y*X2) );
    kernel -= std::pow(r,4.0) / SQR(X*Y*X2*Y2);
    kernel *= std::log( SQR(X*Y2/(X2*Y) ) );

    if (isinf(kernel) or isnan(kernel))          ///TODO: Kill me
        return 0;

    return kernel;
}

/*
 * NLO Kernel, part 2
 *
 */
 
double BKSolver::Kernel_nlo_3(double r, double v, double theta_v, double w, double theta_w)
{
    return 0;// TODO
}







Dipole* BKSolver::GetDipole()
{
    return dipole;
}
