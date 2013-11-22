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

const int RINTPOINTS=5;
const int THETAINTPOINTS = 3;
const double INTACCURACY=0.3;

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

    double y=0; double step = 0.2;  // We have always solved up to y
    int yind=0;

    // Intialize GSL
    DEHelper help; help.solver=this;
    gsl_odeiv_system sys = {Evolve, NULL, vecsize, &help};
        
    const gsl_odeiv_step_type * T = gsl_odeiv_step_rk2; // rkf45 is more accurate 
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
	cout << "Evolving, rapidity " << y << ", rpoints: " << dipole->RPoints() << endl;
	#pragma omp parallel for
    for (unsigned int i=0; i< dipole->RPoints(); i++)
    {
        double lo = par->solver->RapidityDerivative_lo(dipole->RVal(i));;
        double nlo = par->solver->RapidityDerivative_nlo(dipole->RVal(i));
        cout << "*Dipole size " << dipole->RVal(i) << " lo " << lo << " nlo " << nlo << endl;
        dydt[i]= lo + nlo;
    }

    return GSL_SUCCESS;
}



/*
 * *****************
 * LO PART
 * *****************
 */
 
/*
 * Rapidity derivatives
 * compute \partial_y N(r)
 * First leading order contribution (only one 2d integral), then nlo correction
 *
 * Note: here we assume that dipole amplitude is initialized at correct
 * rapidity so that it can just be evaluated
 */

struct Inthelper_nlobk
{
    BKSolver* solver;
    double r;   // = x - y = parent dipole
    double v;   // = y - z = daughter dipole
    double theta_v; // direction of v
    double w;   // = y - z' = daughter dipole 2
    double theta_w; // direction of w
};

double Inthelperf_lo_v(double v, void* p);
double Inthelperf_lo_theta(double theta, void* p);

double BKSolver::RapidityDerivative_lo(double r)
{
    gsl_function fun;
    Inthelper_nlobk helper;
    helper.solver=this;
    helper.r=r;
    
    fun.params = &helper;
    fun.function = Inthelperf_lo_v;

    gsl_integration_workspace *workspace 
     = gsl_integration_workspace_alloc(RINTPOINTS);

    double minlnr = std::log( 0.5*dipole->MinR() );
    double maxlnr = std::log( 2.0*dipole->MaxR() );

    int status; double  result, abserr;
    status=gsl_integration_qag(&fun, minlnr,
            maxlnr, 0, INTACCURACY, RINTPOINTS,
            GSL_INTEG_GAUSS15, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);

    if (status)
    {
        //#pragma omp critical
        //cerr << "RInt failed, r=" << r <<": at " << LINEINFO << ", result " << result << ", relerr "
            //<< std::abs(abserr/result) << endl;
    }

    result *= 0.2/(2.0*M_PI); // alphas*nc/(2pi^2), alphas*Nc/pi=0.2;
    
    return result;
}

double Inthelperf_lo_v(double v, void* p)
{
    Inthelper_nlobk* helper = reinterpret_cast<Inthelper_nlobk*>(p);
    helper->v=std::exp(v);
    //cout << " r " << helper->r << " v " << helper->v << endl;
    gsl_function fun;
    fun.function=Inthelperf_lo_theta;
    fun.params = helper;
    
    gsl_integration_workspace *workspace 
     = gsl_integration_workspace_alloc(THETAINTPOINTS);

    int status; double result, abserr;
    status=gsl_integration_qag(&fun, 0,
            M_PI, 0, INTACCURACY, THETAINTPOINTS,
            GSL_INTEG_GAUSS15, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);

    if (status)
    {
        //#pragma omp critical
        //cerr << "RInt failed, v=" << v <<", r=" << helper->r <<": at " << LINEINFO << ", result " << result << ", relerr "
        //<< std::abs(abserr/result) << endl;
    }

    result *= std::exp(2.0*v);  // Jacobian v^2 dv
    result *= 2.0;  // As integration limits are just [0,pi]
    return result;
}

double Inthelperf_lo_theta(double theta, void* p)
{
        //cout << " theta " << theta << endl;
    Inthelper_nlobk* helper = reinterpret_cast<Inthelper_nlobk*>(p);
    double r = helper->r;   // x - y
    double v = helper->v;

    // X = x - z = r + v
    double X = std::sqrt( r*r + v*v + 2.0*r*v*std::cos(theta));
    // Y = y - z = v
    double Y = v;

    Dipole* d = helper->solver->GetDipole();

    double N_X = d->N(X);    // N(X)
    double N_Y = d->N(Y);
    double N_r = d->N(r);

    //cout << N_X << " " << N_Y << " " << N_r << endl;

    return helper->solver->Kernel_lo(r, v, theta) * ( N_X + N_Y - N_r - N_X*N_Y );
     
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
 * *****************
 * NLO PART
 * *****************
 */

double Inthelperf_nlo_v(double v, void* p);
double Inthelperf_nlo_theta_v(double theta, void* p);
double Inthelperf_nlo_w(double v, void* p);
double Inthelperf_nlo_theta_w(double theta, void* p);
double Inthelperf_nlo(double r, double v, double theta_v, double w, double theta_w, BKSolver* solver);

double BKSolver::RapidityDerivative_nlo(double r)
{
    gsl_function fun;
    Inthelper_nlobk helper;
    helper.solver=this;
    helper.r=r;
    
    fun.params = &helper;
    fun.function = Inthelperf_nlo_v;

    gsl_integration_workspace *workspace 
     = gsl_integration_workspace_alloc(RINTPOINTS);

    double minlnr = std::log( 0.5*dipole->MinR() );
    double maxlnr = std::log( 2.0*dipole->MaxR() );

    int status; double  result, abserr;
    status=gsl_integration_qag(&fun, minlnr,
            maxlnr, 0, INTACCURACY, RINTPOINTS,
            GSL_INTEG_GAUSS15, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);

    if (status)
    {
        //#pragma omp critical
        //cerr << "RInt failed, r=" << r <<": at " << LINEINFO << ", result " << result << ", relerr "
            //<< std::abs(abserr/result) << endl;
    }

    result *= 0.2/(2.0*M_PI); // alphas*nc/(2pi^2), alphas*Nc/pi=0.2;
    
    return result;
}

double Inthelperf_nlo_v(double v, void* p)
{
    Inthelper_nlobk* helper = reinterpret_cast<Inthelper_nlobk*>(p);
    helper->v=std::exp(v);
    gsl_function fun;
    fun.function=Inthelperf_nlo_theta_v;
    fun.params = helper;
    
    gsl_integration_workspace *workspace 
     = gsl_integration_workspace_alloc(THETAINTPOINTS);

    int status; double result, abserr;
    status=gsl_integration_qag(&fun, 0,
            2.0*M_PI, 0, INTACCURACY, THETAINTPOINTS,
            GSL_INTEG_GAUSS15, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);

    if (status)
    {
        //#pragma omp critical
        //cerr << "RInt failed, v=" << v <<", r=" << helper->r <<": at " << LINEINFO << ", result " << result << ", relerr "
        //<< std::abs(abserr/result) << endl;
    }

    result *= std::exp(2.0*v);  // Jacobian v^2 dv
    return result;
}

double Inthelperf_nlo_theta_v(double theta, void* p)
{
    Inthelper_nlobk* helper = reinterpret_cast<Inthelper_nlobk*>(p);
    helper->theta_v=theta;

    gsl_integration_workspace *workspace 
     = gsl_integration_workspace_alloc(RINTPOINTS);

    gsl_function fun;
    fun.function=Inthelperf_nlo_w;
    fun.params=helper;

    double minlnr = std::log( 0.5*helper->solver->GetDipole()->MinR() );
    double maxlnr = std::log( 2.0*helper->solver->GetDipole()->MaxR() );

    int status; double  result, abserr;
    status=gsl_integration_qag(&fun, minlnr,
            maxlnr, 0, INTACCURACY, RINTPOINTS,
            GSL_INTEG_GAUSS15, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);

    if (status)
    {
        //#pragma omp critical
        //cerr << "RInt failed, r=" << r <<": at " << LINEINFO << ", result " << result << ", relerr "
            //<< std::abs(abserr/result) << endl;
    }

    return result;

     
}

double Inthelperf_nlo_w(double w, void* p)
{
    Inthelper_nlobk* helper = reinterpret_cast<Inthelper_nlobk*>(p);
    helper->w=std::exp(w);
    gsl_function fun;
    fun.function=Inthelperf_nlo_theta_w;
    fun.params = helper;
    
    gsl_integration_workspace *workspace 
     = gsl_integration_workspace_alloc(THETAINTPOINTS);

    int status; double result, abserr;
    status=gsl_integration_qag(&fun, 0,
            2.0*M_PI, 0, INTACCURACY, THETAINTPOINTS,
            GSL_INTEG_GAUSS15, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);

    if (status)
    {
        //#pragma omp critical
        //cerr << "RInt failed, v=" << v <<", r=" << helper->r <<": at " << LINEINFO << ", result " << result << ", relerr "
        //<< std::abs(abserr/result) << endl;
    }

    result *= std::exp(2.0*w);  // Jacobian v^2 dv
    return result;
}

double Inthelperf_nlo_theta_w(double theta_w, void* p)
{
    Inthelper_nlobk* helper = reinterpret_cast<Inthelper_nlobk*>(p);
    double result = Inthelperf_nlo( helper->r, helper->v, helper->theta_v, helper->w, theta_w, helper->solver);
    return result;
}


double Inthelperf_nlo(double r, double v, double theta_v, double w, double theta_w, BKSolver* solver)
{
    Dipole* d = solver->GetDipole();
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
    
    double result=0;

    result = solver->Kernel_nlo_1(r,v,theta_v,w,theta_w)
        * (d->S(X) * d->S(z_m_z2) * d->S(Y2) - d->S(X)*d->S(Y));

    result += solver->Kernel_nlo_2(r,v,theta_v,w,theta_w) *
        d->S(X)*d->S(z_m_z2)*d->S(Y2);


    result *= -1.0;     // As NLOBK is writen for S=1-N, and we solve it for N

    result *= 0.2*0.2 / ( 16.0*M_PI);   // alpha_s^2 Nc^2/(16pi^2), alphas*nc/pi=0.2

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
