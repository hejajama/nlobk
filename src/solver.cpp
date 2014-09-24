/*
 * nloBK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2013-2014
 */

#include "solver.hpp"
#include "dipole.hpp"

#include "nlobk_config.hpp"

#include <cmath>
#include <ctime>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h> // odeiv2 Requires GSL 1.15
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_plain.h>

// Integration constants
const double eps = 1e-20;

BKSolver::BKSolver(Dipole* d)
{
    dipole=d;
    tmp_output = "";
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

    cout <<"#### Solving BK equation up to y=" << maxy <<", mcintpoints " << MCINTPOINTS << endl;
    cout << "# Nc=" << NC << ", Nf=" << NF << " alphas(r=1) = " << Alphas(1) << endl;
    

    int vecsize = dipole->RPoints();
    double *ampvec = new double [vecsize];
    dipole->InitializeInterpolation(0); // Initialize interpolation at y=yvals[0]=0
    for (unsigned int rind=0; rind<vecsize; rind++)
    {
        ampvec[rind] = dipole->N( dipole->RVal(rind));
    }

    double y=0; double step = DE_SOLVER_STEP;  // We have always solved up to y
    int yind=0;

    // Intialize GSL
    DEHelper help; help.solver=this;
    gsl_odeiv_system sys = {Evolve, NULL, vecsize, &help};
        
    const gsl_odeiv_step_type * T = gsl_odeiv_step_rk2; // rkf45 is more accurate 
    gsl_odeiv_step * s    = gsl_odeiv_step_alloc (T, vecsize);
    gsl_odeiv_control * c = gsl_odeiv_control_y_new (0.0001, INTACCURACY);    //abserr relerr
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

        if (tmp_output != "")
            dipole->Save(tmp_output);
        
        // Change Dipole interpolator to the new rapidity
        dipole->InitializeInterpolation(yind);
     
        
    } while (y < maxy);

    gsl_odeiv_evolve_free (e);
    gsl_odeiv_control_free (c);
    gsl_odeiv_step_free (s);
    delete[] ampvec;
    return 0;
}

int Evolve(double y, const double amplitude[], double dydt[], void *params)
{
    DEHelper* par = reinterpret_cast<DEHelper*>(params);
    Dipole* dipole = par->solver->GetDipole();
	//cout << "#Evolving, rapidity " << y << ", rpoints: " << dipole->RPoints() << endl;


    // Create interpolators for N(r) and S(r)=1-N(r)
    std::vector<double> rvals,nvals,yvals_s;
    double maxr_interp=-1; // at r>maxr N==1 and S==0
    for (int i=0; i<dipole->RPoints(); i++)
    {
        rvals.push_back(dipole->RVal(i));

        double n = amplitude[i];
        if (n>1.0) n=1.0;
        if (n<0) n=0;
        nvals.push_back(n);

        if (amplitude[i]>0.9999)
            maxr_interp = dipole->RVal(i);
        
        double s = 1.0-amplitude[i];
        if (s<0) s=0;
        if (s>1) s=1.0;
        yvals_s.push_back(s);
    }
    Interpolator interp(rvals,nvals);
    interp.Initialize();
    interp.SetFreeze(true);
    interp.SetUnderflow(0);
    interp.SetOverflow(1.0);

    Interpolator interp_s(rvals,yvals_s);
    interp_s.Initialize();
    interp_s.SetFreeze(true);
    interp_s.SetUnderflow(1.0);
    interp_s.SetOverflow(0.0);


    int ready=0;
	#pragma omp parallel for
    for (unsigned int i=0; i< dipole->RPoints(); i+=1)
    {
        //#pragma omp critical
        //cout <<"# r=" << dipole->RVal(i) << endl;

        if (amplitude[i]==0)
        {
            // An ugly hack: the amplitude is froze to zero, and |dndy| becomes large,
            // thus the de solver tries to make the step size even smaller(???)
            dydt[i]=0;
            continue;
        }


        double lo = par->solver->RapidityDerivative_lo(dipole->RVal(i), &interp);

        double nlo=0;
        if (!LO_BK)
            nlo = par->solver->RapidityDerivative_nlo(dipole->RVal(i), &interp, &interp_s);


        //#pragma omp critical
            //cout << dipole->RVal(i) << " " << lo << " " << nlo << " " << amplitude[i] << endl;
        dydt[i]= lo + nlo;

        #pragma omp critical
        {
            ready++;
            if (ready%50==0)
                cout << "# y=" << y <<", ready " << ready << " / " << dipole->RPoints() << endl;
        }
        
    }
    //exit(1);
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
    double z;   // = y - z = daughter dipole
    double theta_z; // direction of v
    double z2;   // = y - z' = daughter dipole 2
    double theta_z2; // direction of w
    Interpolator* dipole_interp;
    Interpolator* dipole_interp_s;  // interpolates S=1-N
};

double Inthelperf_lo_z(double v, void* p);
double Inthelperf_lo_theta(double theta, void* p);

double BKSolver::RapidityDerivative_lo(double r, Interpolator* dipole_interp)
{
    gsl_function fun;
    Inthelper_nlobk helper;
    helper.solver=this;
    helper.r=r;
    helper.dipole_interp = dipole_interp;
    
    fun.params = &helper;
    fun.function = Inthelperf_lo_z;

    gsl_integration_workspace *workspace 
     = gsl_integration_workspace_alloc(RINTPOINTS);

    double minlnr = std::log( 0.5*dipole->MinR() );
    double maxlnr = std::log( 2.0*dipole->MaxR() );

    int status; double  result, abserr;
    status=gsl_integration_qag(&fun, minlnr,
            maxlnr, 0, INTACCURACY, RINTPOINTS,
            GSL_INTEG_GAUSS21, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);

    if (status)
    {
        //#pragma omp critical
        //cerr << "RInt failed, r=" << r <<": at " << LINEINFO << ", result " << result << ", relerr "
            //<< std::abs(abserr/result) << endl;
    }

    
    return result;
}

double Inthelperf_lo_z(double z, void* p)
{
    Inthelper_nlobk* helper = reinterpret_cast<Inthelper_nlobk*>(p);
    helper->z=std::exp(z);
    //cout << " r " << helper->r << " v " << helper->v << endl;
    gsl_function fun;
    fun.function=Inthelperf_lo_theta;
    fun.params = helper;
    
    gsl_integration_workspace *workspace 
     = gsl_integration_workspace_alloc(THETAINTPOINTS);

    int status; double result, abserr;
    status=gsl_integration_qag(&fun, 0,
            2.0*M_PI, 0, INTACCURACY, THETAINTPOINTS,
            GSL_INTEG_GAUSS21, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);

    if (status)
    {
        //#pragma omp critical
        //cerr << "RInt failed, v=" << v <<", r=" << helper->r <<": at " << LINEINFO << ", result " << result << ", relerr "
        //<< std::abs(abserr/result) << endl;
    }

    result *= std::exp(2.0*z);  // Jacobian v^2 dv
    //result *= 2.0;  // As integration limits are just [0,pi]
    return result;
}

double Inthelperf_lo_theta(double theta, void* p)
{
        //cout << " theta " << theta << endl;
    Inthelper_nlobk* helper = reinterpret_cast<Inthelper_nlobk*>(p);
    double r = helper->r;   // x - y
    double z = helper->z;

    // X = x - z = r - z
    double X = std::sqrt( r*r + z*z - 2.0*r*z*std::cos(theta));
    // Y = y - z = z
    double Y = z;




    double N_X = helper->dipole_interp->Evaluate(X);    // N(X)
    double N_Y = helper->dipole_interp->Evaluate(Y);
    double N_r = helper->dipole_interp->Evaluate(r);

    //cout << N_X << " " << N_Y << " " << N_r << endl;

    //cout << z << " " << theta << " " << helper->solver->Kernel_lo(r, z, theta) * ( N_X + N_Y - N_r - N_X*N_Y ) << " " << helper->solver->Kernel_lo(r, z, theta)  << " " <<  ( N_X + N_Y - N_r - N_X*N_Y ) <<  endl;


    return helper->solver->Kernel_lo(r, z, theta) * ( N_X + N_Y - N_r - N_X*N_Y );
     
}


/*
 * LO BK Kernel
 * Evaluates the kernel at given dipole size
 *
 * Parameters:
 *  Parent dipole size r
 *  Daughter dipole size z
 *  Daughter dipole angle theta [0,2\pi]
 */

double BKSolver::Kernel_lo(double r, double z, double theta)
{
    const double musqr = 0;   // scale in NLO part of LO kernel
    
    // Y = y-z = z
    double Y=z;
    // X = x-z = r - z
    double X = std::sqrt( r*r + z*z - 2.0*r*z*std::cos(theta) );
    
    
    double result=0;

    // Handle divergences
    if (X<eps or Y<eps)
        return 0;

    // N=4 is easy as the coupling does not run
    if (EQUATION == CONFORMAL_N4)
    {
        result = FIXED_AS / (2.0*M_PI*M_PI) * r*r / (X*X * Y*Y);
        if (LO_BK)
            return result;
        else
            result *= (1.0 - FIXED_AS * NC / (4.0*M_PI) * M_PI*M_PI/3.0);
        return result;
    }

    ////// QCD

    // Running coupling
    if (RC_LO == FIXED_LO)
    {
        result = FIXED_AS * NC/(2.0*M_PI*M_PI) *  r*r / (X*X * Y*Y);

        if (LO_BK)
            return result;
        
        if (EQUATION == CONFORMAL_QCD)   ///TODO: should we include ~beta terms here or not?
        {
            result *= (1.0 + FIXED_AS*NC/(4.0*M_PI) * (67.0/9.0 - SQR(M_PI)/3.0) );
        }
        else if (EQUATION == QCD)
        {
            result *= (1.0 + FIXED_AS * NC / (4.0*M_PI)  *
                (
                11.0/3.0*std::log(SQR(r))*musqr
                -11.0/3.0 * (SQR(X) - SQR(Y) ) / SQR(r) * std::log( SQR(X/Y) )
                + 67.0/9.0 - SQR(M_PI)/3.0
                - 2.0 * std::log( SQR(X/r) ) * std::log( SQR(Y/r) )     ///TODO: This has nf=0
                )
            );
        } 
    }
    else if (RC_LO == BALITSKY_LO)
    {
        /*double alphas_y = Alphas(Y);
        double alphas_x = Alphas(X);
        result = 
         NC/(2.0*SQR(M_PI))*Alphas(r)
            * (
            SQR(r) / ( SQR(X) * SQR(Y) )
            + 1.0/SQR(Y)*(alphas_y/alphas_x - 1.0)
            + 1.0/SQR(X)*(alphas_x/alphas_y - 1.0)
            );*/
        cerr << "Balitsky kernel is not supported yet! " << LINEINFO << endl;
    }
    else if (RC_LO == SMALLEST_LO)
    {
        double min_size = r;
        if (X < min_size) min_size = X;
        if (Y < min_size) min_size = Y;
        result = Alphas(min_size)*NC/(2.0*M_PI*M_PI) *  r*r / (X*X * Y*Y);

        if (LO_BK)
            return result;
        
        if (EQUATION == CONFORMAL_QCD)
        {
                result *= 1.0 + Alphas(min_size)*NC / (4.0*M_PI) * (67.0/9.0 - SQR(M_PI)/3.0 - 10.0/9.0*NF/NC); 
        }
        else if (EQUATION == QCD)
        {
            result *= 1.0 + Alphas(min_size) * NC / (4.0*M_PI) * (67.0/9.0 - SQR(M_PI)/3.0 
                    - 10.0/9.0*NF/NC - 2.0 * 2.0*std::log( X/r ) * 2.0*std::log( Y/r ) 
                    ) ;
        }
    }
    else if (RC_LO == PARENT_LO)
    {

        result = Alphas(r)*NC/(2.0*M_PI*M_PI) * r*r / (X*X * Y*Y);
        if (LO_BK)
            return result;

        if (EQUATION == QCD)
            result *= 1.0 + Alphas(r) * NC / (4.0*M_PI) * (67.0/9.0 - SQR(M_PI)/3.0 
                    - 10.0/9.0*NF/NC - 2.0 * 2.0*std::log( X/r ) * 2.0*std::log( Y/r ) 
                    ) ;
        else if (EQUATION == CONFORMAL_QCD)
             result *= 1.0 + Alphas(r)*NC / (4.0*M_PI) * (67.0/9.0 - SQR(M_PI)/3.0 - 10.0/9.0*NF/NC);
        else
            cerr << "Uknwon equation! " << LINEINFO << endl;

    }
    else
    {
        cerr << "Unknown LO kernel RC! " << LINEINFO << endl;
        return -1;
    }

    if (isnan(result) or isinf(result))
    {
        cerr << "infnan " << LINEINFO << ", r=" << r << ", X=" << X << ", Y=" << Y << endl;
        exit(1);
    }
    

    return result;
}




/*
 * *****************
 * NLO PART
 * *****************
 */

double Inthelperf_nlo_z(double v, void* p);
double Inthelperf_nlo_theta_z(double theta, void* p);
double Inthelperf_nlo_z2(double v, void* p);
double Inthelperf_nlo_theta_z2(double theta, void* p);
double Inthelperf_nlo(double r, double z, double theta_z, double z2, double theta_z2, BKSolver* solver, Interpolator* dipole_interp, Interpolator* dipole_interp_s);
double Inthelperf_nlo_mc(double* vec, size_t dim, void* p);

double BKSolver::RapidityDerivative_nlo(double r, Interpolator* dipole_interp, Interpolator* dipole_interp_s)
{
    
    Inthelper_nlobk helper;
    helper.solver=this;
    helper.r=r;
    helper.dipole_interp;
    helper.dipole_interp = dipole_interp;
    helper.dipole_interp_s = dipole_interp_s;
    
    
    double minlnr = std::log( 0.5*dipole->MinR() );
    double maxlnr = std::log( 2.0*dipole->MaxR() );

    int status; double  result, abserr;

    if (INTMETHOD_NLO == MULTIPLE)
    {
        gsl_function fun;
        fun.params = &helper;
        fun.function = Inthelperf_nlo_z;
        
        gsl_integration_workspace *workspace 
        = gsl_integration_workspace_alloc(RINTPOINTS);
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
    }
    else
    {
        size_t dim=4;
        gsl_monte_function fun;
        fun.params=&helper;
        fun.f = Inthelperf_nlo_mc;
        fun.dim=dim;
        double min[4] = {minlnr, minlnr, 0, 0 };
        double max[4] = {maxlnr, maxlnr, 2.0*M_PI, 2.0*M_PI };
        const gsl_rng_type *T;
        gsl_rng *rnd;

        size_t calls = MCINTPOINTS;
        

        gsl_rng_env_setup ();

        T = gsl_rng_default;
        rnd = gsl_rng_alloc (T);

        time_t timer;
        time(&timer);
        int seconds=difftime(timer, 0);
        gsl_rng_set(rnd, seconds);

        
        
        const int maxiter_vegas=10;

        if (INTMETHOD_NLO == VEGAS)
        {
            gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (dim);
            gsl_monte_vegas_integrate (&fun, min, max, dim, calls/5, rnd, s,
                                       &result, &abserr);
            //cout <<"#Warmup result " << result << " error " << abserr << endl; 
            double prevres = result;
            int iters=0;
            do
              {
                gsl_monte_vegas_integrate (&fun, min, max, dim, calls, rnd, s,
                                           &result, &abserr);
                //#pragma omp critical
                //cout << "#Result(r=" << r <<") " << result << " err " << abserr << " relchange " << (result-prevres)/prevres << " chi^2 " << gsl_monte_vegas_chisq (s) << endl;
                prevres=result;
                iters++;
              }
              while ((std::abs( abserr/result) > 0.3 or std::abs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5 ) and iters<maxiter_vegas );
            //while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);
            //#pragma omp critical
            if (iters>=maxiter_vegas)
            {
                cout <<"# Integration failed at r=" << r <<", result->0, bestresult "<< result << " relerr " << abserr/result << " chi^2 "  << gsl_monte_vegas_chisq (s) << endl;
                result=0;
            }
            //else
            //    cout << "Integration finished, r=" << r<< ", result " << result << " relerr " << abserr/result << " chi^2 "  << gsl_monte_vegas_chisq (s) << " (intpoints " << calls << ")" << endl;
            gsl_monte_vegas_free(s);
        
        }
        else if (INTMETHOD_NLO == MISER)
        {    
        
            // plain or miser
            //gsl_monte_plain_state *s = gsl_monte_plain_alloc (4);
            gsl_monte_miser_state *s = gsl_monte_miser_alloc (4);
            int iter=0;
            
            do
            {
                iter++;
                if (iter>=2)
                {
                    cerr << "Mcintegral didn't converge in 2 iterations (r=" << r << "), result->0 " << LINEINFO << endl;
                    return 0;
                }
                //gsl_monte_plain_integrate
                gsl_monte_miser_integrate
                    (&fun, min, max, 4, calls, rnd, s,
                                       &result, &abserr);
                    //if (std::abs(abserr/result)>0.2)
                    //      cerr << "#r=" << r << " misermc integral failed, result " << result << " relerr " << std::abs(abserr/result) << ", again.... (iter " << iter << ")" << endl;
            } while (std::abs(abserr/result)>MCINTACCURACY);
            //gsl_monte_plain_free (s);
            gsl_monte_miser_free(s);
            //cout <<"#Integration finished at r=" << r <<", result " << result << " relerr " << abserr/result << " intpoints " << calls << endl;
            
            gsl_rng_free(rnd);
        }        
           
    }
    
    

    // used for arxiv version evolution for s
	//result *= -SQR(FIXED_AS*NC) / ( 16.0*M_PI*M_PI*M_PI*M_PI);   // alpha_s^2 Nc^2/(16pi^4), alphas*nc/pi=ALPHABAR_s
        // minus sign as we compute here evolution for S=1-N following ref 0710.4330
    
    return result;
}

double Inthelperf_nlo_z(double z, void* p)
{
    Inthelper_nlobk* helper = reinterpret_cast<Inthelper_nlobk*>(p);
    helper->z=std::exp(z);
    gsl_function fun;
    fun.function=Inthelperf_nlo_theta_z;
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

    result *= std::exp(2.0*z);  // Jacobian v^2 dv
    return result;
}

double Inthelperf_nlo_theta_z(double theta, void* p)
{
    Inthelper_nlobk* helper = reinterpret_cast<Inthelper_nlobk*>(p);
    helper->theta_z=theta;

    gsl_integration_workspace *workspace 
     = gsl_integration_workspace_alloc(RINTPOINTS);

    gsl_function fun;
    fun.function=Inthelperf_nlo_z2;
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

double Inthelperf_nlo_z2(double z2, void* p)
{
    Inthelper_nlobk* helper = reinterpret_cast<Inthelper_nlobk*>(p);
    helper->z2=std::exp(z2);
    gsl_function fun;
    fun.function=Inthelperf_nlo_theta_z2;
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

    result *= std::exp(2.0*z2);  // Jacobian v^2 dv
    return result;
}

double Inthelperf_nlo_theta_z2(double theta_z2, void* p)
{
    Inthelper_nlobk* helper = reinterpret_cast<Inthelper_nlobk*>(p);
    double result = Inthelperf_nlo( helper->r, helper->z, helper->theta_z, helper->z2, theta_z2, helper->solver, helper->dipole_interp, helper->dipole_interp_s);
    return result;
}

double Inthelperf_nlo(double r, double z, double theta_z, double z2, double theta_z2, BKSolver* solver, Interpolator* dipole_interp, Interpolator* dipole_interp_s)
{   
    // we choose coordinates s.t. y=0 and x lies on positive x axis
    // X = x-z = -z + r
    double X = std::sqrt(r*r + z*z - 2.0*r*z*std::cos(theta_z));  
    // Y = y-z = z
    double Y = z;
    // X' = x-z' = r - z'
    double X2=std::sqrt(r*r + z2*z2 - 2.0*r*z2*std::cos(theta_z2) );
    // Y' = y-z' = -z'
    double Y2=z2;
    // z - z'
    double z_m_z2 = std::sqrt( z*z + z2*z2 - 2.0*z*z2*std::cos(theta_z - theta_z2) );

    
    double result=0;

    if (EQUATION == QCD)
    {
        double k = solver->Kernel_nlo(r,X,Y,X2,Y2,z_m_z2);
        double kswap = solver->Kernel_nlo(r,X2,Y2,X,Y,z_m_z2);

        
        double dipole = dipole_interp->Evaluate(z_m_z2)
                            - dipole_interp->Evaluate(X)*dipole_interp->Evaluate(z_m_z2)
                            - dipole_interp->Evaluate(z_m_z2)*dipole_interp->Evaluate(Y2)
                            - dipole_interp->Evaluate(X)*dipole_interp->Evaluate(Y2)
                            + dipole_interp->Evaluate(X)*dipole_interp->Evaluate(Y)
                            + dipole_interp->Evaluate(X)*dipole_interp->Evaluate(z_m_z2)*dipole_interp->Evaluate(Y2)
                            + dipole_interp->Evaluate(Y2) - dipole_interp->Evaluate(Y); // This is not part of eq. (136) in PRD
                            
        double dipole_swap = dipole_interp->Evaluate(z_m_z2)
                            - dipole_interp->Evaluate(X2)*dipole_interp->Evaluate(z_m_z2)
                            - dipole_interp->Evaluate(z_m_z2)*dipole_interp->Evaluate(Y)
                            - dipole_interp->Evaluate(X2)*dipole_interp->Evaluate(Y)
                            + dipole_interp->Evaluate(X2)*dipole_interp->Evaluate(Y2)
                            + dipole_interp->Evaluate(X2)*dipole_interp->Evaluate(z_m_z2)*dipole_interp->Evaluate(Y)
                            + dipole_interp->Evaluate(Y) - dipole_interp->Evaluate(Y2);
        
        //result = k*dipole;
        result = (k*dipole + kswap*dipole_swap)/2.0;

        if (NF>0)
        {
            double kernel_f = solver->Kernel_nlo_fermion(r,X,Y,X2,Y2,z_m_z2);
            double kernel_f_swap = solver->Kernel_nlo_fermion(r,X2,Y2,X,Y,z_m_z2);

            //double dipole_f = dipole_interp_s->Evaluate(Y) * ( dipole_interp_s->Evaluate(X2) - dipole_interp_s->Evaluate(X) );
            double dipole_f = dipole_interp->Evaluate(X) - dipole_interp->Evaluate(X2)
                - dipole_interp->Evaluate(X)*dipole_interp->Evaluate(Y) + dipole_interp->Evaluate(X2)*dipole_interp->Evaluate(Y);
            //double dipole_f_swap = dipole_interp_s->Evaluate(Y2) * ( dipole_interp_s->Evaluate(X) - dipole_interp_s->Evaluate(X2) );
            double dipole_f_swap = dipole_interp->Evaluate(X2) - dipole_interp->Evaluate(X)
                - dipole_interp->Evaluate(X2)*dipole_interp->Evaluate(Y2) + dipole_interp->Evaluate(X)*dipole_interp->Evaluate(Y2);

            result += -(kernel_f*dipole_f + kernel_f_swap * dipole_f_swap)/2.0;     // Minus sign as the evolution is written for S and we solve N

        }
    }

    // Evolution for conformal dipole
    else if (EQUATION==CONFORMAL_QCD)
    {
        double k1 = solver->Kernel_nlo_conformal_1(r,X,Y,X2,Y2,z_m_z2);        
        double k2 = solver->Kernel_nlo_conformal_2(r,X,Y,X2,Y2,z_m_z2);

        double dipole1 = dipole_interp_s->Evaluate(X) * dipole_interp_s->Evaluate(z_m_z2) * dipole_interp_s->Evaluate(Y2)
            - dipole_interp_s->Evaluate(X) * dipole_interp_s->Evaluate(Y);
        double dipole2 = dipole_interp_s->Evaluate(X) * dipole_interp_s->Evaluate(z_m_z2) * dipole_interp_s->Evaluate(Y2)
            - dipole_interp_s->Evaluate(X2) * dipole_interp_s->Evaluate(z_m_z2) * dipole_interp_s->Evaluate(Y);

        double k1_swap = solver->Kernel_nlo_conformal_1(r,X2,Y2,X,Y,z_m_z2);        
        double k2_swap = solver->Kernel_nlo_conformal_2(r,X2,Y2,X,Y,z_m_z2);

        double dipole1_swap = dipole_interp_s->Evaluate(X2) * dipole_interp_s->Evaluate(z_m_z2) * dipole_interp_s->Evaluate(Y)
            - dipole_interp_s->Evaluate(X2) * dipole_interp_s->Evaluate(Y2);
        double dipole2_swap = dipole_interp_s->Evaluate(X2) * dipole_interp_s->Evaluate(z_m_z2) * dipole_interp_s->Evaluate(Y)
            - dipole_interp_s->Evaluate(X) * dipole_interp_s->Evaluate(z_m_z2) * dipole_interp_s->Evaluate(Y2);

        result = (k1*dipole1 + k2*dipole2 + k1_swap*dipole1_swap + k2_swap*dipole2_swap)/2.0;

        /// Fermion part
        if (NF > 0)
        {
            double kernel_f = solver->Kernel_nlo_conformal_fermion(r,X,Y,X2,Y2,z_m_z2);
            double dipole_f = dipole_interp_s->Evaluate(Y) * ( dipole_interp_s->Evaluate(X2)
                                                                - dipole_interp_s->Evaluate(X) );

            double kernel_f_swap = solver->Kernel_nlo_conformal_fermion(r,X2,Y2,X,Y,z_m_z2);
            double dipole_f_swap = dipole_interp_s->Evaluate(Y2) * ( dipole_interp_s->Evaluate(X)
                                                                - dipole_interp_s->Evaluate(X2) );

            result += (kernel_f * dipole_f + kernel_f_swap * dipole_f_swap)/2.0;
        }

        result *= -1.0; // Minus sign as the evolution is written for S but we solve N = 1-S

    }

    else if (EQUATION == CONFORMAL_N4)
    {
        result = dipole_interp_s->Evaluate(X) * dipole_interp_s->Evaluate(z_m_z2) * dipole_interp_s->Evaluate(Y2)
            - dipole_interp_s->Evaluate(X) * dipole_interp_s->Evaluate(Y);

        result *= solver->Kernel_nlo_n4_sym(r,X,Y,X2,Y2,z_m_z2);

        result *= -1.0; // Minus sign as the evolution is written for S but we solve N = 1-S
    }

    else
    {
        cerr << "Unknown equation to solve: " << EQUATION << endl;
    }
    


    /////// Coefficients

    // If the alphas scale is set by the smallest dipole, multiply the kernel here by as^2 and
    // other relevant factors
    if (EQUATION == CONFORMAL_N4)
        result *= SQR(FIXED_AS*NC)/(8.0*std::pow(M_PI, 4));
        
    else if (EQUATION == QCD or EQUATION == CONFORMAL_QCD)
    {
        if (RC_NLO == FIXED_NLO)
            result *= SQR(FIXED_AS*NC) / (8.0*std::pow(M_PI,4) );
        else if (RC_NLO == PARENT_NLO)
            result *= SQR( solver->Alphas(r) * NC) / (8.0 * std::pow(M_PI, 4) );
        else if (RC_NLO == SMALLEST_NLO)
        {
            double min_size = r;
            if (X < min_size) min_size = X;
            if (Y < min_size) min_size = Y;
            if (X2 < min_size) min_size=X2;
            if (Y2 < min_size) min_size = Y2;
            if (z_m_z2 < min_size) min_size = z_m_z2;

            result *= SQR( solver->Alphas(min_size) * NC) / (8.0 * std::pow(M_PI, 4) );   
        }
        else 
        {
            cerr << "Unknown NLO kernel alphas! " << LINEINFO << endl;
            return -1;
        }


    }
    else
        cerr << "WTF! " << LINEINFO << endl;
    
    
    
    


  


	if (isnan(result) or isinf(result))
	{
        return 0;    
    }

    return result;

}

/*
 * GSL wrapper for monte carlo integration
 * vec is: [ln u, ln v, theta_u, theta_v]
 */
double Inthelperf_nlo_mc(double* vec, size_t dim, void* p)
{
    Inthelper_nlobk* helper = reinterpret_cast<Inthelper_nlobk*>(p);

    if (dim != 4)
    {
        cerr << "Insane dimension at " << LINEINFO << ": " << dim << endl;
        exit(1);
    }

    double integrand = Inthelperf_nlo(helper->r, std::exp(vec[0]), vec[2], std::exp(vec[1]), vec[3], helper->solver, helper->dipole_interp, helper->dipole_interp_s)
        * std::exp(2.0*vec[0]) * std::exp(2.0*vec[1]);  // Jacobian

    //cout << std::exp(vec[0]) << " " << std::exp(vec[1]) << " " << integrand << endl;

    return integrand; 
    
}

/***************************************************
* Non-conformal kernels
* Note that as^2 nc^2 / (8pi^4) is taken out from the kernels
****************************************************/


/*
 * NLO evolution kernel for non-conformal N
 */
double BKSolver::Kernel_nlo(double r, double X, double Y, double X2, double Y2, double z_m_z2)
{

    double kernel = -2.0/std::pow(z_m_z2,4);

    kernel += (
        ( SQR(X*Y2) + SQR(X2*Y) - 4.0*SQR(r*z_m_z2) ) / ( std::pow(z_m_z2,4) * (SQR(X*Y2) - SQR(X2*Y)) )
        + std::pow(r,4) / ( SQR(X*Y2)*( SQR(X*Y2) - SQR(X2*Y) ) )
        + SQR(r) / ( SQR(X*Y2*z_m_z2) )
        ) * std::log( SQR(X*Y2/(X2*Y)) );

    if (isnan(kernel) or isinf(kernel))
    {
        //cerr << "Kernel " << kernel <<", r=" << r <<", X=" << X << ", Y=" << Y <<", X2=" << X2 <<", Y2=" << Y2 <<", z-z2=" << z_m_z2 << endl;
        return 0;
    }
    
    return kernel;
}

double BKSolver::Kernel_nlo_fermion(double r, double X, double Y, double X2, double Y2, double z_m_z2)
{
    double kernel=0;

    kernel = 2.0 / std::pow(z_m_z2,4.0);

    kernel -= (SQR(X*Y2) + SQR(X2*Y) - SQR(r*z_m_z2) ) / ( std::pow(z_m_z2, 4.0)*(SQR(X*Y2) - SQR(X2*Y)) )
                * 2.0*std::log(X*Y2/(X2*Y));

    kernel *= NF/NC;        // Divided by NC, as in the kernel we have as^2 nc nf/(8pi^4), but
                    // this kernel is multiplied yb as^2 nc^2/(8pi^4)

    if (isinf(kernel) or isnan(kernel))
        return 0;
        

    return kernel;
}

/***************************************************
* Conformal kernels
* Note that as^2 nc^2 / (8pi^2) is taken out from the kernels
****************************************************/
double BKSolver::Kernel_nlo_conformal_1(double r, double X, double Y, double X2, double Y2, double z_m_z2)
{
    double result = 0;


	
    result = ( SQR(r/z_m_z2) * (1.0/SQR(X*Y2) - 1.0/SQR(X2*Y) )
                + std::pow(r,4) / ( SQR(X*Y2) - SQR(X2*Y) ) * ( 1.0/SQR(X*Y2) + 1.0/SQR(X2*Y) )
                    + 2.0*(SQR(X*Y2) + SQR(X2*Y) - 4.0*SQR(r*z_m_z2) )/( std::pow(z_m_z2,4) * (SQR(X*Y2) - SQR(X2*Y) ) )    )
                * 2.0*std::log(X*Y2/(X2*Y)) ;
                
    
    result += -4.0/std::pow(z_m_z2,4) + 2.0*SQR(r/(z_m_z2*X*Y2)) * 2.0*std::log(r*z_m_z2/(X2*Y))
                                      + 2.0*SQR(r/(z_m_z2*X2*Y)) * 2.0*std::log(r*z_m_z2/(X*Y2));
    
    
    if (isinf(result) or isnan(result))
        return 0;

    result /= 2.0;      // Conformal kernel is multiplied by as^2 Nc^2/(8pi^4) in the Inthelperf_nlo function,
                        // but as the gluon part has prefactor as^2 Nc^2/(8pi^4), conformal kernels 1 and 2 are
                        // divided by 2 here
    return result;
}


double BKSolver::Kernel_nlo_conformal_2(double r, double X, double Y, double X2, double Y2, double z_m_z2)
{
    double result = 0;
    result = (SQR(r/(z_m_z2*X*Y2)) + std::pow(r,4)/(SQR(X*Y2)*(SQR(X*Y2) - SQR(X2*Y) ) )  ) * std::log(SQR(X*Y2/(X2*Y)));
    result += 2.0*SQR(r/(z_m_z2*X*Y2)) * 2.0*std::log( r*z_m_z2 / (X2*Y) );

    result /= 2.0;      // Conformal kernel is multiplied by as^2 Nc^2/(8pi^4) in the Inthelperf_nlo function,
                        // but as the gluon part has prefactor as^2 Nc^2/(8pi^4), conformal kernels 1 and 2 are
                        // divided by 2 here
    return result;
}

double BKSolver::Kernel_nlo_conformal_fermion(double r, double X, double Y, double X2, double Y2, double z_m_z2)
{
    double result=0;
    result = 2.0 - ( SQR(X*Y2) + SQR(X2*Y) - SQR(r*z_m_z2) ) / ( SQR(X*Y2) - SQR(X2*Y) )
                    * std::log( SQR( X*Y2/(X2*Y) ) );

    result *= 1.0 / std::pow(z_m_z2, 4);

    result *= NF/NC;        // Divided by NC, as in the kernel we have as^2 nc nf/(8pi^4), but
                    // this kernel is multiplied yb as^2 nc^2/(8pi^4)

    return result;
}



/***************************************************
 * N=4 SYM kernel
 **************************************************/
double BKSolver::Kernel_nlo_n4_sym(double r, double X, double Y, double X2, double Y2, double z_m_z2)
{
    double result=0;

    result = 2.0 * 2.0*std::log( r*z_m_z2/(X2*Y) )
        + (1.0 + SQR(r*z_m_z2) / ( SQR(X*Y2) - SQR(X2*Y) )) * 2.0*std::log(X*Y2/(X2*Y) );
    result *= SQR(r/(X*Y2*z_m_z2));

    return result;

}

////////////////////////////////////////////////
// Running coupling
double BKSolver::Alphas(double r)
{
	double Csqr=1;
	double scalefactor = 4.0*Csqr;
	double rsqr = r*r;
	double maxalphas=0.7;
	double lambdaqcd=LAMBDAQCD;
        
	if (scalefactor/(rsqr*lambdaqcd*lambdaqcd) < 1.0) return maxalphas;
	double alpha = 12.0*M_PI/( (11.0*NC-2.0*NF)*std::log(scalefactor/ (rsqr*lambdaqcd*lambdaqcd) ) );
	if (alpha>maxalphas)
		return maxalphas;
	return alpha;
}


Dipole* BKSolver::GetDipole()
{
    return dipole;
}

void BKSolver::SetTmpOutput(std::string fname)
{
    tmp_output = fname;
}
